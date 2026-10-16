//! Paired gzipped FASTQ reader for split-pool reads, exposed through a C ABI.
//!
//! One background thread per mate decompresses and parses its file, so both
//! files are read ahead while the caller aligns the previous batch. Read 1
//! (cDNA) is returned as 2-bit codes (A, C, G, T = 0..3, anything else 4, as
//! `tealeaf.sc.native_align.CODE`). Read 2 is reduced to the library cell index
//! and the UMI: each barcode segment is looked up in its correction table
//! (`segment_table`, indexed by the segment's 2-bit code with N read as T, as
//! `pack`), the three corrected segments index the cell lookup, and the UMI
//! is valid when it has no N and is not a homopolymer.

use flate2::read::MultiGzDecoder;
use std::ffi::{CStr, CString};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::os::raw::c_char;
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::thread;

const QUEUE: usize = 2;

fn code(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

struct Read1 {
    codes: Vec<u8>,
    n: usize,
    length: usize,
}

struct Read2 {
    cell: Vec<i64>,
    umi: Vec<i64>,
    valid: Vec<u8>,
}

#[derive(Clone)]
struct Barcodes {
    tables: Vec<i64>,
    segment_length: usize,
    segments: [(usize, usize); 3],
    lookup: Vec<i64>,
    size: i64,
    umi: (usize, usize),
}

impl Barcodes {
    fn pack(sequence: &[u8]) -> (i64, bool) {
        let mut value = 0i64;
        let mut has_n = false;
        for &base in sequence {
            let c = code(base);
            has_n |= c > 3;
            value = value * 4 + c.min(3) as i64;
        }
        (value, has_n)
    }

    fn cell(&self, sequence: &[u8]) -> i64 {
        let mut index = 0i64;
        for (i, &(a, b)) in self.segments.iter().enumerate() {
            let (value, _) = Self::pack(&sequence[a..b]);
            let corrected = self.tables[(i << (2 * self.segment_length)) + value as usize];
            if corrected < 0 {
                return -1;
            }
            index = index * self.size + corrected;
        }
        self.lookup[index as usize]
    }

    fn umi(&self, sequence: &[u8]) -> (i64, bool) {
        let (value, has_n) = Self::pack(&sequence[self.umi.0..self.umi.1]);
        let length = (self.umi.1 - self.umi.0) as u32;
        let repeat = (4i64.pow(length) - 1) / 3;
        let homopolymer = value % repeat == 0;
        (value, !has_n && !homopolymer)
    }

    fn needed(&self) -> usize {
        self.segments.iter().map(|s| s.1).max().unwrap().max(self.umi.1)
    }
}

/// Calls `record` with the sequence line of every FASTQ record of a gzipped file.
fn records(path: &str, mut record: impl FnMut(&[u8]) -> Result<bool, String>) -> Result<(), String> {
    let file = File::open(path).map_err(|e| format!("{path}: {e}"))?;
    let mut reader = BufReader::with_capacity(1 << 22, MultiGzDecoder::new(BufReader::with_capacity(1 << 22, file)));
    let mut lines: [Vec<u8>; 4] = Default::default();
    loop {
        for (i, line) in lines.iter_mut().enumerate() {
            line.clear();
            let n = reader.read_until(b'\n', line).map_err(|e| format!("{path}: {e}"))?;
            if n == 0 {
                return if i == 0 { Ok(()) } else { Err(format!("{path}: truncated record")) };
            }
        }
        let sequence = lines[1].strip_suffix(b"\n").unwrap_or(&lines[1]);
        let sequence = sequence.strip_suffix(b"\r").unwrap_or(sequence);
        if !record(sequence)? {
            return Ok(());
        }
    }
}

fn read1_thread(path: String, batch: usize, sender: SyncSender<Result<Read1, String>>) {
    let mut current = Read1 { codes: Vec::new(), n: 0, length: 0 };
    let result = records(&path, |sequence| {
        if current.length == 0 {
            current.length = sequence.len();
            current.codes.reserve(batch * sequence.len());
        } else if sequence.len() != current.length {
            return Err(format!("{path}: variable read lengths are not supported"));
        }
        current.codes.extend(sequence.iter().map(|&b| code(b)));
        current.n += 1;
        if current.n == batch {
            let length = current.length;
            let full = std::mem::replace(&mut current, Read1 { codes: Vec::with_capacity(batch * length), n: 0, length });
            return Ok(sender.send(Ok(full)).is_ok());
        }
        Ok(true)
    });
    let _ = match result {
        Ok(()) if current.n > 0 => sender.send(Ok(current)),
        Ok(()) => Ok(()),
        Err(message) => sender.send(Err(message)),
    };
}

fn read2_thread(path: String, batch: usize, barcodes: Barcodes, sender: SyncSender<Result<Read2, String>>) {
    let needed = barcodes.needed();
    let empty = || Read2 { cell: Vec::with_capacity(batch), umi: Vec::with_capacity(batch), valid: Vec::with_capacity(batch) };
    let mut current = empty();
    let result = records(&path, |sequence| {
        if sequence.len() < needed {
            return Err(format!("{path}: read 2 shorter than the barcode layout ({needed} bases)"));
        }
        let (umi, valid) = barcodes.umi(sequence);
        current.cell.push(barcodes.cell(sequence));
        current.umi.push(umi);
        current.valid.push(valid as u8);
        if current.cell.len() == batch {
            return Ok(sender.send(Ok(std::mem::replace(&mut current, empty()))).is_ok());
        }
        Ok(true)
    });
    let _ = match result {
        Ok(()) if !current.cell.is_empty() => sender.send(Ok(current)),
        Ok(()) => Ok(()),
        Err(message) => sender.send(Err(message)),
    };
}

pub struct Pair {
    read1: Receiver<Result<Read1, String>>,
    read2: Receiver<Result<Read2, String>>,
    pending: Option<(Read1, Read2)>,
    error: CString,
}

impl Pair {
    fn next(&mut self) -> Result<Option<(Read1, Read2)>, String> {
        let first = self.read1.recv().ok().transpose()?;
        let second = self.read2.recv().ok().transpose()?;
        match (first, second) {
            (None, None) => Ok(None),
            (Some(a), Some(b)) if a.n == b.cell.len() => Ok(Some((a, b))),
            _ => Err("read 1 and read 2 have different record counts".to_string()),
        }
    }
}

unsafe fn string(pointer: *const c_char) -> String {
    CStr::from_ptr(pointer).to_string_lossy().into_owned()
}

/// Open a pair of gzipped FASTQ files and start reading them in batches of
/// `batch` records. `tables` holds the three segment correction tables
/// (each 4^`segment_length` entries), `segment_starts`/`segment_ends` the
/// segment positions in read 2, `lookup` the cell index of each corrected
/// segment triple (`size`^3 entries, -1 outside the library) and
/// [`umi_start`, `umi_end`) the UMI.
///
/// # Safety
/// Paths must be NUL-terminated strings and arrays valid for their lengths.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn fastq_pair_open(
    path1: *const c_char, path2: *const c_char, batch: usize,
    tables: *const i64, segment_length: usize, segment_starts: *const usize, segment_ends: *const usize,
    lookup: *const i64, size: i64, umi_start: usize, umi_end: usize,
) -> *mut Pair {
    let starts = std::slice::from_raw_parts(segment_starts, 3);
    let ends = std::slice::from_raw_parts(segment_ends, 3);
    let barcodes = Barcodes {
        tables: std::slice::from_raw_parts(tables, 3 << (2 * segment_length)).to_vec(),
        segment_length,
        segments: [(starts[0], ends[0]), (starts[1], ends[1]), (starts[2], ends[2])],
        lookup: std::slice::from_raw_parts(lookup, (size * size * size) as usize).to_vec(),
        size,
        umi: (umi_start, umi_end),
    };
    let (send1, read1) = sync_channel(QUEUE);
    let (send2, read2) = sync_channel(QUEUE);
    let (path1, path2) = (string(path1), string(path2));
    thread::spawn(move || read1_thread(path1, batch, send1));
    thread::spawn(move || read2_thread(path2, batch, barcodes, send2));
    Box::into_raw(Box::new(Pair { read1, read2, pending: None, error: CString::default() }))
}

/// Wait for the next batch; returns its number of records (0 at the end of
/// both files, -1 on error, see `fastq_pair_error`) and writes the read-1
/// length to `length`. Copy the batch out with `fastq_pair_take`.
///
/// # Safety
/// `pair` must come from `fastq_pair_open`.
#[no_mangle]
pub unsafe extern "C" fn fastq_pair_peek(pair: *mut Pair, length: *mut usize) -> i64 {
    let pair = &mut *pair;
    match pair.next() {
        Ok(Some((a, b))) => {
            let n = a.n as i64;
            *length = a.length;
            pair.pending = Some((a, b));
            n
        }
        Ok(None) => 0,
        Err(message) => {
            pair.error = CString::new(message).unwrap_or_default();
            -1
        }
    }
}

/// Copy the batch returned by the last `fastq_pair_peek`: read-1 codes
/// (n x length), cell index (-1 when not a library cell), UMI code and UMI
/// validity.
///
/// # Safety
/// Output arrays must hold the batch.
#[no_mangle]
pub unsafe extern "C" fn fastq_pair_take(pair: *mut Pair, codes: *mut u8, cell: *mut i64, umi: *mut i64, valid: *mut u8) {
    let pair = &mut *pair;
    if let Some((a, b)) = pair.pending.take() {
        std::ptr::copy_nonoverlapping(a.codes.as_ptr(), codes, a.codes.len());
        std::ptr::copy_nonoverlapping(b.cell.as_ptr(), cell, b.cell.len());
        std::ptr::copy_nonoverlapping(b.umi.as_ptr(), umi, b.umi.len());
        std::ptr::copy_nonoverlapping(b.valid.as_ptr(), valid, b.valid.len());
    }
}

/// Message of the last error.
///
/// # Safety
/// `pair` must come from `fastq_pair_open`.
#[no_mangle]
pub unsafe extern "C" fn fastq_pair_error(pair: *mut Pair) -> *const c_char {
    (*pair).error.as_ptr()
}

/// Stop reading and free the reader.
///
/// # Safety
/// `pair` must come from `fastq_pair_open` and not be used afterwards.
#[no_mangle]
pub unsafe extern "C" fn fastq_pair_close(pair: *mut Pair) {
    drop(Box::from_raw(pair));
}
