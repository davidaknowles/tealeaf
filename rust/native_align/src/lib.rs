//! Rust port of the read alignment kernel of `tealeaf.sc.native_align`
//! (`_align_reads` and its helpers), exposed as `align_reads` through a C ABI.
//! The algorithm, constants, scoring and tie-breaking follow the numba
//! implementation line for line so that both give identical output; see the
//! Python module for the description of the method.

use std::sync::atomic::{AtomicUsize, Ordering};

const MIN_INTRON: i64 = 21;
const MAX_INTRON: i64 = 1_000_000;
const MAX_INDEL: i64 = 6;
const MAX_OCCURRENCES: usize = 32;
const MAX_HITS: usize = 1024;
const MAX_SEGMENTS: usize = 128;
const MAX_BLOCKS: usize = 8;
const MAX_CANDIDATES: usize = 24;
const STRIDE: usize = 4;
const BUCKET_BITS: u64 = 24;

const UNIQUE: i8 = 1;
const MULTI: i8 = 2;
const NOT_GENOME_UNIQUE: i8 = 3;
const TOO_SHORT: i8 = 4;

struct Reference<'a> {
    seq: &'a [u8],
    index_codes: &'a [u64],
    index_positions: &'a [i64],
    buckets: &'a [i64],
    unique: &'a [u8],
    k: usize,
    shift: u64,
    locus_start: &'a [i64],
    locus_end: &'a [i64],
    tx_space_start: i64,
    tx_ref_start: &'a [i64],
    tx_length: &'a [i64],
    exon_offsets: &'a [i64],
    exon_tx_start: &'a [i64],
    exon_ref_start: &'a [i64],
    exon_length: &'a [i64],
    junctions: &'a [i64],
    donors: &'a [i64],
    acceptors: &'a [i64],
}

impl<'a> Reference<'a> {
    #[inline(always)]
    fn base(&self, position: i64) -> u8 {
        self.seq[position as usize]
    }
}

struct Candidates {
    n: usize,
    score: [i64; MAX_CANDIDATES],
    matches: [i64; MAX_CANDIDATES],
    mm: [i64; MAX_CANDIDATES],
    start: [i64; MAX_CANDIDATES],
    end: [i64; MAX_CANDIDATES],
    blocks: [usize; MAX_CANDIDATES],
    rs: [[i64; MAX_BLOCKS]; MAX_CANDIDATES],
    re: [[i64; MAX_BLOCKS]; MAX_CANDIDATES],
    dg: [[i64; MAX_BLOCKS]; MAX_CANDIDATES],
    unique: [bool; MAX_CANDIDATES],
    orient: [i8; MAX_CANDIDATES],
}

impl Candidates {
    fn new() -> Self {
        Candidates {
            n: 0,
            score: [0; MAX_CANDIDATES],
            matches: [0; MAX_CANDIDATES],
            mm: [0; MAX_CANDIDATES],
            start: [0; MAX_CANDIDATES],
            end: [0; MAX_CANDIDATES],
            blocks: [0; MAX_CANDIDATES],
            rs: [[0; MAX_BLOCKS]; MAX_CANDIDATES],
            re: [[0; MAX_BLOCKS]; MAX_CANDIDATES],
            dg: [[0; MAX_BLOCKS]; MAX_CANDIDATES],
            unique: [false; MAX_CANDIDATES],
            orient: [0; MAX_CANDIDATES],
        }
    }
}

struct Work {
    hits_diag: Vec<i64>,
    hits_read: Vec<i64>,
    hit_order: Vec<usize>,
    seg_diag: Vec<i64>,
    seg_start: Vec<i64>,
    seg_end: Vec<i64>,
    rs: [i64; MAX_BLOCKS],
    re: [i64; MAX_BLOCKS],
    dg: [i64; MAX_BLOCKS],
}

impl Work {
    fn new() -> Self {
        Work {
            hits_diag: Vec::with_capacity(MAX_HITS),
            hits_read: Vec::with_capacity(MAX_HITS),
            hit_order: Vec::with_capacity(MAX_HITS),
            seg_diag: Vec::with_capacity(MAX_SEGMENTS),
            seg_start: Vec::with_capacity(MAX_SEGMENTS),
            seg_end: Vec::with_capacity(MAX_SEGMENTS),
            rs: [0; MAX_BLOCKS],
            re: [0; MAX_BLOCKS],
            dg: [0; MAX_BLOCKS],
        }
    }
}

fn lookup(code: u64, r: &Reference) -> (usize, usize) {
    let b = (code >> r.shift) as usize;
    let (mut low, mut high) = (r.buckets[b] as usize, r.buckets[b + 1] as usize);
    let end = high;
    while low < high {
        let middle = (low + high) / 2;
        if r.index_codes[middle] < code { low = middle + 1 } else { high = middle }
    }
    let first = low;
    high = end;
    while low < high {
        let middle = (low + high) / 2;
        if r.index_codes[middle] <= code { low = middle + 1 } else { high = middle }
    }
    (first, low)
}

/// Index of the last element <= value, or -1.
fn find(sorted: &[i64], value: i64) -> i64 {
    let (mut low, mut high) = (0usize, sorted.len());
    while low < high {
        let middle = (low + high) / 2;
        if sorted[middle] <= value { low = middle + 1 } else { high = middle }
    }
    low as i64 - 1
}

fn contains(sorted: &[i64], value: i64) -> bool {
    let i = find(sorted, value);
    i >= 0 && sorted[i as usize] == value
}

fn motif(r: &Reference, donor: i64, acceptor: i64) -> i64 {
    let (a, b, c, d) = (r.base(donor), r.base(donor + 1), r.base(acceptor - 2), r.base(acceptor - 1));
    match (a, b, c, d) {
        (2, 3, 0, 2) => 1,
        (1, 3, 0, 1) => 2,
        (2, 1, 0, 2) => 3,
        (1, 3, 2, 1) => 4,
        (0, 3, 0, 1) => 5,
        (2, 3, 0, 3) => 6,
        _ => 0,
    }
}

fn junction_score(motif: i64, annotated: bool) -> i64 {
    if annotated {
        2
    } else if motif == 1 || motif == 2 {
        0
    } else if motif == 3 || motif == 4 {
        -4
    } else {
        -8
    }
}

#[inline(always)]
fn mismatch(q: &[u8], r: i64, rf: &Reference, diagonal: i64) -> i64 {
    let b = q[r as usize];
    if b > 3 { 0 } else if rf.base(diagonal + r) != b { 1 } else { 0 }
}

fn extend_left(q: &[u8], rf: &Reference, diagonal: i64, start: i64, floor: i64) -> i64 {
    let (mut best, mut score, mut best_start) = (0i64, 0i64, start);
    let mut r = start - 1;
    while r >= floor {
        if rf.base(diagonal + r) > 3 {
            break;
        }
        score += if mismatch(q, r, rf, diagonal) == 1 { -1 } else if q[r as usize] <= 3 { 1 } else { 0 };
        if score > best {
            best = score;
            best_start = r;
        }
        if score < best - 8 {
            break;
        }
        r -= 1;
    }
    best_start
}

fn extend_right(q: &[u8], rf: &Reference, diagonal: i64, end: i64, ceiling: i64) -> i64 {
    let (mut best, mut score, mut best_end) = (0i64, 0i64, end);
    let mut r = end;
    while r < ceiling {
        if rf.base(diagonal + r) > 3 {
            break;
        }
        score += if mismatch(q, r, rf, diagonal) == 1 { -1 } else if q[r as usize] <= 3 { 1 } else { 0 };
        if score > best {
            best = score;
            best_end = r + 1;
        }
        if score < best - 8 {
            break;
        }
        r += 1;
    }
    best_end
}

fn score_blocks(q: &[u8], rf: &Reference, rs: &[i64], re: &[i64], dg: &[i64], n: usize) -> (i64, i64) {
    let (mut matches, mut mismatches) = (0i64, 0i64);
    for i in 0..n {
        for r in rs[i]..re[i] {
            let b = q[r as usize];
            if b > 3 {
                continue;
            }
            if rf.base(dg[i] + r) == b { matches += 1 } else { mismatches += 1 }
        }
    }
    (matches, mismatches)
}

fn genome_unique(q: &[u8], rf: &Reference, rs: &[i64], re: &[i64], dg: &[i64], n: usize) -> bool {
    let k = rf.k as i64;
    for i in 0..n {
        let mut run = 0i64;
        for r in rs[i]..re[i] {
            let b = q[r as usize];
            if b <= 3 && rf.base(dg[i] + r) == b { run += 1 } else { run = 0 }
            if run >= k && rf.unique[(dg[i] + r - k + 1) as usize] != 0 {
                return true;
            }
        }
    }
    false
}

#[allow(clippy::too_many_arguments)]
fn transcript_candidate(q: &[u8], rf: &Reference, diagonal: i64, seed_start: i64, seed_end: i64, tx_start: i64, tx_end: i64, t: usize, w: &mut Work) -> usize {
    let length = q.len() as i64;
    let floor = 0i64.max(tx_start - diagonal);
    let ceiling = length.min(tx_end - diagonal);
    let start = extend_left(q, rf, diagonal, seed_start, floor);
    let end = extend_right(q, rf, diagonal, seed_end, ceiling);
    let mut n = 0usize;
    for e in rf.exon_offsets[t] as usize..rf.exon_offsets[t + 1] as usize {
        let e_start = tx_start + rf.exon_tx_start[e];
        let e_end = e_start + rf.exon_length[e];
        let low = (diagonal + start).max(e_start);
        let high = (diagonal + end).min(e_end);
        if low >= high {
            continue;
        }
        if n >= MAX_BLOCKS {
            return 0;
        }
        w.rs[n] = low - diagonal;
        w.re[n] = high - diagonal;
        w.dg[n] = rf.exon_ref_start[e] - e_start + diagonal;
        n += 1;
    }
    n
}

/// Best split s in [low, high] joining diagonals a and b: (s, motif, annotated), or None.
#[allow(clippy::too_many_arguments)]
fn join(q: &[u8], rf: &Reference, a_diagonal: i64, b_diagonal: i64, low: i64, high: i64, locus_start: i64, locus_end: i64) -> Option<(i64, i64, bool)> {
    let (mut best_s, mut best_mm, mut best_motif, mut best_annotated, mut best_rank) = (-1i64, 1i64 << 30, 0i64, false, -1i64);
    let gap = b_diagonal - a_diagonal;
    for s in low..=high {
        let mut mm = 0;
        for r in low..s {
            mm += mismatch(q, r, rf, a_diagonal);
        }
        for r in s..high {
            mm += mismatch(q, r, rf, b_diagonal);
        }
        let (mut m, mut annotated, mut rank) = (0i64, false, 0i64);
        if gap >= MIN_INTRON {
            let donor = a_diagonal + s;
            let acceptor = b_diagonal + s;
            if donor < locus_start || acceptor > locus_end {
                continue;
            }
            annotated = contains(rf.junctions, (donor << 32) | acceptor);
            m = motif(rf, donor, acceptor);
            if !annotated && m == 0 {
                continue;
            }
            rank = if annotated { 3 } else if m <= 2 { 2 } else { 1 };
        }
        if mm < best_mm || (mm == best_mm && rank > best_rank) {
            best_s = s;
            best_mm = mm;
            best_motif = m;
            best_annotated = annotated;
            best_rank = rank;
        }
    }
    if best_s >= 0 { Some((best_s, best_motif, best_annotated)) } else { None }
}

fn rescue_right(q: &[u8], rf: &Reference, diagonal: i64, block_start: i64, end: i64, locus_end: i64) -> (i64, i64, i64) {
    let length = q.len() as i64;
    let (mut best_split, mut best_acceptor, mut best_delta) = (-1i64, -1i64, 0i64);
    if length - end < 5 {
        return (best_split, best_acceptor, best_delta);
    }
    let first = (find(rf.acceptors, diagonal + end + MIN_INTRON - 1) + 1) as usize;
    for &acceptor in &rf.acceptors[first..] {
        if acceptor > diagonal + end + MAX_INTRON || acceptor + (length - end) > locus_end {
            break;
        }
        for split in (block_start + 1).max(end - 6)..=end {
            let donor = diagonal + split;
            let annotated = contains(rf.junctions, (donor << 32) | acceptor);
            let m = motif(rf, donor, acceptor);
            if !annotated && m == 0 {
                continue;
            }
            let (mut gain, mut lost, mut tail_mm) = (0i64, 0i64, 0i64);
            for r in split..length {
                let b = q[r as usize];
                if b > 3 {
                    continue;
                }
                if rf.base(acceptor + r - split) == b {
                    gain += 1
                } else {
                    gain -= 1;
                    tail_mm += 1
                }
            }
            for r in split..end {
                let b = q[r as usize];
                if b <= 3 {
                    lost += if rf.base(diagonal + r) == b { 1 } else { -1 };
                }
            }
            if tail_mm > 1i64.max((length - split) / 10) {
                continue;
            }
            let delta = gain - lost + junction_score(m, annotated);
            if delta > best_delta {
                best_split = split;
                best_acceptor = acceptor;
                best_delta = delta;
            }
        }
    }
    (best_split, best_acceptor, best_delta)
}

fn rescue_left(q: &[u8], rf: &Reference, diagonal: i64, start: i64, block_end: i64, locus_start: i64) -> (i64, i64, i64) {
    let (mut best_split, mut best_donor, mut best_delta) = (-1i64, -1i64, 0i64);
    if start < 5 {
        return (best_split, best_donor, best_delta);
    }
    let mut i = find(rf.donors, diagonal + start - MIN_INTRON + 1);
    while i >= 0 {
        let donor = rf.donors[i as usize];
        if donor < diagonal + start - MAX_INTRON || donor - start < locus_start {
            break;
        }
        for split in start..=(start + 6).min(block_end - 1) {
            let acceptor = diagonal + split;
            let annotated = contains(rf.junctions, (donor << 32) | acceptor);
            let m = motif(rf, donor, acceptor);
            if !annotated && m == 0 {
                continue;
            }
            let (mut gain, mut lost, mut head_mm) = (0i64, 0i64, 0i64);
            for r in 0..split {
                let b = q[r as usize];
                if b > 3 {
                    continue;
                }
                if rf.base(donor - split + r) == b {
                    gain += 1
                } else {
                    gain -= 1;
                    head_mm += 1
                }
            }
            for r in start..split {
                let b = q[r as usize];
                if b <= 3 {
                    lost += if rf.base(diagonal + r) == b { 1 } else { -1 };
                }
            }
            if head_mm > 1i64.max(split / 10) {
                continue;
            }
            let delta = gain - lost + junction_score(m, annotated);
            if delta > best_delta {
                best_split = split;
                best_donor = donor;
                best_delta = delta;
            }
        }
        i -= 1;
    }
    (best_split, best_donor, best_delta)
}

fn locus_of(rf: &Reference, position: i64) -> i64 {
    let i = find(rf.locus_start, position);
    if i >= 0 && position < rf.locus_end[i as usize] { i } else { -1 }
}

fn add_candidate(q: &[u8], rf: &Reference, w: &Work, n: usize, bonus: i64, c: &mut Candidates) {
    let (matches, mismatches) = score_blocks(q, rf, &w.rs, &w.re, &w.dg, n);
    let span = w.dg[n - 1] + w.re[n - 1] - (w.dg[0] + w.rs[0]);
    let length_penalty = (((span.max(1) as f64).log2() * -0.25 - 0.5).ceil() as i64).min(0);
    let score = matches - mismatches + bonus + length_penalty;
    let (start, end) = (w.dg[0] + w.rs[0], w.dg[n - 1] + w.re[n - 1]);
    for j in 0..c.n {
        if c.start[j] == start && c.end[j] == end && c.blocks[j] == n {
            if score > c.score[j] {
                c.score[j] = score;
            }
            return;
        }
    }
    if c.n >= MAX_CANDIDATES {
        return;
    }
    let j = c.n;
    c.score[j] = score;
    c.matches[j] = matches;
    c.mm[j] = mismatches;
    c.start[j] = start;
    c.end[j] = end;
    c.blocks[j] = n;
    c.rs[j][..n].copy_from_slice(&w.rs[..n]);
    c.re[j][..n].copy_from_slice(&w.re[..n]);
    c.dg[j][..n].copy_from_slice(&w.dg[..n]);
    c.unique[j] = genome_unique(q, rf, &w.rs, &w.re, &w.dg, n);
    c.n += 1;
}

fn align_one(q: &[u8], rf: &Reference, w: &mut Work, c: &mut Candidates) {
    let length = q.len();
    let k = rf.k;
    w.hits_diag.clear();
    w.hits_read.clear();
    let mut p = 0usize;
    while p + k <= length {
        let mut code = 0u64;
        let mut good = true;
        for &base in &q[p..p + k] {
            if base > 3 {
                good = false;
                break;
            }
            code = (code << 2) | base as u64;
        }
        if good {
            let (first, last) = lookup(code, rf);
            if last > first && last - first <= MAX_OCCURRENCES {
                for j in first..last {
                    if w.hits_diag.len() < MAX_HITS {
                        w.hits_diag.push(rf.index_positions[j] - p as i64);
                        w.hits_read.push(p as i64);
                    }
                }
            }
        }
        if p == length - k {
            break;
        }
        p = (p + STRIDE).min(length - k);
    }
    if w.hits_diag.is_empty() {
        return;
    }
    w.hit_order.clear();
    w.hit_order.extend(0..w.hits_diag.len());
    {
        let (hd, hr) = (&w.hits_diag, &w.hits_read);
        w.hit_order.sort_by_key(|&h| hd[h] * 256 + hr[h]);
    }
    w.seg_diag.clear();
    w.seg_start.clear();
    w.seg_end.clear();
    for idx in 0..w.hit_order.len() {
        let h = w.hit_order[idx];
        let (d, r) = (w.hits_diag[h], w.hits_read[h]);
        let n = w.seg_diag.len();
        if n > 0 && w.seg_diag[n - 1] == d {
            w.seg_end[n - 1] = w.seg_end[n - 1].max(r + k as i64);
        } else if n < MAX_SEGMENTS {
            w.seg_diag.push(d);
            w.seg_start.push(r);
            w.seg_end.push(r + k as i64);
        }
    }
    let n_seg = w.seg_diag.len();
    // transcript-space candidates
    for s in 0..n_seg {
        let d = w.seg_diag[s];
        let position = d + w.seg_start[s];
        if position < rf.tx_space_start || c.n >= MAX_CANDIDATES {
            continue;
        }
        let t = find(rf.tx_ref_start, position);
        if t < 0 || position >= rf.tx_ref_start[t as usize] + rf.tx_length[t as usize] {
            continue;
        }
        let t = t as usize;
        let (seg_start, seg_end) = (w.seg_start[s], w.seg_end[s]);
        let n = transcript_candidate(q, rf, d, seg_start, seg_end, rf.tx_ref_start[t], rf.tx_ref_start[t] + rf.tx_length[t], t, w);
        if n == 0 {
            continue;
        }
        let bonus = 2 * (n as i64 - 1);
        add_candidate(q, rf, w, n, bonus, c);
    }
    // locus-space chaining
    let mut by_start: Vec<usize> = (0..n_seg).filter(|&s| w.seg_diag[s] + w.seg_start[s] < rf.tx_space_start).collect();
    let n_locus = by_start.len();
    if n_locus == 0 {
        return;
    }
    by_start.sort_by_key(|&s| w.seg_start[s]);
    let mut dp = vec![0i64; n_locus];
    let mut previous = vec![-1i64; n_locus];
    let mut loci = vec![0i64; n_locus];
    for j in 0..n_locus {
        let sj = by_start[j];
        loci[j] = locus_of(rf, w.seg_diag[sj] + w.seg_start[sj]);
        dp[j] = w.seg_end[sj] - w.seg_start[sj];
        for i in 0..j {
            let si = by_start[i];
            if loci[i] != loci[j] || loci[j] < 0 {
                continue;
            }
            let gap = w.seg_diag[sj] - w.seg_diag[si];
            if gap == 0 || gap < -MAX_INDEL || gap > MAX_INTRON || (gap > 0 && gap < MIN_INTRON && gap > MAX_INDEL) {
                continue;
            }
            if w.seg_start[si] >= w.seg_start[sj] || w.seg_end[si] > w.seg_end[sj] {
                continue;
            }
            let value = dp[i] + w.seg_end[sj] - w.seg_start[sj].max(w.seg_end[si]);
            if value > dp[j] {
                dp[j] = value;
                previous[j] = i as i64;
            }
        }
    }
    let mut used = vec![false; n_locus];
    let mut ends: Vec<usize> = (0..n_locus).collect();
    ends.sort_by_key(|&j| -dp[j]);
    let mut chain: Vec<usize> = Vec::with_capacity(n_locus);
    for &e in &ends {
        if used[e] || c.n >= MAX_CANDIDATES || dp[e] < (2 * k as i64) / 2 {
            continue;
        }
        chain.clear();
        let mut j = e as i64;
        while j >= 0 {
            used[j as usize] = true;
            chain.push(by_start[j as usize]);
            j = previous[j as usize];
        }
        chain.reverse();
        let locus = loci[e] as usize;
        let (l_start, l_end) = (rf.locus_start[locus], rf.locus_end[locus]);
        let mut n = 0usize;
        let mut penalty = 0i64;
        let first = chain[0];
        w.dg[0] = w.seg_diag[first];
        w.rs[0] = extend_left(q, rf, w.seg_diag[first], w.seg_start[first], 0i64.max(l_start - w.seg_diag[first]));
        let mut ok = true;
        for ci in 1..chain.len() {
            let (a, b) = (chain[ci - 1], chain[ci]);
            let gap = w.seg_diag[b] - w.seg_diag[a];
            let mut low = w.seg_end[a].min(w.seg_start[b]);
            let mut high = w.seg_end[a].max(w.seg_start[b]);
            low = (low - 4).max(w.seg_start[a] + 1).max(w.rs[n] + 1);
            high = (high + 4).min(w.seg_end[b] - 1);
            if gap < 0 {
                high = high.min(w.seg_end[b] + gap - 1);
            }
            if low > high {
                ok = false;
                break;
            }
            let joined = join(q, rf, w.seg_diag[a], w.seg_diag[b], low, high, l_start, l_end);
            let (split, m, annotated) = match joined {
                Some(value) if n + 1 < MAX_BLOCKS => value,
                _ => {
                    ok = false;
                    break;
                }
            };
            w.re[n] = split;
            if gap >= MIN_INTRON {
                penalty -= junction_score(m, annotated);
                n += 1;
                w.rs[n] = split;
            } else if gap > 0 {
                penalty += 2 + 2 * gap;
                n += 1;
                w.rs[n] = split;
            } else {
                penalty += 2 - 2 * gap;
                n += 1;
                w.rs[n] = split - gap;
            }
            w.dg[n] = w.seg_diag[b];
        }
        if !ok {
            continue;
        }
        let last = chain[chain.len() - 1];
        w.re[n] = extend_right(q, rf, w.seg_diag[last], w.seg_end[last], (length as i64).min(l_end - w.seg_diag[last]));
        n += 1;
        if n < MAX_BLOCKS {
            let (split, acceptor, delta) = rescue_right(q, rf, w.dg[n - 1], w.rs[n - 1], w.re[n - 1], l_end);
            if delta > 0 {
                w.re[n - 1] = split;
                w.rs[n] = split;
                w.re[n] = length as i64;
                w.dg[n] = acceptor - split;
                n += 1;
            }
        }
        if n < MAX_BLOCKS {
            let (split, donor, delta) = rescue_left(q, rf, w.dg[0], w.rs[0], w.re[0], l_start);
            if delta > 0 {
                for i in (1..=n).rev() {
                    w.rs[i] = w.rs[i - 1];
                    w.re[i] = w.re[i - 1];
                    w.dg[i] = w.dg[i - 1];
                }
                w.rs[1] = split;
                w.rs[0] = 0;
                w.re[0] = split;
                w.dg[0] = donor - split;
                n += 1;
            }
        }
        if (0..n).any(|i| w.re[i] <= w.rs[i]) {
            continue;
        }
        add_candidate(q, rf, w, n, -penalty, c);
    }
}

struct Output {
    status: *mut i8,
    orientation: *mut i8,
    scores: *mut i16,
    matched: *mut i16,
    mismatched: *mut i16,
    n_blocks: *mut i8,
    blocks: *mut i64,
}

unsafe impl Send for Output {}
unsafe impl Sync for Output {}

fn align_read(q: &[u8], reverse: &mut [u8], rf: &Reference, w: &mut Work, out: &Output, read: usize) {
    let length = q.len();
    for i in 0..length {
        let b = q[length - 1 - i];
        reverse[i] = if b <= 3 { 3 - b } else { 4 };
    }
    let mut c = Candidates::new();
    for o in 0..2i8 {
        let before = c.n;
        if o == 0 {
            align_one(q, rf, w, &mut c);
        } else {
            align_one(reverse, rf, w, &mut c);
        }
        for j in before..c.n {
            c.orient[j] = o;
        }
    }
    if c.n == 0 {
        return;
    }
    let mut best = 0usize;
    for j in 1..c.n {
        if c.score[j] > c.score[best] {
            best = j;
        }
    }
    let minimum = (0.66 * (length as f64 - 1.0)) as i64;
    // SAFETY: each read index is written by exactly one thread.
    unsafe {
        if c.score[best] < minimum || c.matches[best] < minimum || c.mm[best] > 10 || (c.mm[best] as f64) > 0.3 * ((c.matches[best] + c.mm[best]) as f64) {
            *out.status.add(read) = TOO_SHORT;
            return;
        }
        let mut multi = false;
        for j in 0..c.n {
            if j != best && c.score[j] >= c.score[best] - 1 && (c.end[j] <= c.start[best] || c.start[j] >= c.end[best] || c.orient[j] != c.orient[best]) {
                multi = true;
            }
        }
        *out.orientation.add(read) = c.orient[best];
        *out.scores.add(read) = c.score[best] as i16;
        *out.matched.add(read) = c.matches[best] as i16;
        *out.mismatched.add(read) = c.mm[best] as i16;
        *out.n_blocks.add(read) = c.blocks[best] as i8;
        for i in 0..c.blocks[best] {
            *out.blocks.add((read * MAX_BLOCKS + i) * 2) = c.dg[best][i] + c.rs[best][i];
            *out.blocks.add((read * MAX_BLOCKS + i) * 2 + 1) = c.dg[best][i] + c.re[best][i];
        }
        *out.status.add(read) = if multi { MULTI } else if !c.unique[best] { NOT_GENOME_UNIQUE } else { UNIQUE };
    }
}

unsafe fn slice<'a, T>(pointer: *const T, n: usize) -> &'a [T] {
    if n == 0 { &[] } else { std::slice::from_raw_parts(pointer, n) }
}

/// Align `n_reads` reads of `length` 2-bit codes (row-major) to the target
/// reference. Output arrays must be zero-initialised: status, orientation
/// (int8), scores, matched, mismatched (int16), n_blocks (int8) and blocks
/// (int64, n_reads x 8 x 2). Returns 0.
///
/// # Safety
/// All pointers must be valid for the given lengths.
#[no_mangle]
#[allow(clippy::too_many_arguments)]
pub unsafe extern "C" fn align_reads(
    n_reads: usize, length: usize, reads: *const u8,
    seq: *const u8, n_seq: usize,
    index_codes: *const u64, index_positions: *const i64, n_index: usize,
    buckets: *const i64, unique: *const u8, k: usize,
    locus_start: *const i64, locus_end: *const i64, n_loci: usize,
    tx_space_start: i64, tx_ref_start: *const i64, tx_length: *const i64, n_tx: usize,
    exon_offsets: *const i64, exon_tx_start: *const i64, exon_ref_start: *const i64, exon_length: *const i64, n_exons: usize,
    junctions: *const i64, n_junctions: usize, donors: *const i64, n_donors: usize, acceptors: *const i64, n_acceptors: usize,
    n_threads: usize, chunk: usize,
    status: *mut i8, orientation: *mut i8, scores: *mut i16, matched: *mut i16, mismatched: *mut i16, n_blocks: *mut i8, blocks: *mut i64,
) -> i32 {
    let rf = Reference {
        seq: slice(seq, n_seq),
        index_codes: slice(index_codes, n_index),
        index_positions: slice(index_positions, n_index),
        buckets: slice(buckets, (1usize << BUCKET_BITS) + 1),
        unique: slice(unique, n_seq),
        k,
        shift: 2 * k as u64 - BUCKET_BITS,
        locus_start: slice(locus_start, n_loci),
        locus_end: slice(locus_end, n_loci),
        tx_space_start,
        tx_ref_start: slice(tx_ref_start, n_tx),
        tx_length: slice(tx_length, n_tx),
        exon_offsets: slice(exon_offsets, n_tx + 1),
        exon_tx_start: slice(exon_tx_start, n_exons),
        exon_ref_start: slice(exon_ref_start, n_exons),
        exon_length: slice(exon_length, n_exons),
        junctions: slice(junctions, n_junctions),
        donors: slice(donors, n_donors),
        acceptors: slice(acceptors, n_acceptors),
    };
    let reads = slice(reads, n_reads * length);
    let out = Output { status, orientation, scores, matched, mismatched, n_blocks, blocks };
    let next = AtomicUsize::new(0);
    let chunk = chunk.max(1);
    std::thread::scope(|scope| {
        for _ in 0..n_threads.max(1) {
            scope.spawn(|| {
                let mut w = Work::new();
                let mut reverse = vec![0u8; length];
                loop {
                    let first = next.fetch_add(chunk, Ordering::Relaxed);
                    if first >= n_reads {
                        break;
                    }
                    for read in first..(first + chunk).min(n_reads) {
                        align_read(&reads[read * length..(read + 1) * length], &mut reverse, &rf, &mut w, &out, read);
                    }
                }
            });
        }
    });
    0
}
