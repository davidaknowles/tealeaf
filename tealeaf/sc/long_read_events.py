"""Assign annotated transcripts to junctions and binary event forms by their own exon chains.

Exon chains are zero-based half-open intervals sorted by genomic start, as
returned by ``read_gtf_exons``. Introns are one-based inclusive pairs
``(first intronic base, last intronic base)``, so an intron between exons
``(s_a, e_a)`` and ``(s_b, e_b)`` is ``(e_a + 1, s_b)``. A transcript counts for
a feature only through its own exon chain, never through an event catalogue's
transcript lists.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np

from .differential import read_gtf_exons


def stable_identifier(value):
    """Drop the GENCODE version suffix."""
    return str(value).split(".", 1)[0]


def chain_introns(exons):
    """Return one-based inclusive introns of a sorted zero-based exon chain."""
    return tuple((int(exons[index][1]) + 1, int(exons[index + 1][0])) for index in range(len(exons) - 1))


def transcript_exon_chains(gtf_path, transcripts=None):
    """Return ``{stable transcript ID: (chromosome, strand, exons)}`` from a GTF.

    ``transcripts`` optionally restricts the output to these (stable or
    versioned) identifiers. Exons are merged and sorted by genomic start
    whatever their order in the GTF.
    """
    wanted = None if transcripts is None else {stable_identifier(value) for value in transcripts}
    chains = {}
    for record in read_gtf_exons(gtf_path).values():
        for transcript, exons in record["transcripts"].items():
            key = stable_identifier(transcript)
            if wanted is None or key in wanted:
                chains[key] = (record["chromosome"], record["strand"], tuple(exons))
    return chains


def chain_junction_rows(chains, rows_by_transcript, wanted):
    """Map each wanted junction ``(chromosome, first, last, strand)`` to source rows.

    ``rows_by_transcript`` maps a stable transcript ID to the row indices of
    its source features. Junctions are the introns of each transcript's sorted
    exon chain, so minus-strand transcripts are handled like plus-strand ones.
    """
    result = defaultdict(set)
    for transcript, rows in rows_by_transcript.items():
        if transcript not in chains:
            continue
        chromosome, strand, exons = chains[transcript]
        for first, last in chain_introns(exons):
            junction = (chromosome, first, last, strand)
            if junction in wanted:
                result[junction].update(int(row) for row in rows)
    return {junction: np.asarray(sorted(rows), dtype=int) for junction, rows in result.items()}


@dataclass(frozen=True)
class EventForm:
    """Exon-chain requirements of one event form.

    ``introns`` must all be introns of the chain; each interval in ``covered``
    (one-based inclusive) must lie inside a single exon; ``first_intron`` and
    ``last_intron`` must be the chain's genomically first or last intron.
    """

    introns: tuple = ()
    covered: tuple = ()
    first_intron: tuple | None = None
    last_intron: tuple | None = None

    def matches(self, exons, introns, intron_set):
        if not all(intron in intron_set for intron in self.introns):
            return False
        for first, last in self.covered:
            if not any(start < first and last <= end for start, end in exons):
                return False
        if self.first_intron is not None and (not introns or introns[0] != self.first_intron):
            return False
        if self.last_intron is not None and (not introns or introns[-1] != self.last_intron):
            return False
        return True

    def coordinates(self):
        values = [value for pair in self.introns + self.covered for value in pair]
        values += [value for pair in (self.first_intron, self.last_intron) if pair is not None for value in pair]
        return values


@dataclass(frozen=True)
class BinaryEvent:
    """Inclusion and exclusion forms of a binary splicing event on one strand."""

    chromosome: str
    strand: str
    inclusion: EventForm
    exclusion: EventForm

    @property
    def span(self):
        values = self.inclusion.coordinates() + self.exclusion.coordinates()
        return min(values), max(values)

    def classify(self, exons):
        """Return 1 for the inclusion form, 0 for the exclusion form, else None."""
        introns = chain_introns(exons)
        intron_set = set(introns)
        included = self.inclusion.matches(exons, introns, intron_set)
        excluded = self.exclusion.matches(exons, introns, intron_set)
        if included == excluded:
            return None
        return 1 if included else 0


def _intron(exon_end, next_exon_start):
    """Intron between a one-based exon end and the next one-based exon start."""
    return int(exon_end) + 1, int(next_exon_start) - 1


def _junction_token(token):
    end, start = token.split("-")
    return _intron(end, start)


def suppa2_event(event_id):
    """Parse a SUPPA2 event identifier into inclusion and exclusion forms.

    Coordinates follow SUPPA2's ``generateEvents`` identifiers (one-based exon
    ends and starts). Inclusion is the form whose junctions are listed first,
    except for SE and RI where it is the exon-inclusion form. AF and AL forms
    also require their junction to be the transcript's terminal junction on the
    alternative side.
    """
    body = str(event_id).split(";", 1)[1]
    fields = body.split(":")
    event_type, chromosome, strand, tokens = fields[0], fields[1], fields[-1], fields[2:-1]
    junctions = [_junction_token(token) for token in tokens if "-" in token]
    if event_type == "SE":
        (e1, s2), (e2, s3) = (token.split("-") for token in tokens)
        inclusion = EventForm(introns=(_intron(e1, s2), _intron(e2, s3)))
        exclusion = EventForm(introns=(_intron(e1, s3),))
    elif event_type in {"A5", "A3"}:
        inclusion, exclusion = EventForm(introns=(junctions[0],)), EventForm(introns=(junctions[1],))
    elif event_type == "MX":
        inclusion, exclusion = EventForm(introns=tuple(junctions[:2])), EventForm(introns=tuple(junctions[2:]))
    elif event_type == "RI":
        intron = junctions[0]
        inclusion, exclusion = EventForm(covered=(intron,)), EventForm(introns=(intron,))
    elif event_type in {"AF", "AL"}:
        side = "first_intron" if (event_type == "AF") == (strand == "+") else "last_intron"
        inclusion, exclusion = EventForm(**{side: junctions[0]}), EventForm(**{side: junctions[1]})
    else:
        raise ValueError(f"unsupported SUPPA2 event type {event_type}")
    return BinaryEvent(chromosome, strand, inclusion, exclusion)


def rmats_event(event_type, row):
    """Parse one rMATS ``fromGTF`` row into inclusion and exclusion forms.

    rMATS starts are zero-based and ends one-based. Inclusion is the
    exon-inclusion form for SE, the long-exon form for A5SS/A3SS and the
    retained-intron form for RI. For MXE it is the first listed (genomically
    upstream) exon on the plus strand and the second listed exon on the minus
    strand, i.e. the 5' exon in transcript order; on events shared with SUPPA2
    the sign of -IncLevelDifference agrees with SUPPA2's MX effect only under
    this orientation.
    """
    value = {key: int(row[key]) for key in row if key not in {"ID", "GeneID", "geneSymbol", "chr", "strand"}}
    chromosome, strand = str(row["chr"]), str(row["strand"])
    if event_type == "SE":
        inclusion = EventForm(introns=((value["upstreamEE"] + 1, value["exonStart_0base"]), (value["exonEnd"] + 1, value["downstreamES"])))
        exclusion = EventForm(introns=((value["upstreamEE"] + 1, value["downstreamES"]),))
    elif event_type in {"A5SS", "A3SS"}:
        if value["flankingES"] >= value["longExonEnd"]:
            inclusion = EventForm(introns=((value["longExonEnd"] + 1, value["flankingES"]),))
            exclusion = EventForm(introns=((value["shortEE"] + 1, value["flankingES"]),))
        else:
            inclusion = EventForm(introns=((value["flankingEE"] + 1, value["longExonStart_0base"]),))
            exclusion = EventForm(introns=((value["flankingEE"] + 1, value["shortES"]),))
    elif event_type == "MXE":
        inclusion = EventForm(introns=((value["upstreamEE"] + 1, value["1stExonStart_0base"]), (value["1stExonEnd"] + 1, value["downstreamES"])))
        exclusion = EventForm(introns=((value["upstreamEE"] + 1, value["2ndExonStart_0base"]), (value["2ndExonEnd"] + 1, value["downstreamES"])))
        if strand == "-":
            inclusion, exclusion = exclusion, inclusion
    elif event_type == "RI":
        intron = (value["upstreamEE"] + 1, value["downstreamES"])
        inclusion, exclusion = EventForm(covered=(intron,)), EventForm(introns=(intron,))
    else:
        raise ValueError(f"unsupported rMATS event type {event_type}")
    return BinaryEvent(chromosome, strand, inclusion, exclusion)


def event_transcript_classes(events, chains):
    """Classify every transcript chain against every event.

    ``events`` is a sequence of ``BinaryEvent`` (or None); returns one pair
    ``(inclusion transcripts, exclusion transcripts)`` per event. Transcripts
    on another chromosome or strand, not overlapping the event, or consistent
    with neither form are excluded.
    """
    index = defaultdict(list)
    for transcript, (chromosome, strand, exons) in chains.items():
        index[(chromosome, strand)].append((exons[0][0] + 1, exons[-1][1], transcript, exons))
    arrays = {}
    for key, values in index.items():
        arrays[key] = (np.asarray([value[0] for value in values]), np.asarray([value[1] for value in values]), [value[2] for value in values], [value[3] for value in values])
    output = []
    for event in events:
        included, excluded = [], []
        if event is not None and (event.chromosome, event.strand) in arrays:
            starts, ends, names, exon_sets = arrays[(event.chromosome, event.strand)]
            low, high = event.span
            for position in np.flatnonzero((starts <= high) & (ends >= low)):
                label = event.classify(exon_sets[position])
                if label == 1:
                    included.append(names[position])
                elif label == 0:
                    excluded.append(names[position])
        output.append((included, excluded))
    return output
