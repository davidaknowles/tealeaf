"""Annotation-only marker scope, distinct from RNA quantification."""

from .local_read_support import alignment_local_signature


def transcript_marker_classes(event, chains, junction_only=False):
    """Return any included/excluded marker opportunity for each transcript.

    Apply the existing marker geometry to whole annotated exon chains. This
    detects potential marker sources, not read-start probabilities, measured
    contamination, or a claim that one read spans every flagged feature.
    """
    masks = (4, 8) if junction_only else (5, 10)
    result = {}
    for transcript, chain in chains.items():
        ordered = tuple(sorted(tuple(interval) for interval in chain))
        if not ordered or any(not 0 <= first < second for first, second in ordered) or any(first[1] > second[0] for first, second in zip(ordered, ordered[1:])):
            raise ValueError('nonempty ordered nonoverlapping positive exon intervals required')
        junctions = {(first[1], second[0]) for first, second in zip(ordered, ordered[1:]) if first[1] < second[0]}
        signature = alignment_local_signature(ordered, junctions, event)
        result[transcript] = tuple(bool(signature & mask) for mask in masks)
    return result
