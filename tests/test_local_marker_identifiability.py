import pytest

from tealeaf.sc.local_marker_identifiability import transcript_marker_classes
from tealeaf.sc.local_read_support import event_local_read_contrast


def cassette():
    chains = dict(included=((10, 20), (30, 40), (50, 60)), excluded=((10, 20), (50, 60)), outside=((10, 20), (30, 40), (70, 80)))
    event = event_local_read_contrast('e', 'g', dict(chromosome='chr1', strand='+', transcripts=chains), ['included'], ['excluded'], (19, 51))
    return event, chains


@pytest.mark.parametrize('junction_only', [False, True])
def test_outside_event_transcript_can_produce_class_specific_markers(junction_only):
    event, chains = cassette()
    classes = transcript_marker_classes(event, chains, junction_only)
    assert classes['included'] == (True, False)
    assert classes['excluded'] == (False, True)
    assert classes['outside'] == (True, False)


def test_whole_transcript_both_marker_classes_is_not_claimed_as_one_read():
    chains = dict(included=((10, 20), (30, 40), (70, 80)), excluded=((10, 20), (50, 60), (70, 80)), both=((10, 20), (30, 40), (50, 60), (70, 80)))
    event = event_local_read_contrast('e', 'g', dict(chromosome='chr1', strand='+', transcripts=chains), ['included'], ['excluded'], (19, 71))
    assert transcript_marker_classes(event, chains)['both'] == (True, True)


def test_unrelated_exons_neither_class_and_genomic_order_is_preserved():
    event, chains = cassette()
    forward = transcript_marker_classes(event, chains)
    backward = transcript_marker_classes(event, {key: chain[::-1] for key, chain in chains.items()})
    assert backward == forward
    assert transcript_marker_classes(event, dict(neither=((90, 100),)))['neither'] == (False, False)


@pytest.mark.parametrize('chain', [(), ((20, 10),), ((0, 20), (10, 30))])
def test_invalid_transcript_annotation_is_rejected(chain):
    event, _ = cassette()
    with pytest.raises(ValueError):
        transcript_marker_classes(event, dict(invalid=chain))
