"""The evaluation invariants. These guard the number the whole comparison reports."""

from __future__ import annotations

import pytest

from panbench.compare import compare, compare_by_stratum
from panbench.strata import Call, Region, assign_stratum, excluded_count, restrict

CONFIDENT = [Region("chr20", 1000, 2000)]


def call(position: int) -> Call:
    return Call(chrom="chr20", position=position, ref="A", alt="G")


def test_calls_outside_confident_regions_are_excluded_not_penalised() -> None:
    """The most common way variant comparisons are quietly wrong."""
    query = [call(1500), call(5000)]
    truth = [call(1500)]
    counts = compare(query, truth, CONFIDENT)
    assert counts.true_positives == 1
    assert counts.false_positives == 0
    assert counts.excluded == 1


def test_exclusions_are_counted_and_reportable() -> None:
    assert excluded_count([call(1500), call(9000)], CONFIDENT) == 1
    assert restrict([call(1500), call(9000)], CONFIDENT) == [call(1500)]


def test_a_genuine_false_positive_inside_the_region_still_counts() -> None:
    counts = compare([call(1500), call(1600)], [call(1500)], CONFIDENT)
    assert counts.false_positives == 1
    assert counts.false_negatives == 0


def test_missed_truth_variants_are_false_negatives() -> None:
    counts = compare([call(1500)], [call(1500), call(1700)], CONFIDENT)
    assert counts.false_negatives == 1
    assert counts.recall == pytest.approx(0.5)


def test_perfect_agreement_scores_one() -> None:
    counts = compare([call(1500)], [call(1500)], CONFIDENT)
    assert counts.precision == 1.0
    assert counts.recall == 1.0
    assert counts.f1 == 1.0


def test_empty_comparison_returns_zero_rather_than_dividing_by_zero() -> None:
    counts = compare([], [], CONFIDENT)
    assert counts.f1 == 0.0


def test_a_site_is_assigned_to_exactly_one_stratum() -> None:
    """Overlapping strata would double-count and inflate whichever is listed later."""
    memberships = {
        "segmental_duplication": [Region("chr20", 1400, 1600)],
        "low_mappability": [Region("chr20", 1400, 1600)],
    }
    assert assign_stratum(call(1500), memberships) == "segmental_duplication"


def test_sites_in_no_stratum_fall_through_to_unique() -> None:
    assert assign_stratum(call(1500), {}) == "unique"


def test_stratified_report_separates_hard_from_easy_sequence() -> None:
    """Aggregate F1 hides the result; the per-stratum table is the deliverable."""
    memberships = {"segmental_duplication": [Region("chr20", 1400, 1600)]}
    query = [call(1500), call(1800)]
    truth = [call(1500), call(1800)]
    by_stratum = compare_by_stratum(query, truth, CONFIDENT, memberships)
    assert set(by_stratum) == {"segmental_duplication", "unique"}
    assert by_stratum["segmental_duplication"].true_positives == 1
    assert by_stratum["unique"].true_positives == 1


def test_reversed_regions_are_rejected() -> None:
    with pytest.raises(ValueError):
        Region("chr20", 2000, 1000)


def test_interval_index_agrees_with_the_linear_scan():
    """The fast path must be a pure optimisation, not a different answer."""
    import random

    from panbench.strata import IntervalIndex, restrict, restrict_indexed

    rng = random.Random(0)
    regions = [Region("chr20", i * 100, i * 100 + 37) for i in range(500)]
    calls = [Call("chr20", rng.randrange(60000), "A", "T") for _ in range(1000)]
    assert restrict(calls, regions) == restrict_indexed(calls, regions)
    # And a prebuilt index gives the same answer as building one per call.
    assert restrict_indexed(calls, IntervalIndex(regions)) == restrict_indexed(calls, regions)


def test_interval_index_merges_overlapping_regions():
    from panbench.strata import IntervalIndex

    index = IntervalIndex([Region("c", 0, 10), Region("c", 5, 20), Region("c", 30, 40)])
    assert len(index) == 2
    assert index.total_bases() == 30
    assert index.contains("c", 15)
    assert not index.contains("c", 25)


def test_interval_index_respects_half_open_bed_convention():
    from panbench.strata import IntervalIndex

    index = IntervalIndex([Region("c", 10, 20)])
    assert index.contains("c", 10)  # start is inclusive
    assert not index.contains("c", 20)  # end is exclusive
    assert not index.contains("other", 15)


def test_a_call_outside_confident_regions_is_excluded_not_counted_wrong():
    """The central guarantee: unevaluable is not the same as wrong."""
    from panbench.compare import compare

    confident = [Region("chr20", 100, 200)]
    truth = [Call("chr20", 150, "A", "T")]
    query = [Call("chr20", 150, "A", "T"), Call("chr20", 900, "G", "C")]
    counts = compare(query, truth, confident)
    assert counts.true_positives == 1
    assert counts.false_positives == 0  # the out-of-region call is not an error
    assert counts.excluded == 1  # it is reported instead
