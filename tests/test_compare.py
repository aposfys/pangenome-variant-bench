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
