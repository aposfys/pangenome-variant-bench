"""Per-stratum precision, recall and F1 from normalised call sets.

Call sets arriving here must already be normalised by an external tool. Comparing raw VCF
positions treats representation differences -- the same indel written two ways -- as
disagreements, which invents a false positive and a false negative for a variant both
callers found.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from panbench.strata import (
    STRATA,
    Call,
    IntervalIndex,
    Region,
    as_index,
    assign_stratum_indexed,
    restrict_indexed,
)


@dataclass(frozen=True)
class Counts:
    """Confusion counts for one stratum."""

    true_positives: int
    false_positives: int
    false_negatives: int
    excluded: int

    @property
    def precision(self) -> float:
        called = self.true_positives + self.false_positives
        return self.true_positives / called if called else 0.0

    @property
    def recall(self) -> float:
        truth = self.true_positives + self.false_negatives
        return self.true_positives / truth if truth else 0.0

    @property
    def f1(self) -> float:
        denominator = self.precision + self.recall
        if denominator == 0.0:
            return 0.0
        return 2 * self.precision * self.recall / denominator


def compare(
    query: Sequence[Call],
    truth: Sequence[Call],
    confident: Sequence[Region] | IntervalIndex,
) -> Counts:
    """Compare one normalised call set against truth, inside confident regions only.

    Restriction happens before any counting, which is the whole point: a call outside the
    confident regions is not a false positive, it is unknown, and scoring it as an error
    systematically penalises whichever caller is more willing to call in hard regions.
    """
    index = as_index(confident)
    kept_query = restrict_indexed(query, index)
    excluded = len(query) - len(kept_query)
    query_set = set(kept_query)
    truth_set = set(restrict_indexed(truth, index))
    return Counts(
        true_positives=len(query_set & truth_set),
        false_positives=len(query_set - truth_set),
        false_negatives=len(truth_set - query_set),
        excluded=excluded,
    )


def compare_by_stratum(
    query: Sequence[Call],
    truth: Sequence[Call],
    confident: Sequence[Region] | IntervalIndex,
    memberships: dict[str, Sequence[Region] | IntervalIndex],
) -> dict[str, Counts]:
    """The primary output: one confusion table per stratum, aggregate reported separately.

    Aggregate F1 is dominated by the easy majority of a chromosome, so a real improvement
    confined to hard regions is diluted into invisibility. Per stratum first, always.
    """
    index = as_index(confident)
    indexed = {name: as_index(regions) for name, regions in memberships.items()}

    # Assign once per call rather than once per call per stratum.
    query_by_stratum: dict[str, list[Call]] = {}
    truth_by_stratum: dict[str, list[Call]] = {}
    for call in query:
        query_by_stratum.setdefault(assign_stratum_indexed(call, indexed), []).append(call)
    for call in truth:
        truth_by_stratum.setdefault(assign_stratum_indexed(call, indexed), []).append(call)

    results: dict[str, Counts] = {}
    for stratum in STRATA:
        stratum_query = query_by_stratum.get(stratum, [])
        stratum_truth = truth_by_stratum.get(stratum, [])
        if not stratum_query and not stratum_truth:
            continue
        results[stratum] = compare(stratum_query, stratum_truth, index)
    return results
