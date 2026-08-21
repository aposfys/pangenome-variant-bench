"""Per-stratum precision, recall and F1 from normalised call sets.

Call sets arriving here must already be normalised by an external tool. Comparing raw VCF
positions treats representation differences -- the same indel written two ways -- as
disagreements, which invents a false positive and a false negative for a variant both
callers found.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from panbench.strata import STRATA, Call, Region, assign_stratum, restrict


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
    confident: Sequence[Region],
) -> Counts:
    """Compare one normalised call set against truth, inside confident regions only."""
    excluded = len(query) - len(restrict(query, confident))
    query_set = set(restrict(query, confident))
    truth_set = set(restrict(truth, confident))
    return Counts(
        true_positives=len(query_set & truth_set),
        false_positives=len(query_set - truth_set),
        false_negatives=len(truth_set - query_set),
        excluded=excluded,
    )


def compare_by_stratum(
    query: Sequence[Call],
    truth: Sequence[Call],
    confident: Sequence[Region],
    memberships: dict[str, Sequence[Region]],
) -> dict[str, Counts]:
    """The primary output: one confusion table per stratum, aggregate reported separately."""
    results: dict[str, Counts] = {}
    for stratum in STRATA:
        stratum_query = [c for c in query if assign_stratum(c, memberships) == stratum]
        stratum_truth = [c for c in truth if assign_stratum(c, memberships) == stratum]
        if not stratum_query and not stratum_truth:
            continue
        results[stratum] = compare(stratum_query, stratum_truth, confident)
    return results
