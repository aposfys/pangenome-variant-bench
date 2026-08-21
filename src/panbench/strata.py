"""Confident-region restriction and stratification.

The single most consequential step in a variant-calling comparison is deciding which calls
are allowed to count. A call outside the truth set's high-confidence regions is not wrong --
it is unevaluable -- and scoring it as a false positive punishes the caller that was braver
in hard sequence. That is the opposite of what a pangenome comparison is trying to measure.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

#: Strata reported separately, because aggregate F1 is dominated by easy sequence.
STRATA = ("segmental_duplication", "homopolymer", "mhc", "low_mappability", "unique")


@dataclass(frozen=True)
class Region:
    """Half-open interval, BED convention: ``start`` inclusive, ``end`` exclusive."""

    chrom: str
    start: int
    end: int

    def __post_init__(self) -> None:
        if self.end <= self.start:
            raise ValueError(f"empty or reversed region {self.chrom}:{self.start}-{self.end}")

    def contains(self, chrom: str, position: int) -> bool:
        return chrom == self.chrom and self.start <= position < self.end


@dataclass(frozen=True)
class Call:
    """One variant call, reduced to what the comparison needs."""

    chrom: str
    position: int
    ref: str
    alt: str


def in_confident_regions(call: Call, regions: Sequence[Region]) -> bool:
    """Whether a call falls inside any high-confidence region."""
    return any(region.contains(call.chrom, call.position) for region in regions)


def restrict(calls: Iterable[Call], regions: Sequence[Region]) -> list[Call]:
    """Drop calls that cannot be evaluated.

    Dropped, not counted. This function exists so that no code path can accidentally treat
    an unevaluable call as an error.
    """
    return [call for call in calls if in_confident_regions(call, regions)]


def excluded_count(calls: Sequence[Call], regions: Sequence[Region]) -> int:
    """How many calls were unevaluable. Reported, never silently discarded."""
    return len(calls) - len(restrict(calls, regions))


def assign_stratum(call: Call, memberships: dict[str, Sequence[Region]]) -> str:
    """Assign a call to the first stratum containing it, defaulting to ``unique``.

    Order matters and follows :data:`STRATA`: a site can be both a segmental duplication
    and low-mappability, and reporting it twice would inflate whichever stratum is listed
    later.
    """
    for stratum in STRATA:
        regions = memberships.get(stratum, ())
        if any(region.contains(call.chrom, call.position) for region in regions):
            return stratum
    return "unique"
