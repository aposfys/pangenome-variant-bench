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


class IntervalIndex:
    """Merged, sorted intervals with binary-search lookup.

    The linear scan in :func:`in_confident_regions` is fine for a test and hopeless for a
    chromosome: 86,000 calls against 10,000 confident regions is nearly a billion
    comparisons. Intervals are merged once on construction -- overlapping BED rows are
    common and a merged set is both smaller and correct to bisect -- and every lookup after
    that is O(log n).
    """

    __slots__ = ("_by_chrom",)

    def __init__(self, regions: Iterable[Region]) -> None:
        grouped: dict[str, list[tuple[int, int]]] = {}
        for region in regions:
            grouped.setdefault(region.chrom, []).append((region.start, region.end))

        self._by_chrom: dict[str, tuple[list[int], list[int]]] = {}
        for chrom, spans in grouped.items():
            spans.sort()
            merged: list[tuple[int, int]] = []
            for start, end in spans:
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            self._by_chrom[chrom] = (
                [span[0] for span in merged],
                [span[1] for span in merged],
            )

    def contains(self, chrom: str, position: int) -> bool:
        import bisect

        found = self._by_chrom.get(chrom)
        if found is None:
            return False
        starts, ends = found
        # Rightmost interval whose start is <= position; it is the only candidate.
        index = bisect.bisect_right(starts, position) - 1
        return index >= 0 and position < ends[index]

    def total_bases(self) -> int:
        return sum(
            end - start
            for starts, ends in self._by_chrom.values()
            for start, end in zip(starts, ends, strict=True)
        )

    def __len__(self) -> int:
        return sum(len(starts) for starts, _ in self._by_chrom.values())


def as_index(regions: Sequence[Region] | IntervalIndex) -> IntervalIndex:
    """Accept either a plain region list or an already-built index."""
    return regions if isinstance(regions, IntervalIndex) else IntervalIndex(regions)


def restrict_indexed(
    calls: Iterable[Call], regions: Sequence[Region] | IntervalIndex
) -> list[Call]:
    """:func:`restrict`, but O(log n) per call. Same answer, different cost."""
    index = as_index(regions)
    return [call for call in calls if index.contains(call.chrom, call.position)]


def assign_stratum_indexed(
    call: Call, memberships: dict[str, Sequence[Region] | IntervalIndex]
) -> str:
    """:func:`assign_stratum`, using indexes. Same first-match ordering."""
    for stratum in STRATA:
        regions = memberships.get(stratum)
        if regions is not None and as_index(regions).contains(call.chrom, call.position):
            return stratum
    return "unique"
