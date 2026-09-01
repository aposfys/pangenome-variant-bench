"""What can be measured on this machine, and what cannot.

**The caller comparison has not been run.** DeepVariant and its pangenome-aware variant
need a container runtime and a reference genome; the Nextflow workflow in ``main.nf`` is
written and unexecuted. No number here compares two callers.

Two things are measured instead, and they are different in kind. Keep them apart:

1. **The real stratification of GIAB chr20.** Published truth data, published stratification
   BEDs, no modelling. This is measurement.
2. **The cost of skipping confident-region restriction.** A controlled experiment with an
   explicit, stated error model applied to the real truth coordinates. This is a
   *simulation*, and its conclusion is about the evaluation procedure -- which is the thing
   this repository is actually about -- not about any caller.
"""

from __future__ import annotations

import json
import random
from dataclasses import asdict, dataclass
from pathlib import Path

from panbench.compare import Counts, compare, compare_by_stratum
from panbench.strata import STRATA, Call, IntervalIndex, Region, assign_stratum_indexed


def read_bed(path: Path) -> list[Region]:
    regions: list[Region] = []
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "track", "browser")):
            continue
        chrom, start, end = line.split("\t")[:3]
        if int(end) > int(start):
            regions.append(Region(chrom=chrom, start=int(start), end=int(end)))
    return regions


def read_vcf(path: Path) -> list[Call]:
    calls: list[Call] = []
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            continue
        # Multi-allelic records are skipped rather than split: splitting without
        # left-normalising would create positions that no normalising comparison would
        # ever agree with.
        if "," in fields[4]:
            continue
        calls.append(
            Call(chrom=fields[0], position=int(fields[1]), ref=fields[3], alt=fields[4])
        )
    return calls


@dataclass
class ErrorModel:
    """A stated, auditable error model for a synthetic caller.

    Every number here is an assumption, not a measurement. They are arguments rather than
    constants so that the simulated conclusion can be re-derived under different ones.
    """

    name: str
    #: Probability of missing a true variant, per stratum.
    miss_rate: dict[str, float]
    #: False positives generated per 100 kb, inside confident regions.
    fp_per_100kb_confident: float
    #: False positives generated per 100 kb, outside confident regions.
    fp_per_100kb_outside: float


#: Stands in for a linear-reference caller: good on unique sequence, poor in hard regions,
#: and conservative about calling where it is unsure.
LINEAR_LIKE = ErrorModel(
    name="linear_like",
    miss_rate={
        "unique": 0.01,
        "homopolymer": 0.12,
        "low_mappability": 0.25,
        "segmental_duplication": 0.40,
        "mhc": 0.30,
    },
    fp_per_100kb_confident=1.0,
    fp_per_100kb_outside=2.0,
)

#: Stands in for a pangenome-aware caller: better in hard regions, and more willing to call
#: there -- which is exactly the behaviour a naive evaluation punishes.
PANGENOME_LIKE = ErrorModel(
    name="pangenome_like",
    miss_rate={
        "unique": 0.01,
        "homopolymer": 0.08,
        "low_mappability": 0.12,
        "segmental_duplication": 0.18,
        "mhc": 0.15,
    },
    fp_per_100kb_confident=1.2,
    fp_per_100kb_outside=8.0,
)


def simulate_caller(
    truth: list[Call],
    model: ErrorModel,
    memberships: dict[str, IntervalIndex],
    confident: IntervalIndex,
    *,
    region: str,
    region_length: int,
    seed: int = 0,
) -> list[Call]:
    """Generate a synthetic call set from real truth coordinates under a stated model."""
    rng = random.Random(seed)
    calls: list[Call] = []

    for call in truth:
        stratum = assign_stratum_indexed(call, memberships)
        if rng.random() >= model.miss_rate.get(stratum, 0.01):
            calls.append(call)

    truth_positions = {call.position for call in truth}
    n_confident = int(confident.total_bases() / 100_000 * model.fp_per_100kb_confident)
    n_outside = int(
        (region_length - confident.total_bases()) / 100_000 * model.fp_per_100kb_outside
    )

    made_confident = made_outside = 0
    attempts = 0
    while (made_confident < n_confident or made_outside < n_outside) and attempts < 2_000_000:
        attempts += 1
        position = rng.randrange(1, region_length)
        if position in truth_positions:
            continue
        inside = confident.contains(region, position)
        if inside and made_confident < n_confident:
            made_confident += 1
        elif not inside and made_outside < n_outside:
            made_outside += 1
        else:
            continue
        calls.append(Call(chrom=region, position=position, ref="A", alt="G"))

    return calls


def compare_unrestricted(query: list[Call], truth: list[Call]) -> Counts:
    """Score without confident-region restriction -- the mistake, measured.

    Every query call outside the confident regions becomes a false positive here, because
    nothing removed it. This function exists to quantify the damage, and is never used for
    a reported result.
    """
    query_set = set(query)
    truth_set = set(truth)
    return Counts(
        true_positives=len(query_set & truth_set),
        false_positives=len(query_set - truth_set),
        false_negatives=len(truth_set - query_set),
        excluded=0,
    )


def run(data_dir: Path, results_dir: Path, *, region: str = "chr20", seed: int = 0) -> dict:
    """Both measurements, written to ``findings.json``."""
    sliced = data_dir / region
    truth = read_vcf(sliced / f"truth.{region}.vcf")
    confident = IntervalIndex(read_bed(sliced / f"confident.{region}.bed"))

    memberships: dict[str, IntervalIndex] = {}
    for stratum, filename in (
        ("segmental_duplication", f"segdup.{region}.bed"),
        ("homopolymer", f"homopolymer.{region}.bed"),
        ("mhc", f"mhc.{region}.bed"),
        ("low_mappability", f"low_mappability.{region}.bed"),
    ):
        path = sliced / filename
        if path.exists():
            memberships[stratum] = IntervalIndex(read_bed(path))

    # chr20 length, GRCh38.
    region_length = 64_444_167

    # ---- measurement 1: the real data -------------------------------------------------
    distribution: dict[str, int] = dict.fromkeys(STRATA, 0)
    for call in truth:
        distribution[assign_stratum_indexed(call, memberships)] += 1
    confident_truth = sum(1 for call in truth if confident.contains(call.chrom, call.position))

    print(f"{len(truth):,} truth variants on {region}")
    for stratum, count in distribution.items():
        covered = memberships[stratum].total_bases() if stratum in memberships else 0
        print(
            f"  {stratum:<24} {count:>7,} variants  "
            f"{covered:>10,} bases ({covered / region_length:.2%} of {region})"
        )
    print(f"  in confident regions:  {confident_truth:,} of {len(truth):,}")

    # ---- measurement 2: the simulation ------------------------------------------------
    arms = {}
    for model in (LINEAR_LIKE, PANGENOME_LIKE):
        calls = simulate_caller(
            truth,
            model,
            memberships,
            confident,
            region=region,
            region_length=region_length,
            seed=seed,
        )
        restricted = compare(calls, truth, confident)
        unrestricted = compare_unrestricted(calls, truth)
        by_stratum = compare_by_stratum(calls, truth, confident, memberships)
        arms[model.name] = {
            "model": asdict(model),
            "n_calls": len(calls),
            "restricted": {
                **asdict(restricted),
                "precision": round(restricted.precision, 4),
                "recall": round(restricted.recall, 4),
                "f1": round(restricted.f1, 4),
            },
            "unrestricted": {
                **asdict(unrestricted),
                "precision": round(unrestricted.precision, 4),
                "recall": round(unrestricted.recall, 4),
                "f1": round(unrestricted.f1, 4),
            },
            "by_stratum": {
                name: {
                    **asdict(counts),
                    "precision": round(counts.precision, 4),
                    "recall": round(counts.recall, 4),
                    "f1": round(counts.f1, 4),
                }
                for name, counts in by_stratum.items()
            },
        }
        print(
            f"\n{model.name}: {len(calls):,} calls\n"
            f"  restricted    P {restricted.precision:.4f}  R {restricted.recall:.4f}  "
            f"F1 {restricted.f1:.4f}  (excluded {restricted.excluded:,})\n"
            f"  unrestricted  P {unrestricted.precision:.4f}  R {unrestricted.recall:.4f}  "
            f"F1 {unrestricted.f1:.4f}"
        )

    findings = {
        "caller_comparison_not_run": (
            "DeepVariant and its pangenome-aware variant need a container runtime and a "
            "reference genome. main.nf is written and unexecuted. No number here compares "
            "two real callers."
        ),
        "region": region,
        "region_length": region_length,
        "truth": {
            "n_variants": len(truth),
            "n_in_confident": confident_truth,
            "confident_bases": confident.total_bases(),
            "confident_fraction": round(confident.total_bases() / region_length, 4),
            "by_stratum": distribution,
            "stratum_bases": {
                name: index.total_bases() for name, index in memberships.items()
            },
        },
        "simulation": {
            "disclaimer": (
                "Synthetic call sets generated from the real truth coordinates under the "
                "stated error models. The conclusion is about the evaluation procedure, "
                "not about any caller."
            ),
            "seed": seed,
            "arms": arms,
        },
    }
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "findings.json").write_text(json.dumps(findings, indent=1))
    return findings
