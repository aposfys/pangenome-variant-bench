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
    """Biallelic records from a VCF, with positions converted to 0-based.

    VCF ``POS`` is 1-based and BED intervals are 0-based half-open. Every containment test
    in :mod:`panbench.strata` uses the BED convention, so ``POS`` is shifted by one here,
    once, at the only place VCF coordinates enter the package.
    """
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
            Call(chrom=fields[0], position=int(fields[1]) - 1, ref=fields[3], alt=fields[4])
        )
    return calls


def read_contig_length(path: Path, region: str) -> int:
    """The length of ``region`` from the VCF ``##contig`` header.

    Read from the data rather than hard-coded, so that a run on any chromosome places its
    synthetic calls on that chromosome and reports its own confident fraction.
    """
    with path.open() as handle:
        for line in handle:
            if not line.startswith("##"):
                break
            if not line.startswith("##contig=<"):
                continue
            fields = dict(
                item.split("=", 1)
                for item in line.strip()[len("##contig=<") : -1].split(",")
                if "=" in item
            )
            if fields.get("ID") == region and "length" in fields:
                return int(fields["length"])
    raise RuntimeError(f"no ##contig length for {region} in the header of {path}")


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
        position = rng.randrange(0, region_length)
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

    Query calls outside the confident regions are scored against the truth VCF as if it
    were complete there. Those matching a truth record count as true positives and the rest
    as false positives. This function exists to quantify that mistake and is never used
    for a reported result.
    """
    query_set = set(query)
    truth_set = set(truth)
    return Counts(
        true_positives=len(query_set & truth_set),
        false_positives=len(query_set - truth_set),
        false_negatives=len(truth_set - query_set),
        excluded=0,
    )


def derive_headline(
    linear: tuple[Counts, Counts, dict[str, Counts]],
    pangenome: tuple[Counts, Counts, dict[str, Counts]],
) -> dict:
    """The differences and ratios the write-up quotes, from unrounded counts.

    Computed here rather than by subtracting the rounded values in ``arms``, which shifts
    the last digit.
    """
    lin_restricted, lin_unrestricted, lin_strata = linear
    pan_restricted, pan_unrestricted, pan_strata = pangenome

    aggregate_gain = pan_restricted.f1 - lin_restricted.f1
    segdup = "segmental_duplication"
    segdup_gain = (
        pan_strata[segdup].f1 - lin_strata[segdup].f1
        if segdup in pan_strata and segdup in lin_strata
        else None
    )
    lin_loss = lin_restricted.precision - lin_unrestricted.precision
    pan_loss = pan_restricted.precision - pan_unrestricted.precision

    def outside(restricted: Counts, unrestricted: Counts) -> dict[str, int]:
        return {
            "calls_outside_confident": restricted.excluded,
            "scored_as_true_positive_unrestricted": (
                unrestricted.true_positives - restricted.true_positives
            ),
            "scored_as_false_positive_unrestricted": (
                unrestricted.false_positives - restricted.false_positives
            ),
        }

    return {
        "aggregate_f1_gain": round(aggregate_gain, 4),
        "segdup_f1_gain": None if segdup_gain is None else round(segdup_gain, 4),
        "f1_gain_by_stratum": {
            name: round(pan_strata[name].f1 - lin_strata[name].f1, 4)
            for name in STRATA
            if name in pan_strata and name in lin_strata
        },
        "segdup_to_aggregate_ratio": (
            None
            if segdup_gain is None or aggregate_gain == 0
            else round(segdup_gain / aggregate_gain, 1)
        ),
        "precision_loss_unrestricted": {
            "linear_like": round(lin_loss, 4),
            "pangenome_like": round(pan_loss, 4),
        },
        "precision_loss_ratio": None if lin_loss == 0 else round(pan_loss / lin_loss, 2),
        "outside_confident": {
            "linear_like": outside(lin_restricted, lin_unrestricted),
            "pangenome_like": outside(pan_restricted, pan_unrestricted),
        },
    }


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

    region_length = read_contig_length(sliced / f"truth.{region}.vcf", region)

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
    scored: dict[str, tuple[Counts, Counts, dict[str, Counts]]] = {}
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
        scored[model.name] = (restricted, unrestricted, by_stratum)
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

    derived = derive_headline(scored[LINEAR_LIKE.name], scored[PANGENOME_LIKE.name])

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
            "derived": derived,
        },
    }
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "findings.json").write_text(json.dumps(findings, indent=1))
    return findings
