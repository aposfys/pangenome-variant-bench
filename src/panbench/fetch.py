"""Fetch the GIAB truth set, its high-confidence regions, and the stratification BEDs.

Everything here is real published data, pinned to a release. Nothing is generated.

The download is restricted to one chromosome on the way in rather than afterwards, because
the whole-genome truth VCF is 149 MB and the point of scoping this project to chr20 is that
it should run on a laptop.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import urllib.request
from dataclasses import dataclass
from pathlib import Path

GIAB_BASE = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38"
)
STRATIFICATION_BASE = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/"
    "genome-stratifications/v3.1/GRCh38"
)

TRUTH_VCF = "HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
CONFIDENT_BED = "HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

#: The strata. Order matters: :func:`panbench.strata.assign_stratum` takes the first match,
#: so the most specific and most interesting regions come first.
STRATIFICATIONS: tuple[tuple[str, str], ...] = (
    ("mhc", "OtherDifficult/GRCh38_MHC.bed.gz"),
    ("segdup", "SegmentalDuplications/GRCh38_segdups.bed.gz"),
    ("homopolymer", "LowComplexity/GRCh38_SimpleRepeat_homopolymer_7to11_slop5.bed.gz"),
    ("low_mappability", "mappability/GRCh38_lowmappabilityall.bed.gz"),
)

#: Pinned so a rerun is comparable. A stratification set is a database, not a constant.
PINNED = {
    "truth": "GIAB HG002 NISTv4.2.1, GRCh38",
    "stratifications": "GIAB genome-stratifications v3.1, GRCh38",
}


@dataclass
class Fetched:
    """Where everything landed, and how much of it there is."""

    region: str
    truth_vcf: Path
    confident_bed: Path
    stratification_beds: dict[str, Path]
    n_truth_variants: int
    n_confident_regions: int
    confident_bases: int


def _download(url: str, path: Path) -> Path:
    if path.exists() and path.stat().st_size > 0:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_suffix(path.suffix + ".part")
    with urllib.request.urlopen(url, timeout=600) as response, partial.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    # Renamed only once complete, so an interrupted download is never mistaken for a
    # finished one on the next run.
    partial.rename(path)
    return path


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def slice_bed(source: Path, destination: Path, region: str) -> tuple[int, int]:
    """Write only ``region`` rows from a BED. Returns (rows, total bases covered)."""
    rows = 0
    bases = 0
    with _open_text(source) as handle, destination.open("w") as out:
        for line in handle:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] != region:
                continue
            out.write(f"{fields[0]}\t{fields[1]}\t{fields[2]}\n")
            rows += 1
            bases += int(fields[2]) - int(fields[1])
    return rows, bases


def slice_vcf(source: Path, destination: Path, region: str) -> int:
    """Write only ``region`` records from a VCF, keeping the header. Returns the count."""
    kept = 0
    with _open_text(source) as handle, destination.open("w") as out:
        for line in handle:
            if line.startswith("#"):
                out.write(line)
                continue
            if line.split("\t", 1)[0] != region:
                continue
            out.write(line)
            kept += 1
    return kept


def fetch(data_dir: Path, *, region: str = "chr20") -> Fetched:
    """Download and slice everything needed for one chromosome."""
    raw = data_dir / "raw"
    sliced = data_dir / region
    sliced.mkdir(parents=True, exist_ok=True)

    print(f"truth VCF ({TRUTH_VCF}, ~149 MB on first run)...", flush=True)
    truth_source = _download(f"{GIAB_BASE}/{TRUTH_VCF}", raw / TRUTH_VCF)
    truth = sliced / f"truth.{region}.vcf"
    n_variants = (
        slice_vcf(truth_source, truth, region)
        if not truth.exists()
        else sum(1 for line in truth.open() if not line.startswith("#"))
    )
    print(f"  {n_variants:,} truth variants on {region}", flush=True)

    print("high-confidence BED...", flush=True)
    confident_source = _download(f"{GIAB_BASE}/{CONFIDENT_BED}", raw / CONFIDENT_BED)
    confident = sliced / f"confident.{region}.bed"
    n_regions, confident_bases = slice_bed(confident_source, confident, region)
    print(f"  {n_regions:,} confident regions, {confident_bases:,} bases", flush=True)

    beds: dict[str, Path] = {}
    for name, remote in STRATIFICATIONS:
        print(f"stratification: {name}...", flush=True)
        source = _download(f"{STRATIFICATION_BASE}/{remote}", raw / Path(remote).name)
        target = sliced / f"{name}.{region}.bed"
        rows, bases = slice_bed(source, target, region)
        beds[name] = target
        print(f"  {rows:,} regions, {bases:,} bases", flush=True)

    return Fetched(
        region=region,
        truth_vcf=truth,
        confident_bed=confident,
        stratification_beds=beds,
        n_truth_variants=n_variants,
        n_confident_regions=n_regions,
        confident_bases=confident_bases,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="panbench.fetch", description=__doc__)
    parser.add_argument("--sample", default="HG002", help="only HG002 is wired up")
    parser.add_argument("--region", default="chr20")
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    args = parser.parse_args(argv)

    if args.sample != "HG002":
        raise SystemExit(f"only HG002 is wired up, not {args.sample!r}")
    result = fetch(args.data_dir, region=args.region)
    print(
        f"\n{result.n_truth_variants:,} variants, {result.confident_bases:,} confident bases"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
