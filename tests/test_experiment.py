"""The readers, the synthetic callers and the experiment, on small hand-built files."""

from __future__ import annotations

import json

import pytest

from panbench.experiment import (
    LINEAR_LIKE,
    PANGENOME_LIKE,
    read_bed,
    read_contig_length,
    read_vcf,
    run,
    simulate_caller,
)
from panbench.strata import IntervalIndex

HEADER = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chrT,length=1000,assembly=test>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)


def vcf_line(pos: int, alt: str = "G") -> str:
    return f"chrT\t{pos}\t.\tA\t{alt}\t50\tPASS\t.\n"


def test_vcf_positions_are_tested_against_bed_intervals_in_one_convention(tmp_path):
    """BED 100-200 covers VCF POS 101 to 200, and neither 100 nor 201."""
    bed = tmp_path / "c.bed"
    bed.write_text("chrT\t100\t200\n")
    vcf = tmp_path / "t.vcf"
    vcf.write_text(HEADER + "".join(vcf_line(pos) for pos in (100, 101, 200, 201)))

    index = IntervalIndex(read_bed(bed))
    inside = {
        pos: index.contains(c.chrom, c.position)
        for pos, c in zip((100, 101, 200, 201), read_vcf(vcf), strict=True)
    }
    assert inside == {100: False, 101: True, 200: True, 201: False}


def test_read_vcf_skips_multiallelic_records(tmp_path):
    vcf = tmp_path / "t.vcf"
    vcf.write_text(HEADER + vcf_line(10) + vcf_line(20, alt="G,T"))
    calls = read_vcf(vcf)
    assert [c.position for c in calls] == [9]


def test_read_bed_skips_headers_and_empty_intervals(tmp_path):
    bed = tmp_path / "c.bed"
    bed.write_text("track name=x\n#comment\nchrT\t5\t5\nchrT\t5\t9\n")
    regions = read_bed(bed)
    assert [(r.start, r.end) for r in regions] == [(5, 9)]


def test_contig_length_comes_from_the_vcf_header(tmp_path):
    vcf = tmp_path / "t.vcf"
    vcf.write_text(HEADER)
    assert read_contig_length(vcf, "chrT") == 1000
    with pytest.raises(RuntimeError, match="chrX"):
        read_contig_length(vcf, "chrX")


def _index(*spans: tuple[int, int]) -> IntervalIndex:
    from panbench.strata import Region

    return IntervalIndex([Region("chrT", start, end) for start, end in spans])


def test_simulate_caller_is_deterministic_and_stays_on_the_contig() -> None:
    from panbench.strata import Call

    truth = [Call("chrT", pos, "A", "G") for pos in range(10, 400_000, 1_000)]
    confident = _index((0, 300_000))
    kwargs = {"region": "chrT", "region_length": 400_000, "seed": 3}

    first = simulate_caller(truth, PANGENOME_LIKE, {}, confident, **kwargs)
    second = simulate_caller(truth, PANGENOME_LIKE, {}, confident, **kwargs)
    assert first == second

    truth_set = set(truth)
    synthetic = [c for c in first if c not in truth_set]
    # int(bases / 100 kb * rate), inside and outside the confident regions.
    assert sum(confident.contains(c.chrom, c.position) for c in synthetic) == int(3 * 1.2)
    assert sum(not confident.contains(c.chrom, c.position) for c in synthetic) == int(1 * 8.0)
    assert all(0 <= c.position < 400_000 for c in first)


def test_run_writes_findings_with_derived_numbers(tmp_path):
    data = tmp_path / "data" / "chrT"
    data.mkdir(parents=True)
    positions = range(1_000, 990_000, 500)
    (data / "truth.chrT.vcf").write_text(
        HEADER.replace("length=1000", "length=1000000")
        + "".join(vcf_line(pos) for pos in positions)
    )
    (data / "confident.chrT.bed").write_text("chrT\t0\t800000\n")
    (data / "segdup.chrT.bed").write_text("chrT\t100000\t300000\n")

    findings = run(tmp_path / "data", tmp_path / "results", region="chrT")

    assert findings["region_length"] == 1_000_000
    assert findings["truth"]["confident_fraction"] == 0.8
    assert json.loads((tmp_path / "results" / "findings.json").read_text()) == findings

    derived = findings["simulation"]["derived"]
    for name, model in (("linear_like", LINEAR_LIKE), ("pangenome_like", PANGENOME_LIKE)):
        outside = derived["outside_confident"][name]
        assert outside["calls_outside_confident"] == (
            outside["scored_as_true_positive_unrestricted"]
            + outside["scored_as_false_positive_unrestricted"]
        )
        # Every synthetic call outside the confident regions becomes a false positive.
        assert outside["scored_as_false_positive_unrestricted"] == int(
            2 * model.fp_per_100kb_outside
        )
