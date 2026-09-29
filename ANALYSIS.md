# Analysis

What was built, why it was built that way, and the line between what was measured and what
was simulated.

## What could not be run

DeepVariant and its pangenome-aware variant are distributed as Linux Docker images and
need a reference genome. There is no Docker on this machine. `main.nf` is an unrun
skeleton. **No number in this repository compares two real callers**, and the README says
so before it gives any result.

What is measured instead is the evaluation layer. The design notes are about how
variant-calling comparisons go wrong, and the callers were only ever the vehicle.

## Two claims, two kinds of evidence

1. **The stratification of GIAB chr20 is measurement.** Published truth data and published
   stratification BEDs, no modelling.
2. **The cost of the two evaluation mistakes is simulation.** Synthetic call sets generated
   from the real truth coordinates under stated error models.

A reader who conflates them would come away thinking this repository benchmarked
DeepVariant. It did not.

## Design decisions, and the reasoning

**Restriction happens before any counting.** A call outside the benchmark regions is not a
false positive. It is unknown. `restrict` drops those calls and `excluded_count` reports
them, so no code path can treat an unevaluable call as an error.

**Per stratum first, aggregate second.** Most truth variants lie outside the difficult
strata, so a real improvement confined to hard regions is diluted in an aggregate.

**First-match stratum assignment, ordered.** A site can be both a segmental duplication and
low-mappability. Counting it twice would inflate whichever stratum is listed later, so
`STRATA` defines a fixed precedence. The fallback stratum, `unique` in the code, means only
"none of the four". It still holds homopolymers longer than 11 bp and other tandem repeats.

**One coordinate convention.** VCF `POS` is 1-based and BED is 0-based half-open.
`read_vcf` converts `POS` to 0-based once, and every containment test uses the BED
convention. The first version skipped this, which put variants on interval boundaries in
the wrong stratum. Fixing it moved stratum counts by up to 106 variants.

**The MHC stratum stays in the table at zero.** The MHC is on chr6, so chr20 has none.
Silently dropping an empty stratum would produce a table whose columns differ between
chromosomes.

**An interval index, added because real data demanded it.** The original lookup was a linear
scan over every confident region for every call. 86,000 calls against 10,000 regions is
nearly a billion comparisons. `IntervalIndex` merges overlapping intervals once and
bisects. The naive implementation is kept as the reference the fast one is tested against.

## What was measured

84,998 biallelic truth variants on chr20. **90.5% lie outside the four strata.** The
benchmark regions cover 56,916,239 bases, 88.3% of the chromosome. The full table is in
[results/RESULTS.md](results/RESULTS.md).

## What was simulated, and what it shows

**Aggregate F1 dilutes the effect.** The pangenome-like caller beats the linear-like one by
0.0041 in aggregate and by **0.1492 in segmental duplications, 36 times more**. A reader
given only the aggregate would see a difference in the third decimal place.

**Skipping restriction costs the pangenome-like caller more precision.** The loss is 0.0015
for the linear-like caller and **0.0068 for the pangenome-like one, 4.4 times more**.

## What depends on the assumed error models

Both the magnitudes and the directions do. The dilution appears because the pangenome-like
model is given lower miss rates only inside the four strata. The restriction penalty falls
on the pangenome-like caller because it is given four times the out-of-region
false-positive rate (8.0 against 2.0 per 100 kb). The 602 and 150 extra false positives
under unrestricted scoring are those rates applied to the 7,527,928 bases outside the
benchmark regions. Swap the rates and the penalty moves to the other caller.

What the real data fixes is the scale on which any such difference plays out. 90.5% of
truth variants lie outside the four strata and 11.7% of chr20 lies outside the benchmark
regions.

## What would change the conclusion

Docker and a reference genome. Running the real callers would replace the simulated error
models with measured behaviour, including whether a pangenome-aware caller really calls
more often outside the benchmark regions. That would settle the direction of the
restriction effect for real callers, which this simulation cannot.

## Prior work

Both measurements are established practice, not new findings. This repository implements
the GA4GH recommendation rather than discovering it.

- Krusche et al., *Nature Biotechnology* 2019 (doi 10.1038/s41587-019-0054-x), the GA4GH
  benchmarking framework. It recommends stratifying performance by variant type and genome
  context, and reports SNV concordance between two methods of 99.7% inside high-confidence
  regions against 76.5% outside, on real callers.
- Wagner et al., *Cell Genomics* 2022 (doi 10.1016/j.xgen.2022.100128), the GIAB NISTv4.2.1
  HG002 benchmark used here.
- Dwarshuis et al., *Nature Communications* 2024 (doi 10.1038/s41467-024-53260-y), the GIAB
  genome stratifications, of which v3.1 is used here.
- Zook et al., *Nature Biotechnology* 2019 (doi 10.1038/s41587-019-0074-6), the earlier GIAB
  benchmark generation, which discusses interpreting results against a truth set that is
  neither perfect nor comprehensive.
- Prodanov et al., *Bioinformatics* 2023 (doi 10.1093/bioinformatics/btad268), variant
  calling in low-copy repeats, comparing ParascopyVC with DeepVariant, GATK and FreeBayes in
  the segmental-duplication territory treated here as the headline stratum.

What this adds is an implementation that runs in one command on a laptop, with the
stratification logic under test, and a controlled demonstration that prices each mistake
against known ground truth. That is a teaching and tooling contribution. **It is not a new
result about variant calling, and the callers here are synthetic.**
