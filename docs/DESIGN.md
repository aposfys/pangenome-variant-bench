# pangenome-variant-bench design notes

Pangenome-aware DeepVariant reports up to **25.5%** fewer errors than the linear-reference
version by adding pangenome haplotypes to the pileup (Asri et al. 2025, bioRxiv,
doi 10.1101/2025.06.05.657102). PanVariants proposes a best-practice pipeline and framework
for pangenome-based calling (Yi et al. 2026, bioRxiv, doi 10.64898/2026.04.22.720142). Both
report aggregate improvements.

The question this repo was started to ask is *where* such an improvement lives. The caller
comparison that would answer it is planned and has not been run. What exists today is the
evaluation layer in `src/panbench`, exercised on real GIAB truth data with simulated
callers (see [results/RESULTS.md](../results/RESULTS.md)).

| | Built | Planned |
| --- | --- | --- |
| **Sample** | GIAB HG002 v4.2.1, chr20 | same |
| **Callers** | two simulated callers from stated error models | DeepVariant (linear) vs pangenome-aware DeepVariant |
| **Truth** | GIAB benchmark VCF and benchmark-region BED | same |
| **Comparison** | exact match on (chrom, pos, ref, alt), biallelic records only, stratified | `hap.py` or `vcfeval` after normalisation, stratified |
| **Strata** | segmental duplications, 7 to 11 bp homopolymers, MHC, low mappability, none of these | same |
| **Workflow** | Python package and CLI | Nextflow DSL2 with containers. `main.nf` is a skeleton that runs no process yet |

## Traps the evaluation is built to avoid

- **A call outside the benchmark regions is not a false positive. It is unknown.** Scoring
  it against a truth set that does not cover it distorts precision, and the distortion
  lands on whichever caller calls more often outside those regions. Restriction happens
  before any counting, and the test suite asserts it. Built.
- **Aggregate F1 is dominated by easy sites.** Most truth variants lie outside the
  difficult strata, so an improvement confined to hard regions is diluted in the aggregate.
  Every number is reported per stratum first. Built.
- **VCF and BED use different coordinates.** VCF `POS` is 1-based and BED is 0-based
  half-open. `read_vcf` converts once, and a test pins the boundary cases. Built.
- **Representation differences are not disagreements.** The same indel can be written
  several ways in VCF, and naive comparison counts a false positive for one caller and a
  false negative for the other for a variant both found. The planned pipeline normalises
  with `bcftools norm` and compares with `hap.py` or `vcfeval`. Planned. The current Python
  comparison matches raw records, which is exact for the simulated callers because they
  copy the truth records, and would not be enough for real ones.
- **The pangenome and the truth set may share samples.** If HG002 contributed haplotypes to
  the pangenome graph, evaluating on HG002 measures memorisation. The graph must exclude the
  evaluation sample, and the pipeline should record which graph build was used. Planned.

## Layout

```
main.nf              DSL2 skeleton for align -> call (both ways) -> normalise -> compare.
                     Milestone stubs only. The workflow block invokes no process yet.
src/panbench/
  fetch.py           GIAB truth, benchmark regions and stratifications, sliced per chromosome
  strata.py          regions, benchmark-region restriction, the interval index
  compare.py         precision, recall and F1 per stratum
  experiment.py      the real stratification, and the simulation
tests/               24 tests, no Nextflow and no network
```

The benchmark-region lookup was originally a linear scan over every region for every call.
86,000 calls against 10,000 regions is nearly a billion comparisons, fine in a test and
hopeless on a chromosome. `IntervalIndex` merges overlapping intervals once and bisects,
and a test asserts it returns exactly what the linear version does.
