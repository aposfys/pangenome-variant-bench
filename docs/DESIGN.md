# pangenome-variant-bench — design notes

Pangenome-aware DeepVariant reports up to **25.5%** fewer errors than the linear-reference
version by adding pangenome haplotypes to the pileup, and the 2026 PanVariants work proposes
a best-practice framework for pangenome-based calling. Both report aggregate improvements.

**The question this repo asks instead:** *where* does the improvement live?

| | |
| --- | --- |
| **Sample** | GIAB HG002, chr20 |
| **Callers** | DeepVariant (linear) vs pangenome-aware DeepVariant |
| **Truth** | GIAB benchmark VCF + high-confidence BED |
| **Comparison** | `hap.py` / `vcfeval`, stratified |
| **Strata** | segmental duplications · homopolymers · MHC · low-mappability · everything else |
| **Workflow** | Nextflow DSL2, containerised, resumable |

## Traps this pipeline is built to avoid

- **A call outside the high-confidence regions is not a false positive — it is unknown.**
  Scoring it as an error is the most common way variant-calling comparisons are quietly
  wrong, and it systematically penalises whichever caller is more willing to call in hard
  regions, which is exactly the caller a pangenome is supposed to enable. Confident-region
  restriction happens before any counting, and the test suite asserts it.
- **Aggregate F1 is dominated by easy sites.** Most of a chromosome is unique,
  high-mappability sequence where both callers agree, so a real improvement in hard regions
  is diluted into invisibility. Every number here is reported per stratum first, and in
  aggregate second.
- **Representation differences are not disagreements.** The same indel can be written
  several ways in VCF; naive comparison counts one caller's FP and the other's FN for a
  variant both found. Comparison is done with a normalising tool, never by comparing
  positions.
- **The pangenome and the truth set may share samples.** If HG002 contributed haplotypes to
  the pangenome graph, evaluating on HG002 measures memorisation. The graph used must
  exclude the evaluation sample, and the pipeline records which graph build was used.

## Layout

```
main.nf              DSL2 workflow: align -> call (both ways) -> normalise -> compare
                     (milestone stubs; the calling arm needs a Linux host)
src/panbench/
  fetch.py           GIAB truth, confident regions and stratifications, sliced per chromosome
  strata.py          regions, confident-region restriction, the interval index
  compare.py         precision / recall / F1 per stratum
  experiment.py      the real stratification, and the simulation
tests/               18 tests, no Nextflow and no network
```

The confident-region lookup was originally a linear scan over every region for
every call — 86,000 calls against 10,000 regions is nearly a billion comparisons,
fine in a test and hopeless on a chromosome. `IntervalIndex` merges overlapping
intervals once and bisects, and a test asserts it returns exactly what the linear
version does. The naive implementation is kept as the reference the fast one is
checked against.
