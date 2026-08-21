# pangenome-variant-bench — where does a pangenome reference actually help?

[![CI](https://github.com/aposfys/pangenome-variant-bench/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/pangenome-variant-bench/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-DSL2-brightgreen)](https://www.nextflow.io/)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Status: skeleton.** Pipeline structure and the evaluation library are in place; no run
> has been executed.

A reproducible comparison of **linear-reference** and **pangenome-aware** variant calling on
a Genome in a Bottle sample, scoped to **chromosome 20** so the whole thing runs on a laptop
in an evening rather than needing a cluster.

Pangenome-aware DeepVariant reports up to **25.5%** fewer errors than the linear-reference
version by adding pangenome haplotypes to the pileup, and the 2026 PanVariants work proposes
a best-practice framework for pangenome-based calling. Both report aggregate improvements.

**The question this repo asks instead:** *where* does the improvement live? Aggregate F1 on a
whole chromosome hides the answer, and the answer is what determines whether it's worth
changing your pipeline.

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
- **Aggregate F1 is dominated by easy sites.** Most of a chromosome is unique, high-mappability
  sequence where both callers agree, so a real improvement in hard regions is diluted into
  invisibility. Every number here is reported per stratum first, and in aggregate second.
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
nextflow.config      profiles: docker, conda, and a laptop profile with capped resources
src/panbench/
  strata.py          region stratification and confident-region restriction
  compare.py         precision / recall / F1 per stratum, from normalised call sets
  report.py          the per-stratum table and figures
```

```bash
nextflow run . -profile docker,laptop --sample HG002 --region chr20
make test            # the evaluation library, no Nextflow required
```
