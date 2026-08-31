# pangenome-variant-bench
Where does a pangenome reference actually help?

[![CI](https://github.com/aposfys/pangenome-variant-bench/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/pangenome-variant-bench/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-DSL2-brightgreen)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Status: skeleton.** Pipeline structure and the evaluation library are in place; no run has been executed.

A reproducible comparison of linear-reference and pangenome-aware variant calling on a Genome in a Bottle sample, scoped to chromosome 20 so the whole thing runs on a laptop in an evening rather than needing a cluster.

Pangenome-aware DeepVariant reports large aggregate improvements. The question here is *where* the improvement lives — aggregate F1 on a whole chromosome hides the answer, and the answer is what determines whether it is worth changing your pipeline. Every number is reported per stratum first, and in aggregate second.

### Running it
```
nextflow run . -profile docker,laptop --sample HG002 --region chr20
make test            # the evaluation library, no Nextflow required
```

### Layout
```
main.nf              DSL2 workflow: align -> call (both ways) -> normalise -> compare
nextflow.config      profiles: docker, conda, and a laptop profile with capped resources
src/panbench/
  strata.py          region stratification and confident-region restriction
  compare.py         precision / recall / F1 per stratum, from normalised call sets
```
Planned: `report.py` (the per-stratum table and figures).

### Design notes
[The strata, the callers, and the traps the pipeline is built to avoid](docs/DESIGN.md)
