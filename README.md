# pangenome-variant-bench
Where does a pangenome reference actually help?

[![CI](https://github.com/aposfys/pangenome-variant-bench/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/pangenome-variant-bench/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-DSL2-brightgreen)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **The caller comparison has not been run.** DeepVariant and its pangenome-aware variant need a container runtime and a reference genome. `main.nf` is written and unexecuted, and **no number here compares two real callers.**

What has been run is the evaluation layer, on real GIAB data — and it quantifies both of the mistakes this repository exists to avoid.

### The real chromosome

GIAB HG002 NISTv4.2.1, GRCh38 chr20, with GIAB genome-stratifications v3.1:

| Stratum | Truth variants | Bases | Share of chr20 |
| --- | ---: | ---: | ---: |
| unique | 76,826 | — | — |
| homopolymer | 3,747 | 957,549 | 1.49% |
| low mappability | 2,709 | 5,161,676 | 8.01% |
| segmental duplication | 1,716 | 2,885,987 | 4.48% |
| MHC | 0 | 0 | 0.00% |

84,998 biallelic truth variants, **90.4% of them in unique sequence**. High-confidence regions cover 88.3% of the chromosome. MHC is empty because the MHC is on chr6 — the stratum stays in the table rather than being quietly dropped, since a benchmark that omits a column on some chromosomes and not others cannot be compared across runs.

### Both failure modes, priced

Two synthetic callers built from the real truth coordinates under stated error models. **This part is a simulation** — its conclusions are about the evaluation procedure, not about DeepVariant.

**Aggregate F1 hides the effect.** The pangenome-like caller beats the linear-like one by 0.0043 in aggregate and by **0.1492 in segmental duplications — 35× larger**:

| Stratum | linear-like | pangenome-like | Δ F1 |
| --- | ---: | ---: | ---: |
| segmental duplication | 0.7432 | 0.8924 | +0.1492 |
| low mappability | 0.8550 | 0.9298 | +0.0748 |
| homopolymer | 0.9362 | 0.9569 | +0.0207 |
| unique | 0.9914 | 0.9907 | −0.0007 |
| **aggregate** | **0.9815** | **0.9858** | **+0.0043** |

**Skipping confident-region restriction penalises the better caller more.** Precision loss from forgetting to restrict: −0.0015 for the linear-like caller, **−0.0068 for the pangenome-like one — 4.5× worse**. The bias has a direction, and it points against the caller a pangenome exists to enable.

### Running it

```
make install
make data      # GIAB truth VCF, confident BED, stratification BEDs for chr20
make test      # the evaluation library, no Nextflow required

python3 -m panbench.experiment    # the two measurements above
nextflow run . -profile docker,laptop --sample HG002 --region chr20   # needs Docker
```

`make data` downloads ~180 MB once and slices everything to one chromosome on the way in.

### One change real data forced

The confident-region lookup was a linear scan over every region for every call — 86,000 calls against 10,000 regions is nearly a billion comparisons, fine in a test and hopeless on a chromosome. `IntervalIndex` merges overlapping intervals once and bisects, and a test asserts it returns exactly what the linear version does. The naive implementation is kept as the reference the fast one is checked against.

### Layout

```
main.nf              DSL2 workflow: align -> call (both ways) -> normalise -> compare
src/panbench/
  fetch.py           GIAB truth, confident regions and stratifications, sliced per chromosome
  strata.py          regions, confident-region restriction, the interval index
  compare.py         precision / recall / F1 per stratum
  experiment.py      the real stratification, and the simulation
tests/               14 tests, no Nextflow and no network
```

### More

- [Full results, including what the simulation does and does not establish](results/RESULTS.md)
- [The strata, the callers, and the traps this pipeline avoids](docs/DESIGN.md)
