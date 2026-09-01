# pangenome-variant-bench
What does a variant-calling comparison hide if you evaluate it wrong?

[![CI](https://github.com/aposfys/pangenome-variant-bench/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/pangenome-variant-bench/actions/workflows/ci.yml)
[![Nextflow](https://img.shields.io/badge/nextflow-DSL2-brightgreen)](https://www.nextflow.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

A stratified evaluation harness for variant callers, built on real GIAB data and
used to price two evaluation mistakes. **Scope: this is the evaluation layer, not
a caller benchmark** — see [what is and isn't measured](#scope) below.

```
make install
panbench fetch --region chr20   # GIAB truth VCF, confident BED, stratification BEDs
panbench experiment             # the two measurements below
make test                       # 18 tests, no Nextflow and no network
```

### Two mistakes, priced

On GIAB HG002 chr20 — 84,998 biallelic truth variants, 90.4% of them in unique
sequence, high-confidence regions covering 88.3% of the chromosome:

**Aggregate F1 hides the effect.** A pangenome-like caller beats a linear-like one
by 0.0043 in aggregate and by **0.1492 in segmental duplications — 35× larger**.
An aggregate number is dominated by the 90% of the genome where the two methods
cannot differ.

**Skipping confident-region restriction has a direction.** Precision loss from
forgetting to restrict is −0.0015 for the linear-like caller and **−0.0068 for the
pangenome-like one — 4.5× worse**. The bias points against the caller a pangenome
exists to enable.

### Scope

The stratification, the truth set and the interval logic run on real GIAB data.
**The two callers are synthetic**, built from the real truth coordinates under
stated error models, so these conclusions are about the evaluation procedure and
not about DeepVariant.

`main.nf` sketches the real comparison, and it has never been run: DeepVariant is
distributed linux-64 only and `vg` has no macOS build at all, so the calling arm
needs a Linux host with a container runtime. Its process bodies are milestones,
not implementations. **No number here compares two real callers.**

### More

- [Results](results/RESULTS.md) — full tables, and what the simulation does and does not establish
- [Analysis](ANALYSIS.md) — what was done and why it was done that way
- [Design](docs/DESIGN.md) — the strata, the callers, and the traps this avoids
