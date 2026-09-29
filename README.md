# pangenome-variant-bench

A stratified evaluation harness for variant callers on GIAB HG002 chr20, used with two
simulated callers to show what aggregate F1 and unrestricted scoring hide.

[![CI](https://github.com/aposfys/pangenome-variant-bench/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/pangenome-variant-bench/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

```
make install
panbench fetch --region chr20   # GIAB truth VCF, benchmark BED, stratifications (about 191 MB)
panbench experiment             # writes results/findings.json
make test                       # 24 tests, no Nextflow and no network
```

## Scope

**The linear versus pangenome-aware DeepVariant comparison has not been run, and both
callers here are synthetic.** The truth set, the stratification and the interval logic run
on real GIAB data. The two call sets are generated from the real truth coordinates under
error models I chose, so every caller-level number below follows from those models and says
nothing about DeepVariant. `main.nf` is an unrun skeleton. DeepVariant and pangenome-aware
DeepVariant are distributed as Linux Docker images, and this machine has no Docker.

## Results

GIAB HG002 v4.2.1 chr20 has 84,998 biallelic truth variants. 90.5% of them lie outside the
four strata used here (segmental duplications, low mappability, 7 to 11 bp homopolymers,
MHC), and the benchmark regions cover 88.3% of the chromosome.

- **Aggregate F1 dilutes a gain confined to hard regions.** The simulated pangenome-like
  caller was given lower miss rates only inside the four strata. It beats the linear-like
  one by 0.0041 F1 in aggregate and by 0.1492 in segmental duplications, 36 times more.
- **Unrestricted scoring turns out-of-region calls into precision loss.** Scoring without
  the benchmark regions costs the linear-like caller 0.0015 precision and the pangenome-like
  one 0.0068, 4.4 times more. The direction comes from the error models, which give the
  pangenome-like caller four times the false-positive rate outside the benchmark regions.

Both points are established practice (GA4GH, Krusche et al. 2019). This repo implements them
and prices them on real truth coordinates. It is not a new result about variant calling.

## More

- [Results](results/RESULTS.md) with the full tables and what the simulation does and does not show
- [Analysis](ANALYSIS.md) with the design decisions and prior work
- [Design](docs/DESIGN.md) with the planned caller comparison
