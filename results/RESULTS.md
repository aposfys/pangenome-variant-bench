# Results

First run 2026-09-01. Rerun 2026-09-29 after fixing a VCF/BED coordinate off-by-one, which
moved some stratum counts by up to 106 variants. Every number below is in `findings.json`,
and the differences and ratios are in its `simulation.derived` block, computed from
unrounded counts.

## What was run, and what was not

**The caller comparison has not been run.** DeepVariant and its pangenome-aware variant are
distributed as Linux Docker images and need a reference genome, and there is no Docker
here. `main.nf` is an unrun skeleton. **No number below compares two real callers.**

Two things were measured, and they are different in kind.

1. **The real stratification of GIAB chr20.** Published truth data and published
   stratification BEDs, no modelling.
2. **The cost of two evaluation mistakes on simulated callers.** Synthetic call sets built
   from the real truth coordinates under stated error models. The magnitudes, and the
   direction of the second effect, follow from those models.

## 1. The real data

GIAB HG002 NISTv4.2.1 on GRCh38 (Wagner et al. 2022), chr20. Stratifications from GIAB
genome-stratifications v3.1 (Dwarshuis et al. 2024).

| Stratum | Truth variants | Bases | Share of chr20 |
| --- | ---: | ---: | ---: |
| none of the four (`unique` in the code) | 76,927 | n/a | n/a |
| homopolymer, 7 to 11 bp | 3,752 | 957,549 | 1.49% |
| low mappability | 2,603 | 5,161,676 | 8.01% |
| segmental duplication | 1,716 | 2,885,987 | 4.48% |
| MHC | 0 | 0 | 0.00% |

84,998 biallelic truth variants. 81,851 of them (96.3%) fall inside the benchmark regions,
which cover 56,916,239 bases, **88.3% of the chromosome**. The other 11.7% (7,527,928 bases)
is outside them.

A site in more than one stratum is assigned to the first in the order segmental
duplication, homopolymer, MHC, low mappability. The first row is whatever is left. It is
not "unique sequence" in a strict sense, because it still holds homopolymers longer than
11 bp, other tandem repeats and any difficult region the four BEDs do not cover.

**MHC is empty because the MHC is on chr6.** The stratum stays in the table so that its
columns are the same on every chromosome.

**90.5% of chr20's truth variants lie outside the four strata.** That is why an aggregate
number can hide a difference confined to them.

## 2. The simulation

Two synthetic callers, generated from the real truth coordinates under the error models
`LINEAR_LIKE` and `PANGENOME_LIKE` in `experiment.py`. Both miss 1% of variants outside the
four strata. Inside them the pangenome-like model misses fewer (18% against 40% in
segmental duplications). Outside the benchmark regions it makes 8.0 false calls per 100 kb
against 2.0 for the linear-like model.

### Aggregate F1 dilutes the effect

| | Aggregate F1 | Segmental-duplication F1 |
| --- | ---: | ---: |
| linear-like | 0.9817 | 0.7432 |
| pangenome-like | 0.9858 | 0.8924 |
| **Difference** | **+0.0041** | **+0.1492** |

**The improvement is 36 times larger in segmental duplications than in the aggregate.** A
reader given only the aggregate would see a difference in the third decimal place. The
per-stratum table follows.

| Stratum | linear-like F1 | pangenome-like F1 | Difference |
| --- | ---: | ---: | ---: |
| segmental duplication | 0.7432 | 0.8924 | +0.1492 |
| low mappability | 0.8571 | 0.9309 | +0.0738 |
| homopolymer, 7 to 11 bp | 0.9365 | 0.9565 | +0.0200 |
| none of the four | 0.9914 | 0.9907 | -0.0007 |

The gain sits in the hard strata because the error models put it there. The small loss in
the last row comes from the pangenome-like model's higher false-positive rate inside the
benchmark regions (1.2 against 1.0 per 100 kb).

### Skipping benchmark-region restriction

| | Precision, restricted | Precision, unrestricted | Loss |
| --- | ---: | ---: | ---: |
| linear-like | 0.9929 | 0.9914 | -0.0015 |
| pangenome-like | 0.9916 | 0.9848 | **-0.0068** |

**The pangenome-like caller loses 4.4 times more precision** when scoring skips the
restriction.

3,660 of the pangenome-like caller's calls fall outside the benchmark regions. Restricted
scoring drops them as unevaluable and reports the count. Unrestricted scoring counts 3,058
of them as true positives, because they match truth records that lie outside the benchmark
regions, and the other 602 as false positives (682 rise to 1,284). For the linear-like
caller the split is 3,139 calls, 2,989 true positives and 150 false positives.

The 602 and 150 are the simulated out-of-region false calls, 8.0 and 2.0 per 100 kb over
the 7,527,928 bases outside the benchmark regions. The 4.4-fold gap is therefore close to
the 4-fold gap in those assumed rates. Swap the rates and the penalty falls on the
linear-like caller instead.

## What this does and does not establish

It shows that the evaluation code behaves as designed. Per-stratum reporting recovers a
difference the aggregate hides, and restriction to the benchmark regions keeps
out-of-region calls from being scored against a truth set that does not cover them.

It establishes **nothing about DeepVariant**, pangenome-aware or otherwise, and the
direction of both effects is an input. What the real data does fix is the scale. 90.5% of
truth variants lie outside the four strata, so a gain confined to them reaches the
aggregate scaled down roughly in proportion to their 9.5% share. 11.7% of chr20 lies outside the benchmark regions, so any
difference between two callers' out-of-region call rates reaches unrestricted precision.
Whether a real pangenome-aware caller calls more often there is not measured here.
