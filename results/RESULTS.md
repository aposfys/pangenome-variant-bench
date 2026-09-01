# Results

Run 2026-09-01.

## What was run, and what was not

**The caller comparison has not been run.** DeepVariant and its pangenome-aware variant need
a container runtime and a reference genome; neither is available here. `main.nf` is written
and unexecuted. **No number below compares two real callers.**

Two things were measured, and they are different in kind:

1. **The real stratification of GIAB chr20** — published truth data, published
   stratification BEDs, no modelling.
2. **The cost of the two evaluation mistakes this repository is about** — a controlled
   simulation with stated error models, applied to the real truth coordinates. Its
   conclusions are about the *evaluation procedure*, which is what the repo is actually
   for.

## 1. The real data

GIAB HG002 NISTv4.2.1 on GRCh38, chr20; stratifications from GIAB genome-stratifications
v3.1.

| Stratum | Truth variants | Bases | Share of chr20 |
| --- | ---: | ---: | ---: |
| unique | 76,826 | — | — |
| homopolymer | 3,747 | 957,549 | 1.49% |
| low mappability | 2,709 | 5,161,676 | 8.01% |
| segmental duplication | 1,716 | 2,885,987 | 4.48% |
| MHC | 0 | 0 | 0.00% |

84,998 biallelic truth variants; 81,850 of them (96.3%) fall inside the high-confidence
regions, which cover 56,916,239 bases — **88.3% of the chromosome**.

**MHC is empty because the MHC is on chr6.** The stratum stays in the table rather than
being quietly dropped: a benchmark that silently omits a stratum on some chromosomes and not
others is one whose columns cannot be compared across runs.

**90.4% of chr20's truth variants are in unique sequence.** That single number is why
aggregate F1 is the wrong headline, and the simulation below prices it.

## 2. The simulation

Two synthetic callers, generated from the real truth coordinates under explicitly stated
error models (`LINEAR_LIKE` and `PANGENOME_LIKE` in `experiment.py`). The pangenome-like
model misses less in hard regions and is more willing to call outside the confident
regions — the behaviour a naive evaluation punishes.

### Aggregate F1 hides the effect almost entirely

| | Aggregate F1 | Segmental-duplication F1 |
| --- | ---: | ---: |
| linear-like | 0.9815 | 0.7432 |
| pangenome-like | 0.9858 | 0.8924 |
| **Difference** | **+0.0043** | **+0.1492** |

**The improvement is 35 times larger in segmental duplications than in the aggregate.** A
reader given only the aggregate would see a difference in the fourth decimal place and
conclude the two callers are equivalent. Per stratum:

| Stratum | linear-like F1 | pangenome-like F1 | Δ |
| --- | ---: | ---: | ---: |
| segmental duplication | 0.7432 | 0.8924 | +0.1492 |
| low mappability | 0.8550 | 0.9298 | +0.0748 |
| homopolymer | 0.9362 | 0.9569 | +0.0207 |
| unique | 0.9914 | 0.9907 | −0.0007 |

The gain is concentrated exactly where a pangenome is supposed to help, and is absent — very
slightly negative — in the 90% of the chromosome that is easy.

### Skipping confident-region restriction penalises the better caller more

| | Precision, restricted | Precision, unrestricted | Loss |
| --- | ---: | ---: | ---: |
| linear-like | 0.9929 | 0.9914 | −0.0015 |
| pangenome-like | 0.9916 | 0.9848 | **−0.0068** |

**The caller more willing to call in hard regions loses 4.5× more precision** from an
evaluation that forgets to restrict. That is not noise; it is a systematic bias with a
direction, and the direction is against the caller a pangenome exists to enable. It is the
most common way a variant-calling comparison is quietly wrong.

3,657 of the pangenome-like caller's calls fall outside the confident regions. Restricted
scoring drops them as unevaluable and reports the count; unrestricted scoring counts every
one as a false positive, having no basis for that.

## What this does and does not establish

It establishes that the **evaluation procedure** in this repository behaves as designed: the
per-stratum reporting recovers an effect the aggregate hides, and confident-region
restriction removes a bias that would otherwise favour the more conservative caller.

It establishes **nothing about DeepVariant**, pangenome-aware or otherwise. The magnitudes
above follow from the assumed error models. What does not depend on those assumptions is the
*direction* of both effects, which follows from the structure of the chromosome — 90% easy
sequence, and 12% of it outside the confident regions — and that structure is real,
measured, and in section 1.
