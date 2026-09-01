# Analysis

What was built, why it was built that way, and the line between what was measured and what
was simulated.

## What could not be run

DeepVariant and its pangenome-aware variant need a container runtime and a reference
genome. `main.nf` is written and unexecuted. **No number in this repository compares two
real callers**, and the README says so before it says anything else.

What is measured instead is the evaluation layer — which is what the repository is actually
about. The design notes are entirely about how variant-calling comparisons go wrong; the
callers were only ever the vehicle.

## Two claims, two kinds of evidence

Keeping these apart is the point of the write-up:

1. **The stratification of GIAB chr20 is measurement.** Published truth data, published
   stratification BEDs, no modelling.
2. **The cost of the two evaluation mistakes is simulation.** Synthetic call sets generated
   from the real truth coordinates under stated error models.

A reader who conflates them would come away thinking this repository benchmarked
DeepVariant. It did not.

## Design decisions, and the reasoning

**Restriction happens before any counting.** A call outside the high-confidence regions is
not a false positive — it is unknown. `restrict` drops those calls and `excluded_count`
reports them, so no code path can treat an unevaluable call as an error.

**Per stratum first, aggregate second.** Most of a chromosome is easy unique sequence, so a
real improvement in hard regions is diluted into invisibility by an aggregate.

**First-match stratum assignment, ordered.** A site can be both a segmental duplication and
low-mappability. Counting it twice would inflate whichever stratum is listed later, so
`STRATA` defines a fixed precedence.

**The MHC stratum stays in the table at zero.** The MHC is on chr6, so chr20 has none.
Silently dropping an empty stratum would produce a table whose columns differ between
chromosomes and cannot be compared across runs.

**An interval index, added because real data demanded it.** The original lookup was a linear
scan over every confident region for every call — 86,000 calls against 10,000 regions is
nearly a billion comparisons. `IntervalIndex` merges overlapping intervals once and bisects.
The naive implementation is kept as the reference the fast one is tested against.

## What was measured

84,998 biallelic truth variants on chr20. **90.4% are in unique sequence.** High-confidence
regions cover 56,916,239 bases — 88.3% of the chromosome.

| Stratum | Variants | Bases | Share of chr20 |
| --- | ---: | ---: | ---: |
| unique | 76,826 | — | — |
| homopolymer | 3,747 | 957,549 | 1.49% |
| low mappability | 2,709 | 5,161,676 | 8.01% |
| segmental duplication | 1,716 | 2,885,987 | 4.48% |
| MHC | 0 | 0 | 0.00% |

## What was simulated, and what it shows

**Aggregate F1 hides the effect.** The pangenome-like caller beats the linear-like one by
0.0043 in aggregate and by **0.1492 in segmental duplications — 35× larger**. A reader given
only the aggregate would see a fourth-decimal difference and call the two equivalent.

**Skipping restriction penalises the better caller more.** Precision loss from forgetting to
restrict: −0.0015 for the linear-like caller, **−0.0068 for the pangenome-like one, 4.5×
worse**. The bias has a direction and it points against the caller a pangenome exists to
enable.

## What depends on the assumed error models, and what does not

The magnitudes do. The **direction** of both effects does not: it follows from the structure
of the chromosome — 90% easy sequence, 12% of it outside the confident regions — and that
structure is real, measured, and in the section above.

## What would change the conclusion

Docker and a reference genome. Running the real callers would replace the simulated
magnitudes with measured ones. It would not change the evaluation-procedure findings, which
are properties of the chromosome and of the scoring rules, not of any caller.
