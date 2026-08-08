# dockconsensusR

`dockconsensusR` summarises agreement among precomputed docking rankings while
keeping score scales separate. It calculates within-method ranks, coverage,
mean normalized ranks, Borda scores, and pairwise Spearman agreement.

It does not run docking or establish molecular binding, biological activity,
efficacy, safety, or prospective screening performance.

```r
library(dockconsensusR)

rankings <- data.frame(
  ligand_id = rep(c("A", "B", "C"), 2),
  method_id = rep(c("receptor_1", "receptor_2"), each = 3),
  score = c(-9, -8, -7, -7, -9, -8)
)

summarise_consensus(rankings)
summarise_rank_agreement(rankings)
```

## Current validation status

A limited fixed-seed multiblanco consensus pilot has been run against
precomputed rankings from the companion study. The protocol and conservative
result statement are in [VALIDATION-PILOT.md](VALIDATION-PILOT.md) and
[PILOT-RESULTS.md](PILOT-RESULTS.md). It is not a biological or prospective
screening validation.
