# Multitarget seed-consensus pilot

This validation applies `dockconsensusR` to three fixed-seed, precomputed
score tables from the companion `dockclaimR` DUD-E pilot. For each target, the
three seeds are treated as separate ranking methods. Scores are ranked within
seed; raw score scales are not pooled.

## Fixed contract

Each target must have three seed scenarios and 24 unique ligands per scenario.
The workflow reports every pairwise Spearman rank agreement and one consensus
summary per ligand, including the fraction of seed methods in which it was
observed. A ligand absent from a seed is not assigned a poor rank.

This is a ranking-agreement exercise only. It does not use or assess reference
activity labels and cannot demonstrate binding, biological activity, efficacy,
safety, or prospective screening performance. Source data are local and not
redistributed. DUD-E target sources are documented at
<https://dude.docking.org/targets/ada>,
<https://dude.docking.org/targets/ampc>, and
<https://dude.docking.org/targets/comt>.
