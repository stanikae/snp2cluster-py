# SNP2Cluster-Py v1.1.1

SNP2Cluster-Py v1.1.1 is a focused visualisation hotfix for high-cardinality minimum spanning trees.

## Fixed

- Increased the supported MST sequence-type visual vocabulary from 33 to 44 unique combinations.
- Added a fourth dash-dot border-pattern family.
- Preserved node fill colour exclusively for Core SNP cluster annotation.
- Prevented MST generation from failing when a displayed analysis group contains between 34 and 44 sequence-type categories.

## Validation

The hotfix was validated using a Ghana pneumococcal dataset in which the Abura community MST contained 36 sequence-type categories.

All 36 categories were assigned unique node shape-and-border-pattern combinations while preserving Core SNP cluster fill colours and exact sequence-type values in node tooltips.

## Current upper bound

The visual encoding supports up to 44 unique sequence-type categories using 11 node shapes combined with four border-pattern variants.

Datasets exceeding 44 categories remain subject to a clear rendering error in this patch release. A non-fatal overflow strategy is planned for a subsequent feature release.
