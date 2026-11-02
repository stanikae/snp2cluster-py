# Changelog

All notable changes to this project will be documented in this file.

## [1.1.0] - 2026-11-02
### Added
- Structured run provenance in the standalone HTML report, including input files, effective analytical parameters, configured columns, validation results, platform information, and dependency versions.
- Embedded SNP2Cluster-Py branding in standalone HTML reports.
- Sticky report navigation and clustered-isolate summary percentages.
- Separate report sections for the Core SNP heatmap, transmission scatterplots, and collapsible K-selection diagnostics.
- Descriptions for machine-readable output tables.
- High-cardinality MST encoding using deterministic node-shape and border-pattern combinations, supporting up to 33 distinct sequence-type categories without using the Core SNP cluster colour channel.
- Adaptive multi-column heatmap legends for datasets containing many sequence types and Core SNP clusters.
### Fixed
- Correctly select initial K=1 when an eligible context contains multiple isolates with one distinct SNP-distance profile.
- Allow valid Core-clustering contexts and final Core SNP clusters containing exactly two isolates.
- Continue valid K=1 groups through SNP-cutoff refinement.
- Prevent all-NA silhouette diagnostics from terminating the pipeline.
- Prevent silent reuse of identical MST symbols across different sequence types.
- Prevent high-cardinality categorical heatmap legends from overlapping.
### Changed
- Clarified the two-stage Core SNP clustering workflow in the HTML report: initial data-driven K clustering followed by SNP-cutoff refinement.
- Improved scientific figure labels, result interpretation guidance, and provenance presentation.
- Final Core SNP clusters continue to require at least two isolates.
### Validation
- Validated on the bundled example dataset and a 172-isolate real-world multi-facility dataset.
- Confirmed that five clonal ST377 isolates with zero within-context SNP distances form one valid Core SNP cluster.
- Added regression tests for K=1 clustering, two-isolate contexts, diagnostic plotting, HTML reporting, high-cardinality MST symbols, and adaptive heatmap legends.

## [1.0.0] - 2026-09-29
### Added
- First public release of SNP2Cluster-Py.
- Core SNP and SNP-Epi transmission clustering.
- R-compatible analytical behavior and K-selection diagnostics.
- Ordered context-aware Core clustering (Global, ST, MainVar, SecondaryVar).
- Facility, Area, and Community analysis modes.
- Interactive minimum spanning trees and transmission networks.
- Annotated Core SNP heatmaps and transmission scatterplots.
- Machine-readable result tables, audit records, and open scientific HTML reports.
