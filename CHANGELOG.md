# Changelog

All notable changes to this project will be documented in this file.

## [1.1.2] - 2026-11-16
### Fixed
- Made the configured SNP cutoff authoritative for Core SNP-cluster membership by constructing threshold-graph connected components within each analytical context.
- Prevented exploratory K-means partitions and silhouette thresholds from splitting, excluding, merging, or otherwise determining final Core SNP clusters.
- Fixed valid Core clusters disappearing under `core_cluster_context: [ST]` and `[MainVar, ST]`.
- Restored consistent Hospital C cluster membership between single-facility, multi-facility, and Global analyses.
- Correctly retained all six Hospital C ST307 isolates in one Core SNP cluster at the 11-SNP cutoff.
- Correctly retained BG315, BG259, BG258, and BG254 in one four-isolate transmission cluster at the 14-day cutoff while keeping BG32 and BG317 transmission-unclustered.
- Replaced pandas date inference with explicit configured parsing for `ymd`, `ydm`, `dmy`, and `mdy`.
- Added support for slash, dash, and dot date separators and for two- and four-digit years.
- Routed transmission-network collection dates through the shared configured date parser, eliminating date-inference warnings.

### Changed
- K-means assignments and silhouette widths are now exploratory diagnostics only and no longer control authoritative Core SNP-cluster membership.
- Core SNP clusters are now defined as connected components of a graph in which isolate pairs are linked when their SNP distance is at or below the configured cutoff.
- Pairwise SNP chains are supported through single-linkage connectivity, including clusters whose indirectly connected members may exceed the cutoff.

### Validation
- Validated Hospital C under `[ST]`, `[MainVar, ST]`, and `[Global]` Core contexts.
- Validated Hospital C results in both single-facility and 172-isolate multi-facility analyses.
- Confirmed matching Hospital C transmission membership for ST152, ST17, ST307, and ST45 between single- and multi-facility inputs.
- Added regression coverage for threshold-connected components, pairwise SNP chains, disconnected components, singleton exclusion, input-order invariance, the ST307 temporal chain, explicit date parsing, and network-date parsing.
- Final automated test suite: 60 passing tests with no date-inference warnings.

## [1.1.1] - 2026-11-03
### Fixed
- Extended high-cardinality minimum-spanning-tree sequence-type encoding from 33 to 44 unique visual symbols.
- Added a fourth dash-dot border-pattern family while preserving node fill colour for Core SNP cluster annotation.
- Prevented valid MST generation from failing for datasets containing between 34 and 44 sequence-type categories in one displayed analysis group.
- Added regression coverage for 36-category MSTs and the new 44-category upper bound.
### Validation
- Validated using a Ghana pneumococcal dataset containing 36 sequence-type categories in the Abura community MST.
- Confirmed that all 36 sequence types receive unique shape-and-border combinations.

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
