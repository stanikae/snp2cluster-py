# SNP2Cluster-Py v1.1.0

SNP2Cluster-Py v1.1.0 improves the handling, visualisation, reporting, and auditability of context-aware Core SNP and transmission-cluster analyses.

## Highlights

- Correct handling of valid K=1 Core-clustering solutions
- Support for eligible contexts and final clusters containing exactly two isolates
- Robust K-selection diagnostics for clonal and identical SNP profiles
- Structured run provenance in standalone HTML reports
- Embedded SNP2Cluster-Py branding and improved report navigation
- Reorganised scientific figures and enhanced interpretation guidance
- Unique MST shape-and-border encoding for up to 33 sequence-type categories
- Adaptive multi-column heatmap legends for high-cardinality datasets

## Core-clustering corrections

Eligible contexts containing multiple isolates with one distinct SNP-distance profile now return an initial K of 1 and proceed through SNP-cutoff refinement. This supports biologically valid clonal groups without requiring an artificial multi-cluster partition.

Contexts containing exactly two isolates are now evaluated rather than skipped. A two-isolate group is retained as a Core SNP cluster only when it satisfies the configured SNP-cutoff refinement. Final reported Core SNP clusters continue to require at least two isolates.

The K-selection visualisation now handles valid K=1 results without failing on an all-NA silhouette series. A numerical silhouette width is not assigned artificially to one-cluster solutions.

## Visualisation improvements

### Minimum spanning trees

MST node fill colour remains reserved for Core SNP cluster membership.

Sequence type is now encoded using deterministic combinations of node shape and border pattern. This prevents silent reuse of identical symbols in datasets containing many sequence types and supports up to 33 unique sequence-type categories.

### Annotated Core SNP heatmaps

Categorical legends now wrap across multiple columns and adaptive side panels. This prevents overlap when datasets contain many sequence types, Core SNP clusters, transmission clusters, or other annotation categories.

## HTML reporting improvements

The standalone HTML report now includes:

- Structured run provenance
- Input files and effective analytical parameters
- Configured columns and dependency versions
- Validation status
- Embedded SNP2Cluster-Py branding
- Sticky report navigation
- Core- and transmission-clustered isolate percentages
- Separate sections for the Core SNP heatmap, transmission scatterplots, and K-selection diagnostics
- Improved methods and interpretation guidance
- Descriptions of machine-readable output tables

## Validation

This release was validated using the bundled example data and a 172-isolate multi-facility real-world dataset.

The validation included a biologically valid clonal context containing five ST377 isolates with zero within-context Core SNP distances. The five isolates were correctly retained in one Core SNP cluster after initial K=1 selection and SNP-cutoff refinement.

Regression coverage includes:

- K=1 clustering
- Identical SNP-distance profiles
- Two-isolate contexts
- SNP-cutoff refinement
- K-selection diagnostic plotting
- HTML report generation
- Structured provenance
- High-cardinality MST symbols
- Adaptive heatmap legends
- Python 3.11 compatibility

## Interpretation note

SNP2Cluster-Py identifies genomic and epidemiologically informed clusters. Minimum-spanning-tree and transmission-network edges visualise genetic relationships between isolates. These edges must not be interpreted as confirmed direct transmission events.

## Compatibility

The release retains the existing configuration-driven Core and Transmission analysis workflows and is backward compatible with Version 1 configurations.
