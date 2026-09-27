from __future__ import annotations

from importlib import metadata
from pathlib import Path
import pandas as pd

from .audit import create_run_dirs, write_manifest, write_validation_log
from .clustering import choose_optimal_k, assign_kmeans_clusters, generate_snp_chains
from .config import SNP2ClusterConfig
from .io import read_table, read_snp_matrix, write_csv
from .mst import build_minimum_spanning_tree, mst_edges_table
from .report import write_html_report
from .transmission import build_transmission_clusters
from .validation import validate_all, ValidationError
from sklearn.preprocessing import StandardScaler
# from .r_compat import get_core_snp_clusters_r_compatible
from .r_compat import determine_optimal_k_r_style, matrix_to_mdf, get_core_snp_clusters_r_compatible
from .snp_epi import (
    run_core_snp_cluster_analysis_r_compatible,
    calculate_snp_epi_clusters_r_compatible,
)
# from .snp_epi import calculate_snp_epi_clusters_r_compatible
from .visualization.diagnostics import (
    plot_k_selection,
)
from .visualization.scatterplots import (
    generate_transmission_scatter_from_config,
)
from .visualization.heatmaps import (
    generate_annotated_heatmap_from_config,
)
from .visualization.mst import (
    generate_core_mst_from_config,
)
from .visualization.networks import (
    generate_transmission_networks_from_config,
)


# ---------------------------------------------------------
# Core context hierarchy builder
# ---------------------------------------------------------

def build_core_context_column(
    metadata: pd.DataFrame,
    context_config: list[str],
    main_var: str,
    var_01: str,
    sequence_type_col: str,
) -> pd.Series:
    """
    Build hierarchical core clustering context.

    The order supplied in context_config
    determines the hierarchy.
    """

    if (
        len(context_config) == 1
        and context_config[0] == "Global"
    ):

        return pd.Series(
            ["ALL_SAMPLES"] * len(metadata),
            index=metadata.index,
        )

    context_parts = []

    for level in context_config:

        if level == "MainVar":

            context_parts.append(
                metadata[main_var].astype(str)
            )

        elif level == "SecondaryVar":

            context_parts.append(
                metadata[var_01].astype(str)
            )

        elif level == "ST":

            context_parts.append(
                metadata[
                    sequence_type_col
                ].astype(str)
            )

        else:

            raise ValueError(
                f"Unsupported core context level: {level}"
            )

    return (
        pd.concat(
            context_parts,
            axis=1,
        )
        .astype(str)
        .agg(
            "__".join,
            axis=1,
        )
    )


def run_pipeline(config: SNP2ClusterConfig) -> dict[str, Path]:
    validation_messages = validate_all(config, strict_files=True)
    dirs = create_run_dirs(config)
    write_validation_log(dirs, validation_messages)
    write_manifest(config, dirs, validation_messages)

    if validation_messages:
        raise ValidationError("Input validation failed. See audit/validation_log.txt")

    # ---------------------------------------------------------
    # Read configuration
    # ---------------------------------------------------------

    analysis = config.analysis
    engine = analysis.get("engine", "native").lower()

    columns = config.columns
    variables = config.variables

    cluster_type = str(
        analysis.get(
            "cluster_type",
            "Core",
        )
    ).strip()

    random_seed = int(
        analysis.get(
            "random_seed",
            42,
        )
    )

    snp_cutoff = int(
        analysis.get(
            "snp_cutoff",
            20,
        )
    )

    max_k = analysis.get("max_k")
    max_k = int(max_k) if max_k is not None else None

    # ---------------------------------------------------------
    # Variables
    # ---------------------------------------------------------

    main_var = variables["main_var"]

    var_01 = variables.get("var_01")

    if cluster_type == "Transmission":

        try:
            var_02 = variables["var_02"]

        except KeyError:

            raise KeyError(
                "The configuration 'var_02' "
                "(Collection Date) is mandatory "
                "when cluster_type is 'Transmission'."
            )

    elif cluster_type == "Core":

        var_02 = variables.get("var_02")

    else:

        raise ValueError(
            f"Unknown cluster_type: "
            f"'{cluster_type}'. "
            f"Expected 'Core' or 'Transmission'."
        )

    # ---------------------------------------------------------
    # Inputs
    # ---------------------------------------------------------

    snp = read_snp_matrix(
        config.inputs["snp_matrix"]
    )

    # metadata = None

    # if cluster_type == "Transmission":

    #     metadata = read_table(
    #         config.inputs["metadata"]
    #     )


    metadata = None

    metadata_path = config.inputs.get("metadata")

    if metadata_path:
        metadata = read_table(metadata_path)

    facility_col = main_var

    facility_col = main_var
    date_col = var_02

    all_k_diagnostics = []
    all_k_assignments = []
    all_snp_clusters = []
    all_transmission_clusters = []



    mlst = None

    mlst_path = config.inputs.get("mlst_profile")

    if mlst_path:

        mlst = read_table(
            mlst_path
        )[
            [
                columns["st_sample_id"],
                columns["sequence_type"],
            ]
        ].copy()

        mlst[
            columns["st_sample_id"]
        ] = mlst[
            columns["st_sample_id"]
        ].astype(str)


    # ---------------------------------------------------------
    # Build canonical metadata object
    # ---------------------------------------------------------

    if metadata is not None:

        metadata = metadata.copy()

        metadata[
            columns["sample_id"]
        ] = (
            metadata[
                columns["sample_id"]
            ]
            .astype(str)
            .str.strip()
        )

    # ---------------------------------------------------------
    # Enrich metadata with MLST information
    # ---------------------------------------------------------

    if (
        metadata is not None
        and mlst is not None
    ):

        sequence_type_col = columns["sequence_type"]

        # Only merge if ST is not already present
        if sequence_type_col not in metadata.columns:

            metadata = metadata.merge(
                mlst,
                left_on=columns["sample_id"],
                right_on=columns["st_sample_id"],
                how="left",
            )

            if columns["st_sample_id"] in metadata.columns:

                metadata = metadata.drop(
                    columns=[columns["st_sample_id"]]
                )

        # Defensive cleanup if merge produced ST_x/ST_y
        x_col = f"{sequence_type_col}_x"
        y_col = f"{sequence_type_col}_y"

        if x_col in metadata.columns:

            metadata[sequence_type_col] = metadata[x_col]

            drop_cols = [
                c
                for c in [x_col, y_col]
                if c in metadata.columns
            ]

            metadata = metadata.drop(
                columns=drop_cols
            )


    # ---------------------------------------------------------
    # Normalize configured metadata fields
    # ---------------------------------------------------------

    metadata_columns_to_normalize = [
        columns["sample_id"],
        main_var,
    ]

    if var_01:
        metadata_columns_to_normalize.append(
            var_01
        )

    if var_02:
        metadata_columns_to_normalize.append(
            var_02
        )

    if columns["sequence_type"] in metadata.columns:

        metadata_columns_to_normalize.append(
            columns["sequence_type"]
        )

    for col in metadata_columns_to_normalize:

        if (
            metadata is not None
            and col in metadata.columns
        ):

            metadata[col] = (
                metadata[col]
                .astype(str)
                .str.strip()
            )


    # ---------------------------------------------------------
    # Core context hierarchy
    # ---------------------------------------------------------

    core_cluster_context = (
        analysis.get(
            "core_cluster_context",
            ["Global"]
        )
    )

    if isinstance(
        core_cluster_context,
        str,
    ):
        core_cluster_context = [
            core_cluster_context
        ]


    metadata["_core_context"] = (
        build_core_context_column(
            metadata=metadata,
            context_config=core_cluster_context,
            main_var=main_var,
            var_01=var_01,
            sequence_type_col=columns[
                "sequence_type"
            ],
        )
    )


    # ---------------------------------------------------------
    # Define grouping level
    # ---------------------------------------------------------

    if cluster_type == "Core":

        group_values = ["ALL_SAMPLES"]

    else:

        group_values = (
            metadata[facility_col]
            .dropna()
            .unique()
            .tolist()
        )


    for group_value in group_values:

        # facility_meta = metadata[
        #     metadata[facility_col] == facility
        # ].copy()

        print(f"Processing main transmission cluster group: {group_value}")

        if cluster_type == "Core":

            facility_meta = None

        else:

            facility_meta = metadata[
                metadata[facility_col] == group_value
            ].copy()

        # ---------------------------------------------------------
        # Determine metadata scope for Core Context
        # ---------------------------------------------------------

        if cluster_type == "Core":

            context_source = metadata

        else:

            context_source = facility_meta

        # ---------------------------------------------------------
        # Determine Core Context Groups
        # ---------------------------------------------------------

        if core_cluster_context == ["Global"]:

            core_context_groups = [
                "ALL_SAMPLES"
            ]

        else:

            core_context_groups = (
                context_source["_core_context"]
                .dropna()
                .unique()
                .tolist()
            )


        print(
            f"\nCore-Clustering Context Groups: "
            f"{core_context_groups}"
        )



        # ---------------------------------------------------------
        # Group-level samples
        # Used later for SNP-Epi transmission matrix construction
        # ---------------------------------------------------------


        context_cluster_parts = []
        context_k_parts = []
        context_diag_parts = []


        group_samples = [
            s
            for s in (
                context_source[
                    columns["sample_id"]
                ]
                .astype(str)
                .unique()
                .tolist()
            )
            if s in snp.index
        ]


        for context_value in core_context_groups:

            if context_value == "ALL_SAMPLES":

                context_meta = context_source.copy()

            else:

                context_meta = context_source[
                    context_source["_core_context"]
                    == context_value
                ].copy()

            facility_samples = [
                s
                for s in (
                    context_meta[
                        columns["sample_id"]
                    ]
                    .astype(str)
                    .unique()
                    .tolist()
                )
                if s in snp.index
            ]


            # Match R behaviour:
            # skip facilities with <=2 samples
            if len(facility_samples) <= 2:
                continue

            # ------------------------------------------------------------------
            # Keep TWO matrices:
            #
            # facility_snp_raw
            #     Original SNP distances
            #     Used for SNP cutoff filtering (e.g. <= 20 SNPs)
            #
            # facility_snp_scaled
            #     Scaled matrix
            #     Used for KMeans and silhouette selection
            # ------------------------------------------------------------------

            facility_snp_raw = snp.loc[
                facility_samples,
                facility_samples
            ].copy()

            scaled_values = StandardScaler().fit_transform(
                facility_snp_raw.astype(float)
            )

            facility_snp_scaled = pd.DataFrame(
                scaled_values,
                index=facility_snp_raw.index,
                columns=facility_snp_raw.columns,
            )

            if engine == "r_compatible":

                # sigN, sigNN = determine_optimal_k_r_style(facility_snp_scaled)
                sigN, sigNN, k_selection_diagnostics = (
                    determine_optimal_k_r_style(
                        facility_snp_scaled,
                        random_seed=random_seed,
                    )
                )


                plot_k_selection(
                    diagnostics_df=k_selection_diagnostics,
                    output_dir=dirs["figures"],
                    # title=f"{context_value}: Optimal number of clusters",
                    title=f"Core context {context_value}",
                    filename_prefix=f"01_k_selection_{context_value}",
                )


                core_result = run_core_snp_cluster_analysis_r_compatible(
                    mat=facility_snp_raw,
                    sigN=sigN,
                    sigNN=sigNN,
                    snpco=snp_cutoff,
                    random_seed=random_seed,
                )

                if core_result is None:
                    continue

                facility_clusters = (
                    core_result.snp_clust
                    .rename(
                        columns={
                            "name": "sample_id",
                            "cluster": "snp_cluster_id",
                            "km_cluster": "kmeans_cluster",
                        }
                    )
                    .copy()
                )

                facility_clusters["snp_cutoff"] = snp_cutoff
                facility_clusters["_core_context"] = context_value


                k_assignments = (
                    facility_clusters[
                        ["sample_id","kmeans_cluster"]
                    ]
                    .drop_duplicates()
                )

                k_diagnostics = core_result.core_result.cluster_avg_widths.copy()

            else:

                best_k, k_diagnostics = choose_optimal_k(
                    facility_snp,
                    max_k=max_k,
                    random_seed=random_seed,
                )

                k_assignments = assign_kmeans_clusters(
                    facility_snp,
                    k=best_k,
                    random_seed=random_seed,
                )

                facility_clusters = generate_snp_chains(
                    facility_snp,
                    k_assignments,
                    snp_cutoff=snp_cutoff,
                )


            if cluster_type == "Transmission":
                facility_clusters[main_var] = group_value
                k_assignments[main_var] = group_value
                k_diagnostics[main_var] = group_value


            context_diag_parts.append(k_diagnostics)
            context_k_parts.append(k_assignments)
            context_cluster_parts.append(facility_clusters)



        if len(context_cluster_parts):

            facility_clusters = pd.concat(
                context_cluster_parts,
                ignore_index=True,
            )

            k_assignments = pd.concat(
                context_k_parts,
                ignore_index=True,
            )

            k_diagnostics = pd.concat(
                context_diag_parts,
                ignore_index=True,
            )

        else:

            continue


        all_k_diagnostics.append(k_diagnostics)
        all_k_assignments.append(k_assignments)
        all_snp_clusters.append(facility_clusters)


        if str(
            analysis.get(
                "cluster_type",
                "Core"
            )
        ).lower() == "transmission":


            # -----------------------------------------
            # R-compatible SNP-Epi clustering
            # -----------------------------------------


            epiwkDF = facility_meta.copy()

            mdf = matrix_to_mdf(
                snp.loc[
                    group_samples,
                    group_samples,
                ]
            )


            snpClust = (
                facility_clusters.rename(
                    columns={
                        "sample_id": "name",
                        "snp_cluster_id": "cluster",
                    }
                )[
                    ["name", "cluster"]
                ]
                .copy()
            )

            transmission = calculate_snp_epi_clusters_r_compatible(
                snpClust=snpClust[
                    ["name", "cluster"]
                ].copy(),
                epiwkDF=epiwkDF,
                mdf=mdf,
                sample_id_col=columns["sample_id"],
                date_col=var_02,
                st_col=columns["sequence_type"],
                date_format=analysis.get(
                    "date_format",
                    "ymd",
                ),
                daysco=int(
                    analysis.get(
                        "days_cutoff",
                        45,
                    )
                ),
                excl_vec=[
                    "Date2",
                    "epicumsum",
                    "CG",
                    "num",
                    "name",
                    "cluster",
                    "km_cluster",
                ],
            )

            if transmission is not None and len(transmission):

                # transmission[main_var] = facility
                transmission[main_var] = group_value
                all_transmission_clusters.append(
                    transmission
                )

    k_diagnostics = pd.concat(
        all_k_diagnostics,
        ignore_index=True,
    )

    k_assignments = pd.concat(
        all_k_assignments,
        ignore_index=True,
    )

    snp_clusters = pd.concat(
        all_snp_clusters,
        ignore_index=True,
    )

    # Canonical context-aware Core SNP cluster identifier.
    required_core_columns = {"_core_context", "snp_cluster_id"}
    missing_core_columns = required_core_columns - set(snp_clusters.columns)
    if missing_core_columns:
        raise ValueError(
            "Cannot construct core_snp_cluster. Missing columns: "
            f"{sorted(missing_core_columns)}"
        )

    snp_clusters["core_snp_cluster"] = pd.NA
    assigned_core_cluster = (
        snp_clusters["_core_context"].notna()
        & snp_clusters["snp_cluster_id"].notna()
    )
    cluster_numbers = (
        pd.to_numeric(
            snp_clusters.loc[assigned_core_cluster, "snp_cluster_id"],
            errors="raise",
        )
        .astype("Int64")
        .astype(str)
    )
    snp_clusters.loc[
        assigned_core_cluster,
        "core_snp_cluster",
    ] = (
        snp_clusters.loc[assigned_core_cluster, "_core_context"]
        .astype(str)
        .str.strip()
        + "__"
        + cluster_numbers
    )

    transmission = (
        pd.concat(
            all_transmission_clusters,
            ignore_index=True
        )
        if all_transmission_clusters
        else None
    )

    write_csv(
        k_diagnostics,
        dirs["tables"] / "k_selection_diagnostics.csv",
    )

    write_csv(
        k_assignments,
        dirs["tables"] / "kmeans_cluster_assignments.csv",
    )

    # Complete Core annotation population. Start from metadata so samples
    # without a retained Core cluster remain available to all outputs.
    core_metadata = (
        metadata.drop_duplicates(
            subset=[columns["sample_id"]],
            keep="first",
        )
        .copy()
    )
    core_metadata[columns["sample_id"]] = (
        core_metadata[columns["sample_id"]].astype(str).str.strip()
    )

    core_assignments = (
        snp_clusters[
            ["sample_id", "snp_cluster_id", "core_snp_cluster"]
        ]
        .drop_duplicates(subset=["sample_id"], keep="first")
        .copy()
    )
    core_assignments["sample_id"] = (
        core_assignments["sample_id"].astype(str).str.strip()
    )

    core_annotations = core_metadata.merge(
        core_assignments,
        left_on=columns["sample_id"],
        right_on="sample_id",
        how="left",
        validate="one_to_one",
    )
    if columns["sample_id"] != "sample_id":
        core_annotations = core_annotations.drop(columns=["sample_id"])

    # Compatibility alias for components that still use the old name.
    core_annotations["core_cluster_label"] = (
        core_annotations["core_snp_cluster"]
    )

    core_output_columns = [
        columns["sample_id"],
        "_core_context",
        "snp_cluster_id",
        "core_snp_cluster",
    ]
    for column in (
        columns.get("sequence_type"),
        main_var,
        var_01,
        var_02,
    ):
        if (
            column
            and column in core_annotations.columns
            and column not in core_output_columns
        ):
            core_output_columns.append(column)

    write_csv(
        core_annotations.loc[:, core_output_columns].copy(),
        dirs["tables"] / "core_snp_clusters.csv",
    )

    if transmission is not None:
        transmission_clusters = transmission.merge(
            metadata,
            left_on="sample_id",
            right_on=columns["sample_id"],
            how="left",
        )
        if columns["sample_id"] in transmission_clusters.columns:
            transmission_clusters = transmission_clusters.drop(
                columns=[columns["sample_id"]]
            )
        duplicate_columns = [
            column
            for column in transmission_clusters.columns
            if column.endswith("_x")
        ]
        if duplicate_columns:
            transmission_clusters = transmission_clusters.drop(
                columns=duplicate_columns
            )
        transmission_clusters.columns = [
            column.replace("_y", "")
            for column in transmission_clusters.columns
        ]

        transmission_assignments = (
            transmission[
                [
                    "sample_id",
                    "Clusters",
                    "Days",
                    "SNPs",
                    "Cluster_Cases_count",
                ]
            ]
            .drop_duplicates(subset=["sample_id"], keep="first")
            .copy()
        )
        transmission_assignments["sample_id"] = (
            transmission_assignments["sample_id"].astype(str).str.strip()
        )
        transmission_context_df = core_annotations.merge(
            transmission_assignments,
            left_on=columns["sample_id"],
            right_on="sample_id",
            how="left",
            validate="one_to_one",
        )
        if columns["sample_id"] != "sample_id":
            transmission_context_df = transmission_context_df.drop(
                columns=["sample_id"]
            )

        write_csv(
            transmission_clusters,
            dirs["tables"] / "transmission_clusters.csv",
        )
        write_csv(
            transmission,
            dirs["tables"] / "isolate_cluster_assignments.csv",
        )
        write_csv(
            transmission_context_df,
            dirs["tables"] / "transmission_context.csv",
        )
    else:
        write_csv(
            snp_clusters,
            dirs["tables"] / "isolate_cluster_assignments.csv",
        )

    if config.outputs.get("generate_figures", True):


        # ---------------------------------------------------------
        # Transmission visualizations
        # ---------------------------------------------------------

        if transmission is not None:

            print(
                "Generating transmission scatter plots "
                f"({len(transmission_context_df):,} analysed isolates)"
            )

            transmission_level = str(
                analysis.get(
                    "transmission_level",
                    "Facility",
                )
            ).strip().casefold()

            # -----------------------------------------------------
            # Facility / Area modes
            #
            # One plot per main_var value.
            #
            # Facility:
            #   main_var = FacilityName
            #   y-axis = ST
            #   shape = var_01, normally WardType
            #
            # Area:
            #   main_var = Area
            #   y-axis = ST
            #   shape = var_01, normally FacilityName
            # -----------------------------------------------------

            if transmission_level in {
                "facility",
                "area",
            }:

                for group_value, group_df in (
                    transmission_context_df.groupby(
                        main_var,
                        dropna=False,
                    )
                ):

                    if pd.isna(group_value):
                        continue

                    safe_name = (
                        str(group_value)
                        .strip()
                        .replace(" ", "_")
                        .replace("/", "_")
                        .replace("\\", "_")
                        .replace(":", "_")
                    )

                    print(
                        "Generating transmission scatter plot: "
                        f"{group_value}"
                    )

                    generate_transmission_scatter_from_config(
                        transmission_df=group_df,
                        config=config,
                        output_dir=dirs["figures"],
                        cluster_column="Clusters",
                        filename_prefix=(
                            f"02_transmission_scatter_{safe_name}"
                        ),
                        title=str(group_value),
                    )

            # -----------------------------------------------------
            # Community mode
            #
            # main_var must be ST.
            # One plot is produced per ST.
            # y-axis = var_01, normally FacilityName.
            # Shape encoding is disabled by the config wrapper.

            # if y_axis_label is None:
            #     y_axis_label = "Facility"
            # -----------------------------------------------------

            elif transmission_level == "community":

                for st_value, st_df in (
                    transmission_context_df.groupby(
                        main_var,
                        dropna=False,
                    )
                ):

                    if pd.isna(st_value):
                        continue

                    safe_st = (
                        str(st_value)
                        .strip()
                        .replace(" ", "_")
                        .replace("/", "_")
                        .replace("\\", "_")
                        .replace(":", "_")
                    )

                    print(
                        "Generating community transmission "
                        f"scatter plot: ST {st_value}"
                    )

                    generate_transmission_scatter_from_config(
                        transmission_df=st_df,
                        config=config,
                        output_dir=dirs["figures"],
                        cluster_column="Clusters",
                        filename_prefix=(
                            f"02_transmission_scatter_ST_{safe_st}"
                        ),
                        title=f"ST {st_value}",
                    )

            else:

                raise ValueError(
                    "Unsupported transmission_level: "
                    f"{analysis.get('transmission_level')}. "
                    "Expected Facility, Area, or Community."
                )

        # -----------------------------------------------------
        # Core SNP visualizations
        # -----------------------------------------------------
        heatmap_annotations = (
            transmission_context_df.copy()
            if transmission is not None
            else core_annotations.copy()
        )

        print(
            "Generating annotated Core SNP heatmap "
            f"({len(heatmap_annotations):,} isolates; "
            f"{heatmap_annotations['core_snp_cluster'].notna().sum():,} "
            "Core-clustered isolates)"
        )
        generate_annotated_heatmap_from_config(
            snp_matrix=snp,
            annotation_df=heatmap_annotations,
            config=config,
            output_dir=dirs["figures"],
            core_cluster_column="core_snp_cluster",
            transmission_cluster_column="Clusters",
            filename_prefix="03_core_heatmap",
            title="Core SNP distance heatmap",
        )

        # ---------------------------------------------------------
        # Core SNP minimum spanning trees
        # ---------------------------------------------------------
        if main_var not in heatmap_annotations.columns:
            raise ValueError(
                f"Cannot generate Core MSTs: {main_var} is missing."
            )

        mst_sample_id_column = columns["sample_id"]
        mst_annotations = heatmap_annotations.copy()
        mst_annotations[mst_sample_id_column] = (
            mst_annotations[mst_sample_id_column].astype(str).str.strip()
        )
        snp_sample_ids = (
            set(snp.index.astype(str))
            & set(snp.columns.astype(str))
        )
        mst_annotations = (
            mst_annotations[
                mst_annotations[mst_sample_id_column].isin(snp_sample_ids)
            ]
            .drop_duplicates(subset=[mst_sample_id_column], keep="first")
            .copy()
        )
        if mst_annotations.empty:
            raise ValueError(
                "Cannot generate Core MSTs: no annotation IDs match the SNP matrix."
            )

        is_global_core_context = (
            len(core_cluster_context) == 1
            and str(core_cluster_context[0]).strip().casefold() == "global"
        )
        sequence_type_column = columns.get("sequence_type")
        all_mst_edges = []

        def generate_one_core_mst(
            group_value: str,
            group_annotations: pd.DataFrame,
            output_filename: str,
            split_label: str,
            shape_column: str | None,
        ) -> None:
            group_sample_ids = (
                group_annotations[mst_sample_id_column]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )
            group_sample_ids = [
                sample_id
                for sample_id in group_sample_ids
                if sample_id in snp.index and sample_id in snp.columns
            ]
            if len(group_sample_ids) < 2:
                print(
                    f"Skipping Core MST for {split_label}={group_value}: "
                    "fewer than two eligible isolates"
                )
                return

            group_tree = build_minimum_spanning_tree(
                snp.loc[group_sample_ids, group_sample_ids].copy()
            )
            selected_annotations = group_annotations[
                group_annotations[mst_sample_id_column]
                .astype(str)
                .isin(group_sample_ids)
            ].copy()
            selected_annotations[mst_sample_id_column] = (
                selected_annotations[mst_sample_id_column].astype(str)
            )

            context_map = (
                selected_annotations[
                    [mst_sample_id_column, "_core_context"]
                ]
                .drop_duplicates(subset=[mst_sample_id_column])
                .set_index(mst_sample_id_column)["_core_context"]
                .to_dict()
            )
            group_edges = mst_edges_table(group_tree)
            group_edges["mst_scope"] = split_label
            group_edges["mst_scope_value"] = str(group_value)

            def edge_context(row: pd.Series) -> object:
                context_1 = context_map.get(str(row["sample_id_1"]))
                context_2 = context_map.get(str(row["sample_id_2"]))
                if context_1 == context_2:
                    return context_1
                return "BETWEEN_CORE_CONTEXTS"

            group_edges["_core_context"] = group_edges.apply(
                edge_context,
                axis=1,
            )
            all_mst_edges.append(
                group_edges[
                    [
                        "mst_scope",
                        "mst_scope_value",
                        "_core_context",
                        "sample_id_1",
                        "sample_id_2",
                        "snp_distance",
                    ]
                ]
            )

            print(
                "Generating Core SNP minimum spanning tree: "
                f"{split_label}={group_value} "
                f"({group_tree.number_of_nodes():,} isolates)"
            )
            generate_core_mst_from_config(
                tree=group_tree,
                annotation_df=selected_annotations,
                config=config,
                output_path=dirs["figures"] / output_filename,
                split_value=str(group_value),
                split_label=split_label,
                shape_column=shape_column,
                core_cluster_column="core_snp_cluster",
                core_context_column="_core_context",
                title="Core SNP minimum spanning tree",
            )

        if is_global_core_context:
            generate_one_core_mst(
                group_value="ALL_SAMPLES",
                group_annotations=mst_annotations,
                output_filename="04_mst_ALL_SAMPLES.html",
                split_label="Core context",
                shape_column=main_var,
            )

        for main_value, main_df in mst_annotations.groupby(
            main_var,
            dropna=False,
            sort=False,
        ):
            if pd.isna(main_value):
                continue
            safe_main_value = str(main_value).strip()
            for character in (" ", "/", "\\", ":"):
                safe_main_value = safe_main_value.replace(character, "_")
            safe_main_value = safe_main_value or "UNKNOWN"

            shape_column = (
                var_01
                if (
                    sequence_type_column
                    and str(main_var) == str(sequence_type_column)
                )
                else sequence_type_column
            )
            generate_one_core_mst(
                group_value=str(main_value),
                group_annotations=main_df,
                output_filename=f"04_mst_{safe_main_value}.html",
                split_label=str(main_var),
                shape_column=str(shape_column) if shape_column else None,
            )

        mst_edges = (
            pd.concat(all_mst_edges, ignore_index=True)
            if all_mst_edges
            else pd.DataFrame(
                columns=[
                    "mst_scope",
                    "mst_scope_value",
                    "_core_context",
                    "sample_id_1",
                    "sample_id_2",
                    "snp_distance",
                ]
            )
        )
        write_csv(
            mst_edges,
            dirs["tables"] / "minimum_spanning_tree_edges.csv",
        )

        # ---------------------------------------------------------
        # Transmission-cluster networks
        # ---------------------------------------------------------

        if transmission is not None:

            print(
                "Generating transmission-cluster networks "
                f"({len(transmission_context_df):,} analysed isolates)"
            )

            transmission_network_outputs = (
                generate_transmission_networks_from_config(
                    transmission_df=transmission_context_df,
                    snp_matrix=snp,
                    config=config,
                    output_dir=dirs["figures"],
                    cluster_column="Clusters",
                    filename_prefix=(
                        "05_transmission_network"
                    ),
                )
            )

            transmission_network_nodes = (
                transmission_network_outputs[
                    "nodes"
                ]
            )

            transmission_network_edges = (
                transmission_network_outputs[
                    "edges"
                ]
            )

            if not transmission_network_nodes.empty:

                write_csv(
                    transmission_network_nodes,
                    dirs["tables"]
                    / "transmission_network_nodes.csv",
                )

            else:

                print(
                    "No clustered isolates were available "
                    "for transmission-network node output"
                )

            if not transmission_network_edges.empty:

                write_csv(
                    transmission_network_edges,
                    dirs["tables"]
                    / "transmission_network_edges.csv",
                )

            else:

                print(
                    "No transmission-network edges were generated"
                )

    # ---------------------------------------------------------
    # Open scientific HTML report
    # ---------------------------------------------------------

    if config.outputs.get(
        "generate_html_report",
        True,
    ):

        report_transmission = (
            transmission_context_df
            if transmission is not None
            else None
        )

        write_html_report(
            dirs["report"]
            / "snp2cluster_report.html",
            project_name=config.project_name,
            assignments=core_annotations,
            transmission=report_transmission,
            config=config,
            figures_dir=dirs["figures"],
            tables_dir=dirs["tables"],
            # software_version=__version__, # NameError: name '__version__' is not defined
            # run_manifest=run_manifest, # NameError: name 'run_manifest' is not defined
        )


    return dirs
