from __future__ import annotations

import pandas as pd


def build_transmission_clusters(
    snp_clusters: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_col: str,
    date_col: str,
    main_var: str | None,
    sequence_type_col: str | None,
    days_cutoff: int,
) -> pd.DataFrame:
    """Generate transmission clusters from SNP clusters + epidemiological metadata.

    Conservative starter rule:
    - join SNP clusters to metadata
    - keep analyses within sequence type if available
    - split by facility/main variable if available
    - sort by collection date
    - assign transmission cluster if cases in same SNP cluster/group are within the configured window

    IMPORTANT: Exact R business rules should be ported and validated from the original functions before
    production use.
    """
    meta = metadata.copy()
    meta[sample_col] = meta[sample_col].astype(str)
    meta[date_col] = pd.to_datetime(meta[date_col], errors="coerce")

    df = snp_clusters.merge(meta, left_on="sample_id", right_on=sample_col, how="left")

    group_cols = ["snp_cluster_id"]
    if sequence_type_col and sequence_type_col in df.columns:
        group_cols.append(sequence_type_col)
    if main_var and main_var in df.columns:
        group_cols.append(main_var)

    rows = []
    for group_key, sub in df.sort_values(date_col).groupby(group_cols, dropna=False):
        sub = sub.copy().sort_values(date_col)
        current_index = 1
        cluster_start_date = None

        for _, row in sub.iterrows():
            current_date = row[date_col]
            if pd.isna(current_date):
                transmission_cluster = f"TC_UNDATED_{current_index}"
            else:
                if cluster_start_date is None:
                    cluster_start_date = current_date
                elif (current_date - cluster_start_date).days > days_cutoff:
                    current_index += 1
                    cluster_start_date = current_date
                transmission_cluster = f"TC{current_index}"

            out = row.to_dict()
            out["transmission_cluster_index"] = current_index
            out["transmission_cluster_id"] = f"{row['snp_cluster_id']}_{transmission_cluster}"
            out["days_cutoff"] = days_cutoff
            out["cluster_rule_applied"] = "same_snp_cluster_same_ST_main_var_within_days_cutoff"
            rows.append(out)

    return pd.DataFrame(rows)
