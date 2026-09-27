from __future__ import annotations

from importlib import metadata
from pathlib import Path
from pyexpat.errors import messages
from typing import Any
import pandas as pd

from .config import SNP2ClusterConfig
from .io import read_table, read_snp_matrix


class ValidationError(Exception):
    """Raised when input validation fails."""


def validate_config_structure(config: SNP2ClusterConfig) -> list[str]:
    messages: list[str] = []
    required_sections = ["project", "inputs", "columns", "variables", "analysis", "outputs"]
    for section in required_sections:
        if section not in config.raw:
            messages.append(f"Missing config section: {section}")

    required_inputs = ["metadata", "snp_matrix"]
    for key in required_inputs:
        if key not in config.inputs:
            messages.append(f"Missing required input path: inputs.{key}")

    required_variables = ["main_var"] #, "var_01", "var_02"]
    for key in required_variables:
        if key not in config.variables:
            messages.append(f"Missing required variable: variables.{key}")

    required_columns = ["sample_id"] #, "collection_date"]
    for key in required_columns:
        if key not in config.columns:
            messages.append(f"Missing required column mapping: columns.{key}")

    return messages


def validate_input_files_exist(config: SNP2ClusterConfig) -> list[str]:
    messages: list[str] = []
    for label, path_value in config.inputs.items():
        if path_value and not Path(path_value).exists():
            messages.append(f"Input file not found for {label}: {path_value}")
    return messages



def validate_loaded_inputs(config: SNP2ClusterConfig) -> list[str]:
    messages: list[str] = []
    
    metadata_path = config.inputs.get("metadata")
    snp_path = config.inputs.get("snp_matrix")
    mlst_path = config.inputs.get("mlst_profile")

    # --- 1. Validate SNP Matrix (Required for BOTH modes) ---
    if not snp_path:
        messages.append("Configuration missing: 'snp_matrix' is required.")
        return messages 
    
    snp = read_snp_matrix(snp_path)

    # --- 2. Extract Cluster Type ---
    cluster_type = config.analysis.get("cluster_type", "Core") 
    if cluster_type not in {"Core", "Transmission"}:
        messages.append(f"Unsupported cluster_type: {cluster_type}")
        return messages
    
    # --- 2b. Validate date format configuration ---

    supported_date_formats = {
        "ymd",
        "ydm",
        "dmy",
        "mdy",
    }

    date_format = config.analysis.get(
        "date_format",
        "ymd",
    )

    if (
            cluster_type == "Transmission"
            and date_format not in supported_date_formats
        ):

        messages.append(
            f"Unsupported date_format: '{date_format}'. "
            f"Supported values: "
            f"{sorted(supported_date_formats)}"
        )

    # --- 3. Mode-Specific File Requirements & Loading ---
    metadata_df = None
    mlst_df = None

    if cluster_type == "Transmission":
        missing_transmission_files = False
        
        if not metadata_path:
            messages.append("Configuration missing: 'metadata' is mandatory for Transmission analysis.")
            missing_transmission_files = True
        elif Path(metadata_path).exists():
            metadata_df = read_table(metadata_path)
            
        if not mlst_path:
            messages.append("Configuration missing: 'mlst_profile' is mandatory for Transmission analysis.")
            missing_transmission_files = True
        elif Path(mlst_path).exists():
            mlst_df = read_table(mlst_path)
            
        if missing_transmission_files:
            return messages

    elif cluster_type == "Core":
        if metadata_path and Path(metadata_path).exists():
            metadata_df = read_table(metadata_path)
        if mlst_path and Path(mlst_path).exists():
            mlst_df = read_table(mlst_path)

    # --- 4. Extract Configurations ---
    sample_col = config.columns.get("sample_id")
    seq_col = config.columns.get("sequence_type")
    mlst_sample_col = config.columns.get("st_sample_id")
    
    main_var = config.variables.get("main_var")
    var_01 = config.variables.get("var_01")
    var_02 = config.variables.get("var_02")

    core_cluster_context = config.analysis.get(
        "core_cluster_context",
        ["Global"],
    )

    if isinstance(core_cluster_context, str):
        core_cluster_context = [
            core_cluster_context
        ]

    # ---------------------------------------------------------
    # Core Context Validation
    # ---------------------------------------------------------

    allowed_contexts = {
        "Global",
        "MainVar",
        "SecondaryVar",
        "ST",
    }

    invalid_contexts = [
        c
        for c in core_cluster_context
        if c not in allowed_contexts
    ]

    if invalid_contexts:

        messages.append(
            "Unsupported core context level(s): "
            f"{invalid_contexts}. "
            f"Supported values are: "
            f"{sorted(allowed_contexts)}"
        )
    
    if len(core_cluster_context) != len(
        set(core_cluster_context)
    ):
        messages.append(
            "Duplicate core context levels are not allowed."
        )


    if (
        "Global" in core_cluster_context
        and len(core_cluster_context) > 1
    ):
        messages.append(
            "Global cannot be combined with "
            "other core context levels."
        )

    if "MainVar" in core_cluster_context:

        if not main_var:

            messages.append(
                "MainVar core context requested "
                "but variables.main_var "
                "is not configured."
            )

    if "SecondaryVar" in core_cluster_context:

        if not var_01:

            messages.append(
                "SecondaryVar core context requested "
                "but variables.var_01 "
                "is not configured."
            )

    if "ST" in core_cluster_context:

        if not seq_col:

            messages.append(
                "ST core context requested "
                "but columns.sequence_type "
                "is not configured."
            )

    contexts_requiring_metadata = {
        "MainVar",
        "SecondaryVar",
    }

    if any(
        c in contexts_requiring_metadata
        for c in core_cluster_context
    ):

        if metadata_df is None:

            messages.append(
                "Selected core_cluster_context "
                f"{core_cluster_context} "
                "requires metadata."
            )

    if "ST" in core_cluster_context:

        st_in_metadata = (
            metadata_df is not None
            and seq_col
            and seq_col in metadata_df.columns
        )

        st_in_mlst = (
            mlst_df is not None
            and seq_col
            and seq_col in mlst_df.columns
        )

        if not (st_in_metadata or st_in_mlst):

            messages.append(
                "ST core context requested but "
                "sequence type information "
                "is not available in metadata "
                "or MLST profile."
            )


    # --- 5. YAML Configuration Validation ---
    # if not main_var:
    #     messages.append("Configuration missing: 'main_var' is required in the variables block.")
    if cluster_type == "Transmission":

        if not main_var:
            messages.append(
                "Configuration missing: 'main_var' is required "
                "for Transmission analysis."
            )

    if cluster_type == "Transmission" and not var_02:
        messages.append("Configuration missing: 'var_02' (collection date) is mandatory for Transmission analysis.")

    # --- 6. Metadata CSV Column Validation ---
    if metadata_df is not None:
        if not sample_col:
            messages.append("Configuration missing: 'sample_id' is required when passing metadata.")
        elif sample_col not in metadata_df.columns: 
            messages.append(f"Metadata sample ID column not found: {sample_col}")

        # if main_var and main_var not in metadata_df.columns:
        #     messages.append(f"Metadata main_var column not found: {main_var}")

        if main_var:

            main_var_in_metadata = (
                main_var in metadata_df.columns
            )

            main_var_is_sequence_type = (
                seq_col is not None
                and main_var == seq_col
            )

            main_var_in_mlst = (
                main_var_is_sequence_type
                and mlst_df is not None
                and seq_col in mlst_df.columns
            )

            if not (
                main_var_in_metadata
                or main_var_in_mlst
            ):
                messages.append(
                    "Configured main_var column not found "
                    "in metadata or available through "
                    f"MLST enrichment: {main_var}"
                )

        if var_01 and var_01 not in metadata_df.columns:
            messages.append(f"Metadata var_01 column not found: {var_01}")
            
        if var_02 and var_02 not in metadata_df.columns:
            messages.append(f"Metadata var_02 (collection date) column not found: {var_02}")

    # --- 7. MLST Column Validation ---
    if mlst_df is not None:
        # Check if the user defined the MLST sample ID column in YAML
        if not mlst_sample_col:
            messages.append("Configuration missing: 'st_sample_id' is required when providing an MLST profile.")
        elif mlst_sample_col not in mlst_df.columns:
            messages.append(f"MLST sample ID column not found in file: {mlst_sample_col}")
            
        # Check if the user defined the Sequence Type column in YAML
        if not seq_col:
            messages.append("Configuration missing: 'sequence_type' is required when providing an MLST profile.")
        else:
            in_mlst = seq_col in mlst_df.columns
            in_meta = (metadata_df is not None) and (seq_col in metadata_df.columns)
            
            if not (in_mlst or in_meta):
                messages.append(f"Sequence type column '{seq_col}' not found in MLST file or metadata.")

    # --- 8. SNP Matrix Validation ---
    if snp.shape[0] != snp.shape[1]:
        messages.append("SNP matrix must be square.")

    if set(snp.index) != set(snp.columns):
        messages.append("SNP matrix row IDs and column IDs do not match exactly.")

    try:
        snp.astype(float)
    except Exception:
        messages.append("SNP matrix contains non-numeric values outside expected sample labels.")

    # --- 9. Cross-Validation (Metadata vs SNP) ---
    if metadata_df is not None and sample_col and sample_col in metadata_df.columns:
        metadata_ids = set(metadata_df[sample_col].astype(str))
        snp_ids = set(snp.index.astype(str))

        ignored_samples = {"reference"}

        missing_in_snp = sorted(metadata_ids - snp_ids)
        missing_in_meta = sorted(snp_ids - metadata_ids - ignored_samples)

        if missing_in_snp:
            messages.append(f"Metadata samples missing from SNP matrix: {missing_in_snp[:10]}")

        if missing_in_meta:
            messages.append(f"SNP matrix samples missing from metadata: {missing_in_meta[:10]}")

    return messages

def validate_all(config: SNP2ClusterConfig, strict_files: bool = True) -> list[str]:
    messages = []
    messages.extend(validate_config_structure(config))
    if strict_files:
        messages.extend(validate_input_files_exist(config))
        if not messages:
            messages.extend(validate_loaded_inputs(config))
    return messages
