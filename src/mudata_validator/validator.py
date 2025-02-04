import os
import warnings

import muon as mu
import numpy as np
import pandas as pd
import scipy.sparse


def check_duplicate_objects(
    data: pd.DataFrame, error_messages: list, modality_name: str = None
):
    """Check for duplicate object IDs in the data."""
    if len(set(data.index)) == data.shape[0]:
        return
    counts = data.index.value_counts()
    duplicates = counts[counts > 1]
    message_pieces = [
        f"Found duplicate object IDs in {modality_name if modality_name else 'data'}:",
        *(f"\t{i}\t({count} occurrences)" for i, count in duplicates.items()),
    ]
    error_messages.append("\n".join(message_pieces))
    warnings.warn(
        f"If this data is from multiple datasets, you must prepend the barcode with the plate/well, HuBMAP ID, or the HuBMAP UUID."
    )


def check_sparsity(matrix, matrix_name: str):
    """Check the sparsity of a matrix and warn if it's too dense."""
    if isinstance(matrix, np.ndarray):
        sparsity = scipy.sparse.csr_matrix(matrix).nnz / np.prod(matrix.shape)
        if sparsity > 0.3:
            warnings.warn(
                f"{matrix_name} is a dense matrix with sparsity {sparsity:.4f}. It is recommended to store this as a sparse matrix.",
                UserWarning,
            )


def validate_modality(adata, modality_name, error_messages):
    """Validate a single modality (AnnData object)."""
    print(f"Validating modality: {modality_name}")

    # REQUIRED: Check for duplicate values in the index
    # Change to obs_names?
    print(
        "The values in AnnData.obs.index will be used as the objects' unique identifiers. They look like:"
    )
    print(adata.obs.head().index)
    check_duplicate_objects(adata.obs, error_messages, modality_name)

    # Validate `.obs` fields
    if "original_obs_id" not in adata.obs.columns:
        error_messages.append(
            f"`{modality_name}.obs` must contain a column named 'original_obs_id' containing the original barcode or unique identifier."
        )

    if "object_type" in adata.obs.columns:
        allowed_obj_types = ["cell", "nuclei", "ftu", "spot"]
        invalid_values = set(adata.obs["object_type"].unique()) - set(allowed_obj_types)
        if invalid_values:
            error_messages.append(
                f"'{modality_name}.obs['object_type']' contains invalid values: {', '.join(invalid_values)}. "
                f"Allowed values are: {', '.join(allowed_obj_types)}."
            )
    else:
        error_messages.append(
            f"`{modality_name}.obs` must contain a column named 'object_type' containing the observation type ontology ID (cell, nuclei, ftu, spot)."
        )

    # !!TODO!! Check var values
    print(
        "The HUGO symbol should be included as an annotation for genes and the Uniprot ID should be included as an annotation for proteins."
    )

    if "protocol" not in adata.uns_keys() or adata.uns["protocol"] == None:
        error_messages.append(
            f"`{modality_name}.uns` must contain a key 'protocol' with a valid Protocol DOI."
        )
    valid_analyte_classes = [
        "DNA",
        "RNA",
        "Endogenous fluorophore",
        "Lipid",
        "Metabolite",
        "Polysaccharide",
        "Protein",
        "Nucleic acid + protein",
        "N-glycan",
        "DNA + RNA",
        "Chromatin",
        "Collagen",
        "Fluorochrome",
        "Lipid + metabolite",
        "Peptide",
        "Saturated lipid",
        "Unsaturated lipid",
    ]
    if "analyte_class" not in adata.uns_keys():
        error_messages.append(
            ".uns must contain a key called 'analyte_class' that references a known analyte class defined in 'valid_analyte_classes.txt'."
        )
    elif adata.uns.get("analyte_class") not in valid_analyte_classes:
        error_messages.append(
            "The value in .uns['analyte_class'] must reference a known analyte class defined in 'valid_analyte_classes.txt'."
        )

    # Check sparsity for all matrices

    for layer, key_set in [
        (adata.layers, set()),
        (adata.obsm, set()),
        (adata.obsp, set()),
        (adata.varm, set()),
        (adata.varp, set()),
    ]:
        if hasattr(layer, "keys"):
            for key in layer.keys():
                key_set.add(key)
                check_sparsity(layer[key], f"{modality_name}[{key}]")

    print(
        "Standard plots are expected to be stored in .obsm['X_umap'], .obsm['X_harmony'], .obsm['X_tsne'] and .obsm['X_pca']"
    )
    print("If this is spatial data, coordinates should go in .obsm['X_spatial']")


def validate_annotations(adata, modality_name, error_messages):
    if "annotation" in adata.obsm:
        if "annotation_methods" not in adata.uns:
            error_messages.append(
                f"`{modality_name}.obsm['annotation']` exists, but `{modality_name}.uns['annotation_methods']` is missing."
            )
    else:
        warnings.warn(
            f"It is recommended to use `{modality_name}.obsm['annotation']` for general annotation storage.",
            UserWarning,
        )


def validate_analyses(adata, modality_name, error_messages):
    check_sparsity(adata.X, f"{modality_name}.X")


def validate_mudata(input_data):
    """
    Validates a MuData object or an H5mu file.

    Parameters:
    - input_data: str or muon.MuData
      Either a path to an H5mu file or a MuData object.

    Raises:
    - ValueError: If validation fails with error messages.
    - Warnings for non-critical issues.

    Returns:
    - None: Prints success if validation passes.
    """
    error_messages = []

    if isinstance(input_data, mu.MuData):
        mdata = input_data
    else:
        mdata = mu.read_h5mu(input_data)

    print("Validating overall MuData object...")

    if mdata.uns.get("epic_type") == {"annotations"}:
        for modality_name, adata in mdata.mod.items():
            validate_annotations(adata, modality_name, error_messages)
            validate_modality(adata, modality_name, error_messages)
    elif mdata.uns.get("epic_type") == {"analyses"}:
        for modality_name, adata in mdata.mod.items():
            validate_analyses(adata, modality_name, error_messages)
            validate_modality(adata, modality_name, error_messages)
    elif mdata.uns.get("epic_type") == {"annotations", "analyses"}:
        for modality_name, adata in mdata.mod.items():
            validate_analyses(adata, modality_name, error_messages)
            validate_annotations(adata, modality_name, error_messages)
            validate_modality(adata, modality_name, error_messages)
    else:
        error_messages.append(
            "MuData.uns must contain a key called 'epic_type' with at least one valid epic type: annotations, analyses"
        )

    # Raise an error if validation fails
    if error_messages:
        formatted_errors = "\n- ".join(error_messages)
        raise ValueError(
            f"Validation failed with the following issues:\n- {formatted_errors}"
        )

    print("Validation passed!")
