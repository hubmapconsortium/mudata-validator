import os
import warnings

import anndata as ad
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
        if matrix.ndim != 2:
            return
        if not np.issubdtype(matrix.dtype, np.number):
            return
        sparsity = scipy.sparse.csr_matrix(matrix).nnz / np.prod(matrix.shape)
        if sparsity > 0.3:
            warnings.warn(
                f"{matrix_name} is a dense matrix with sparsity {sparsity:.4f}. It is recommended to store this as a sparse matrix.",
                UserWarning,
            )


def validate_obsm_x_spatial(
    matrix: np.ndarray, modality_name: str, error_messages: list
):
    if matrix.ndim < 2:
        error_messages.append(
            f"Only one dimension found in {modality_name}.obsm['X_spatial']; Expecting at least (y,x)."
        )


def validate_modality(adata: ad.AnnData, modality_name: str, error_messages: list):
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

    # !!TODO!! Check var values?
    print(
        "The HUGO symbol should be included as an annotation for genes and the Uniprot ID should be included as an annotation for proteins."
    )

    if "protocol" not in adata.uns_keys() or adata.uns["protocol"] == None:
        error_messages.append(
            f"`{modality_name}.uns` must contain a key 'protocol' with a valid Protocol DOI."
        )

    # Check analyte class
    try:

        common_assay_fields = pd.read_csv(
            "https://docs.google.com/spreadsheets/d/1oKBb0Elie4wNzjvqqQEVpH7I4-8phKyPn2oyvG3rNv4/export?gid=0&format=csv",
            header=1,
        )
        valid_analyte_classes = common_assay_fields["analyte class"].dropna().to_list()
    except Exception as e:
        print(
            f"Error fetching Common Assay Fields data: {e}. Falling back to stored analyte classes."
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
            f"The value in `{modality_name}.uns['analyte_class']` must reference a known analyte class. \n"
            "The known analyte classes are: \n"
            f"{valid_analyte_classes}"
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

    if "X_spatial" in adata.obsm_keys():
        validate_obsm_x_spatial(adata.obsm["X_spatial"], modality_name, error_messages)


def validate_annotations(adata, modality_name):
    if "annotation" in adata.obsm:
        pass
    else:
        warnings.warn(
            f"It is recommended to use `{modality_name}.obsm['annotation']` for general annotation storage.",
            UserWarning,
        )


def validate_analyses(adata, modality_name):
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

    # Normalize epic_type to a set for reliable comparison
    epic_type = mdata.uns.get("epic_type")
    if isinstance(epic_type, (list, set, np.ndarray)):
        epic_type = set(epic_type)

    if epic_type == {"annotations"}:
        for modality_name, adata in mdata.mod.items():
            validate_annotations(adata, modality_name)
            validate_modality(adata, modality_name, error_messages)
    elif epic_type == {"analyses"}:
        for modality_name, adata in mdata.mod.items():
            validate_analyses(adata, modality_name)
            validate_modality(adata, modality_name, error_messages)
    elif isinstance(epic_type, set) and {"annotations", "analyses"}.issubset(epic_type):
        for modality_name, adata in mdata.mod.items():
            validate_analyses(adata, modality_name)
            validate_annotations(adata, modality_name)
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
