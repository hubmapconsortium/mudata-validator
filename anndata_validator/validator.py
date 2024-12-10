import anndata
import os
import pandas as pd
import scipy.sparse
import numpy as np
import warnings


def check_duplicate_objects(data: pd.DataFrame, error_messages: list):
    """Check for duplicate object IDs in the data."""
    if len(set(data.index)) == data.shape[0]:
        return
    counts = data.index.value_counts()
    duplicates = counts[counts > 1]
    message_pieces = [
        "Found duplicate object IDs:",
        *(f"\t{i}\t({count} occurrences)" for i, count in duplicates.items()),
    ]
    error_messages.append("\n".join(message_pieces))
    warnings.warn("If this data is from multiple datasets, you must prepend the barcode with the plate/well, HuBMAP ID, or the HuBMAP UUID.")


def check_sparsity(matrix, matrix_name: str):
    """Check the sparsity of a matrix and warn if it's too dense."""
    if isinstance(matrix, np.ndarray):
        sparsity = (scipy.sparse.csr_matrix(matrix).nnz / np.prod(matrix.shape))
        if sparsity > 0.3:
            warnings.warn(f"{matrix_name} is a dense matrix with sparsity {sparsity:.4f}. It is recommended to store this as a sparse matrix.", UserWarning)


def validate_anndata(input_data):
    """
    Validates an AnnData object or an H5AD file.

    Parameters:
    - input_data: str or anndata.AnnData
      Either a path to an H5AD file or an AnnData object.

    Raises:
    - ValueError: If validation fails with error messages.
    - Warnings for non-critical issues.

    Returns:
    - None: Prints success if validation passes.
    """
    error_messages = []

    if isinstance(input_data, anndata.AnnData):
        adata = input_data
    else:
        adata = anndata.read_h5ad(input_data)

    # Track accessed columns and keys
    accessed_obs_columns = set()
    accessed_obsm_keys = set()

    # REQUIRED: Check for duplicate values in the index
    print("The values in AnnData.obs.index will be used as the cells' barcodes. They look like:")
    print(adata.obs.head().index)
    check_duplicate_objects(adata.obs, error_messages)
    
    # Accessed columns
    accessed_obs_columns.add("original_obs_id")
    accessed_obs_columns.add("object_type")

    # There must be a column with the original barcodes, even if it is the same as the index
    if "original_obs_id" not in adata.obs.columns:
        error_messages.append("`.obs` must contain a column named 'original_obs_id' containing the original barcode or unique identifier, even if no appending/transformation necessary.")
    
    # There must be a column in .obs with the observation type ontology ID called 'object_type'
    if "object_type" not in adata.obs.columns:
        error_messages.append("`.obs` must contain a column named 'object_type' containing the observation type ontology ID (cell/nucleus).")

        # Check for Protocol DOI in `.uns['protocol']`
    if 'protocol' not in adata.uns or not adata.uns['protocol']:
        error_messages.append("`.uns` must contain a 'protocol' key with a valid Protocol DOI.")
    
    # Recommended: Annotation storage in .obsm['annotation']
    if 'annotation' not in adata.obsm:
        warnings.warn("It is recommended to use `.obsm['annotation']` for general annotation storage.", UserWarning)
    elif 'annotation' in adata.obsm:
        if 'annotation_methods' not in adata.uns:
            error_messages.append("`.obsm['annotation']` exists, but `.uns['annotation_methods']` is missing.")
    
    # Check sparsity for all matrices
    check_sparsity(adata.X, ".X")

    if hasattr(adata, 'layers') and adata.layers:
        for key in adata.layers:
            check_sparsity(adata.layers[key], f".layers['{key}']")
    
    if hasattr(adata, 'obsm') and adata.obsm:
        for key in adata.obsm:
            if key in ['X_spatial', 'X_embedding']:  # Add keys accessed during validation
                accessed_obsm_keys.add(key)
            check_sparsity(adata.obsm[key], f".obsm['{key}']")
    
    if hasattr(adata, 'obsp') and adata.obsp:
        for key in adata.obsp:
            check_sparsity(adata.obsp[key], f".obsp['{key}']")
    
    if hasattr(adata, 'varm') and adata.varm:
        for key in adata.varm:
            check_sparsity(adata.varm[key], f".varm['{key}']")
    
    if hasattr(adata, 'varp') and adata.varp:
        for key in adata.varp:
            check_sparsity(adata.varp[key], f".varp['{key}']")

    # Print all unused .obs columns and .obsm keys
    unused_obs_columns = [col for col in adata.obs.columns if col not in accessed_obs_columns]
    unused_obsm_keys = [key for key in adata.obsm.keys() if key not in accessed_obsm_keys]
    if unused_obs_columns:
        print(f"Unused .obs columns: {', '.join(unused_obs_columns)}")
    if unused_obsm_keys:
        print(f"Unused .obsm keys: {', '.join(unused_obsm_keys)}")

    # If there are any error messages, raise an exception with all of them
    if error_messages:
        formatted_errors = "\n- ".join(error_messages)
        raise ValueError(f"Validation failed with the following issues:\n- {formatted_errors}")
    
    print("Validation passed!")
