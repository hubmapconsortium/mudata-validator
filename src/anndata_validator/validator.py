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

    # Track accessed columns and keys dynamically
    accessed_obs_columns = set()
    accessed_obsm_keys = set()
    accessed_uns_keys = set()

    # REQUIRED: Check for duplicate values in the index
    print("The values in AnnData.obs.index will be used as the cells' barcodes. They look like:")
    print(adata.obs.head().index)
    check_duplicate_objects(adata.obs, error_messages)
    
    # Validate `.obs` fields
    if "original_obs_id" in adata.obs.columns:
        accessed_obs_columns.add("original_obs_id")
    else:
        error_messages.append("`.obs` must contain a column named 'original_obs_id' containing the original barcode or unique identifier.")

    if "object_type" in adata.obs.columns:
        accessed_obs_columns.add("object_type")
    else:
        error_messages.append("`.obs` must contain a column named 'object_type' containing the observation type ontology ID (cell/nucleus).")

    # Validate `.uns` for protocol DOI
    if "protocol" in adata.uns and adata.uns["protocol"]:
        accessed_uns_keys.add("protocol")
    else:
        error_messages.append("`.uns` must contain a key 'protocol' with a valid Protocol DOI.")
    
    # Recommended: Annotation storage in `.obsm['annotation']`
    if "annotation" in adata.obsm:
        accessed_obsm_keys.add("annotation")
        if "annotation_methods" not in adata.uns:
            error_messages.append("`.obsm['annotation']` exists, but `.uns['annotation_methods']` is missing.")
        accessed_uns_keys.add("annotation_methods")
    else:
        warnings.warn("It is recommended to use `.obsm['annotation']` for general annotation storage.", UserWarning)

    # Check sparsity for all matrices
    check_sparsity(adata.X, ".X")

    for layer, key_set in [
        (adata.layers, set()),
        (adata.obsm, accessed_obsm_keys),
        (adata.obsp, set()),
        (adata.varm, set()),
        (adata.varp, set()),
    ]:
        if hasattr(layer, "keys"):
            for key in layer.keys():
                key_set.add(key)
                check_sparsity(layer[key], f"{layer}[{key}]")
    
    print("If this is spatial data, coordinates should go in .obsm['X_spatial']")
    if 'X_spatial' in adata.obsm:
        accessed_obsm_keys.add('X_spatial')

    # Check for embedding coordinates
    if 'X_embedding' not in adata.obsm:
        warnings.warn("The `.obsm` does not contain an entry called 'X_embedding'. Any coordinates for display in Vitessce should go here.", UserWarning)
        missing_keys = []
        if 'X_umap' in adata.obsm:
            missing_keys.append('X_umap')
        if 'X_tsne' in adata.obsm:
            missing_keys.append('X_tsne')
        if 'X_harmony' in adata.obsm:
            missing_keys.append('X_harmony')
        
        if missing_keys:
            warnings.warn(f"Found the following embeddings but not 'X_embedding': {', '.join(missing_keys)}. Consider copying these matrices to .obsm['X_embedding'].", UserWarning)
    else:
        accessed_obsm_keys.add('X_embedding')

    # Print all unused `.obs` columns and `.obsm` keys
    unused_obs_columns = [col for col in adata.obs.columns if col not in accessed_obs_columns]
    unused_obsm_keys = [key for key in adata.obsm.keys() if key not in accessed_obsm_keys]
    unused_uns_keys = [key for key in adata.uns.keys() if key not in accessed_uns_keys]

    if unused_obs_columns:
        warnings.warn(f"Unused .obs columns: {', '.join(unused_obs_columns)}", UserWarning)
    if unused_obsm_keys:
        warnings.warn(f"Unused .obsm keys: {', '.join(unused_obsm_keys)}", UserWarning)
    if unused_uns_keys:
        warnings.warn(f"Unused .uns keys: {', '.join(unused_uns_keys)}", UserWarning)

    # Raise an error if validation fails
    if error_messages:
        formatted_errors = "\n- ".join(error_messages)
        raise ValueError(f"Validation failed with the following issues:\n- {formatted_errors}")
    
    print("Validation passed!")
