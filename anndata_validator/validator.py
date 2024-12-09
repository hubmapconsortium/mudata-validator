import anndata
import os
import scipy.sparse
import warnings

def validate_anndata(input_data):
    """
    Validates an AnnData object or an H5AD file.
    
    Parameters:
    - input_data: str or anndata.AnnData
      Either a path to an H5AD file or an AnnData object.
    
    Raises:
    - ValueError: If the validation fails.
    - Warnings for other things
    
    Returns:
    - None: Prints success if validation passes.
    """
    if isinstance(input_data, str):
        if not os.path.isfile(input_data):
            raise ValueError(f"Provided path '{input_data}' is not a valid file.")
        adata = anndata.read_h5ad(input_data)
    elif isinstance(input_data, anndata.AnnData):
        adata = input_data
    else:
        raise TypeError("Input must be a path to an H5AD file or an AnnData object.")

    # Validation logic
    # REQUIRED: 
    # It will be assumed that the index contains the barcodes/unique identifiers Check for dupliocates.
    print("The values in AnnData.obs.index will be used as the cells' barcodes. They look like:")
    print(adata.obs.head().index)
    if adata.obs.index.duplicated().any():
        warnings.warn("The `.obs` index contains duplicate values, which may lead to issues. Is this data from multiple datasets?", UserWarning)
        warnings.warn("If this data is from multiple datasets, consider prepending the barcode with the plate/well, HuBMAP ID, or the HuBMAP UUID.")
    
    # There must be a column with the original barcodes even if it is the same as the index. 
    if "original_obs_id" not in adata.obs.columns:
        raise ValueError("`.obs` must contain a column named 'original_obs_id' containing the original barcode or unique identifier even if no appending/transformation necessary.")
    
    # There must be a column in .obs with the observation type ontology ID called 'object_type'
    if "object_type" not in adata.obs.columns:
        raise ValueError("`.obs` must contain a column named 'object_type' containing the observation type ontology ID (cell/nucleus).")
    
    # !!TODO!! ID for measured variable (what's being measured: gene, m/z ratio)

    # Main data matrix .X (observation-by-variable)
    # !!TODO!! Do I need to check anything for this besides if it's sparse?
    # Check if .X is sparse and put in CSR if sparse.
    if scipy.sparse.issparse(adata.X):
        if not isinstance(adata.X, scipy.sparse.csr_matrix):
            adata.X = adata.X.tocsr()
    
    print("If this is spatial data, coordinates should go in .obsm['X_spatial']")
    # !!TODO!! Value somewhere in the object for what type of embedding? (UMAP, tSNE, Harmony, etc.)
    # !!TODO!! Should I throw a warning if no X_embedding, an error, or nothing?
    if 'X_embedding' not in adata.obsm:
        warnings.warn("The `.obsm` does not contain an entry called 'X_embedding'. Any coordinates for display in Vitessce should go here.", UserWarning)
    # Additional embeddings can be added as arbitrary other keys in .obsm !!TODO!! What to do about this?

    # Check for Protocol DOI in `.uns['protocol']`
    if 'protocol' not in adata.uns or not adata.uns['protocol']:
        raise ValueError("`.uns` must contain a 'protocol' key with a valid Protocol DOI.")
    
    # !!TODO!! If there are parts of the object that look like they should be named something else for us to read them, give a warning ? 

    # Recommended:
    # General annotation storage
    if 'annotation' not in adata.obsm:
        warnings.warn("It is recommmended to use `.obsm['annotation]` for general annotation storage.")
    elif 'annotation' in adata.obsm:
        if 'annotation' not in adata.uns:
            raise ValueError("`.obsm['annotation']` exists, but `.uns['annotation_methods']` is missing.")
    # Annotation details (prediction score, multiple levels, etc.) from each method can be stored in a separate .obsm key (e.g., .obsm['azimuth'] or .obsm['cluster'])



    # Matt said he'll give me a more detailed bulleted list than what was in the document, follow up
    
    print("Validation passed!")