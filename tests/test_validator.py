import pytest
import anndata
import numpy as np
import pandas as pd
import scipy.sparse
from anndata_validator import validate_anndata  

def create_test_anndata():
    """
    Creates a valid AnnData object for testing.
    """
    obs_data = pd.DataFrame({
        "original_obs_id": ["cell1", "cell2", "cell3"],
        "object_type": ["cell", "nucleus", "cell"]
    }, index=["cell1", "cell2", "cell3"])
    
    var_data = pd.DataFrame(index=["gene1", "gene2", "gene3"])
    X_data = scipy.sparse.csr_matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    obsm_data = {"X_embedding": np.array([[1, 2], [3, 4], [5, 6]])}
    uns_data = {"protocol": "doi:10.1234/protocol"}

    return anndata.AnnData(X=X_data, obs=obs_data, var=var_data, obsm=obsm_data, uns=uns_data)

def test_valid_anndata():
    """
    Test that a valid AnnData object passes validation.
    """
    adata = create_test_anndata()
    validate_anndata(adata)  # Should not raise any exceptions

def test_missing_original_obs_id():
    """
    Test that missing 'original_obs_id' raises a ValueError.
    """
    adata = create_test_anndata()
    adata.obs = adata.obs.drop(columns=["original_obs_id"])
    
    with pytest.raises(ValueError, match=r"must contain a column named 'original_obs_id'"):
        validate_anndata(adata)

def test_missing_object_type():
    """
    Test that missing 'object_type' raises a ValueError.
    """
    adata = create_test_anndata()
    adata.obs = adata.obs.drop(columns=["object_type"])
    
    with pytest.raises(ValueError, match=r"must contain a column named 'object_type'"):
        validate_anndata(adata)

def test_duplicate_obs_index():
    """
    Test that duplicate indices in .obs produce a warning.
    """
    adata = create_test_anndata()
    adata.obs.index = ["cell1", "cell1", "cell3"]  # Duplicate index
    
    with pytest.warns(UserWarning, match=r"index contains duplicate values"):
        validate_anndata(adata)

def test_missing_X_embedding():
    """
    Test that missing 'X_embedding' in .obsm produces a warning.
    """
    adata = create_test_anndata()
    del adata.obsm["X_embedding"]
    
    with pytest.warns(UserWarning, match=r".obsm does not contain an entry called 'X_embedding'"):
        validate_anndata(adata)

def test_sparse_matrix_conversion():
    """
    Test that sparse .X data is converted to CSR format.
    """
    adata = create_test_anndata()
    adata.X = scipy.sparse.csc_matrix(adata.X)  # Use a non-CSR format initially
    
    validate_anndata(adata)
    assert isinstance(adata.X, scipy.sparse.csr_matrix), "X should be converted to CSR format"

def test_missing_protocol():
    """
    Test that missing 'protocol' in .uns raises a ValueError.
    """
    adata = create_test_anndata()
    del adata.uns["protocol"]
    
    with pytest.raises(ValueError, match=r".uns.*must contain a 'protocol' key"):
        validate_anndata(adata)
