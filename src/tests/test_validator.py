import anndata
import muon as mu
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from mudata_validator import validate_mudata


@pytest.fixture
def create_mudata():
    """Fixture to create a mock MuData object for testing."""
    # Modality 1
    obs1 = pd.DataFrame(
        {"original_obs_id": ["A", "B", "C"], "object_type": ["cell", "cell", "cell"]},
        index=["A", "B", "C"],
    )
    var1 = pd.DataFrame(index=["gene1", "gene2", "gene3"])
    X1 = sp.random(3, 3, density=0.1, format="csr")
    adata1 = anndata.AnnData(X=X1, obs=obs1, var=var1)
    adata1.uns["protocol"] = "DOI:whatever/protocol"
    adata1.uns["analyte_class"] = "RNA"

    # Modality 2
    obs2 = pd.DataFrame(
        {"original_obs_id": ["X", "Y", "Z"], "object_type": ["cell", "cell", "cell"]},
        index=["X", "Y", "Z"],
    )
    var2 = pd.DataFrame(index=["geneA", "geneB", "geneC"])
    X2 = sp.random(3, 3, density=0.1, format="csr")
    adata2 = anndata.AnnData(X=X2, obs=obs2, var=var2)
    adata2.uns["protocol"] = "DOI:whatever/protocol"
    adata2.uns["analyte_class"] = "RNA"

    # Combine into MuData
    mdata = mu.MuData({"modality1": adata1, "modality2": adata2})
    mdata.uns["class_types"] = "cell"
    mdata.uns["epic_type"] = {"analyses", "annotations"}
    return mdata


def test_validate_mudata_valid(create_mudata):
    """Test that a valid MuData object passes validation."""
    mdata = create_mudata
    try:
        validate_mudata(mdata)
    except ValueError as e:
        pytest.fail(f"Unexpected ValueError: {e}")


def test_missing_protocol_raises_error(create_mudata):
    """Missing protocol in .uns should raise ValueError."""
    mdata = create_mudata
    del mdata.mod["modality1"].uns["protocol"]

    with pytest.raises(
        ValueError,
        match=r"`modality1.uns` must contain a key 'protocol' with a valid Protocol DOI",
    ):
        validate_mudata(mdata)


def test_missing_original_obs_id_raises_error(create_mudata):
    """Missing original_obs_id column in .obs should raise ValueError."""
    mdata = create_mudata
    mdata.mod["modality2"].obs = mdata.mod["modality2"].obs.drop(
        columns=["original_obs_id"]
    )

    with pytest.raises(
        ValueError,
        match=r"`modality2.obs` must contain a column named 'original_obs_id' containing the original barcode or unique identifier",
    ):
        validate_mudata(mdata)


def test_duplicate_indices_raise_error(create_mudata):
    """Duplicate index values in .obs should raise ValueError."""
    mdata = create_mudata
    mdata.mod["modality2"].obs.index = ["X", "X", "Z"]

    with pytest.raises(ValueError, match=r"Found duplicate object IDs in modality2"):
        validate_mudata(mdata)


def test_dense_matrix_warns(create_mudata):
    """Dense X should trigger a warning."""
    mdata = create_mudata
    mdata.mod["modality1"].X = np.ones((3, 3))

    with pytest.warns(
        UserWarning,
        match=r"modality1.X is a dense matrix with sparsity 1.0000. It is recommended to store this as a sparse matrix.",
    ):
        validate_mudata(mdata)


def test_annotation_obsm_pass(create_mudata):
    """Presence of annotation matrix should not raise an error."""
    mdata = create_mudata
    mdata.mod["modality1"].obsm["annotation"] = np.array(
        ['Immune', 'Immune', 'Immune']
    )
    mdata.mod["modality1"].uns['annotation_methods'] = 'manual'

    try:
        validate_mudata(mdata)
    except ValueError as e:
        pytest.fail(f"Unexpected ValueError: {e}")


def test_spatial_coords_pass(create_mudata):
    """Presence of X_spatial should not raise an error."""
    mdata = create_mudata
    mdata.mod["modality1"].obsm["X_spatial"] = np.array(
        [[0.1, 0.2], [0.3, 0.4], [0.5, 0.6]]
    )

    try:
        validate_mudata(mdata)
    except ValueError as e:
        pytest.fail(f"Unexpected ValueError: {e}")


def test_spatial_coords_raises_error(create_mudata):
    """Only one dimension in X_Spatial should raise ValueError."""
    mdata = create_mudata
    mdata.mod["modality1"].obsm["X_spatial"] = np.array(
        [0.1, 0.2, 0.3]
    )
    
    with pytest.raises(
        ValueError,
        match=r"Only one dimension found in modality1.obsm\['X_spatial'\]; Expecting at least \(y,x\).",
    ):
        validate_mudata(mdata)


def test_invalid_object_type_raises_error(create_mudata):
    """Invalid value in object_type column should raise ValueError."""
    mdata = create_mudata
    mdata.mod["modality1"].obs["object_type"] = "invalid_value"

    with pytest.raises(
        ValueError,
        match=r"'modality1.obs\['object_type'\]' contains invalid values: invalid_value. Allowed values are: cell, nuclei, ftu, spot.",
    ):
        validate_mudata(mdata)


def test_missing_epic_type_raises_error(create_mudata):
    """Missing 'epic_type' in mdata.uns should raise ValueError."""
    mdata = create_mudata
    del mdata.uns["epic_type"]

    with pytest.raises(
        ValueError,
        match=r"MuData.uns must contain a key called 'epic_type' with at least one valid epic type: annotations, analyses",
    ):
        validate_mudata(mdata)


def test_invalid_analyte_class_raises_error(create_mudata):
    """Invalid analyte_class should raise ValueError."""
    mdata = create_mudata
    mdata.mod["modality1"].uns["analyte_class"] = "invalid"

    with pytest.raises(
        ValueError,
        match=r"The value in `modality1.uns\['analyte_class'\]` must reference a known analyte class.",
    ):
        validate_mudata(mdata)
