import pytest
from anndata_validator import validate_anndata
import anndata
import pandas as pd

def test_validate_anndata_with_object():
    # Create anndata object that will not raise exception
    pass

def test_validate_anndata_with_path():
    # Create path anndata object that will not raise exception
    pass
    
def test_missing_doi():
    # Create anndata object with missing doi from uns
    pass

def test_misnamed_obms():
    # Create anndata object without proper obsm name to be read for hubmap
    pass
