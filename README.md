# MuData Validator

A Python package to validate mudata objects based on specific criteria for HuBMAP.

## Installation

Install the package:
```bash
pip install hubmap-mudata-validator
```

Use the package:
```python
from mudata_validator import validate_mudata
validate_mudata(mudata_object)
# OR
validate_mudata(path_to_h5mu)
```

Example data:
An example H5MU file is available at data/example.h5mu.

Local validation:
To run the HuBMAP validator locally, open and follow the steps in validate_example.ipynb, which walks through validating the example dataset.