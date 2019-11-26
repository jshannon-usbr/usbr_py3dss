"""
Summary
-------
This Python module initialization file imports functions from
"dss3_functions_reference.py" for use in the `usbr_py3dss` library.

"""
# %% Import libraries.
# Import custom libraries.
from .dss3_functions_reference import *


# %% Execute script.
if __name__ == '__main__':
    msg = ('This module is intended to be imported for use into another '
           'module. It is not intended to be run as a __main__ file.')
    print(msg)
