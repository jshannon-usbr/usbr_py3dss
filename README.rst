# USBR DSS Query Tool

Summary
-------
This module reads regular time series data & meta-data from a [HEC-DSS](https://www.hec.usace.army.mil/software/hec-dss/)
file into Python objects.

Acronyms
--------
| Abbreviation | Name                                                
|--------------|-------------------------------------------------------------|
|DSS           | Data Storage System                                         |
|DWR           | State of California Department of Water Resources           |
|HEC           | U.S. Army Corps of Engineers' Hydrologic Engineering Center |
|USBR          | U.S. Department of the Interior, Bureau of Reclamation      |

Dependencies
------------
This module relies on the functions compiled in `bin/heclib_x64.dll`, which was
obtained from HEC.

Acknowledgments
---------------
James Gilbert is the original author of this module and developed a majority
of the code while with USBR. His work is version controlled under the branch
`JMG`.

Contact Information
-------------------
This tool is maintained by USBR in Sacramento, CA. Contact Jim Shannon at
jshannon@usbr.gov.

Future Development
------------------
Future development tasks will focus on the following priorities:
    1. Remove dependency on `numpy`;
    2. Expand documentation;
    3. Address documentation wishlists in `dss3_functions_reference.py`.
    
Alternatives
------------
`pyhecdess` is a well documented and stable library to query data from DSS into
a `pandas` DataFrame. DWR maintains the library [here](https://github.com/CADWRDeltaModeling/pyhecdss)
with online documentation hosted [here](https://cadwrdeltamodeling.github.io/pyhecdss/html/index.html).
