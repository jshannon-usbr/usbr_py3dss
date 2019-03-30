# -*- coding: utf-8 -*-
r"""
Summary
-------
Python wrapper for DSS HEC-LIB functions
    - Originally developed summer 2017
    - Functions added as needed, although not all are guaranteed to work as
      intended/needed

Updated 2018.04.05-2018.04.10 to include write capabilities, additional catalog
                              functions

TODO:  Error catching with helpful return info (based on HECLIB documentation)
TODO:  Sort out catalog creation and write functionality
TODO:  Add/retrieve more metadata with each record on read/write
TODO:  Sort out single/double precision consistency on read/write operations

@author: jmgilbert

JAS Wishlist
------------
1. Suppress messages when querying data unless critical errors occur.
2. Delete *.dsc and *.dsk when finished querying.
3. Cleaner code and documentation.

"""
# %% Import libraries.
# Import standard libraries.
import ctypes as C
import os
# Import third party libraries.
import numpy as np
# %% Initialize variables.
# TODO: Explain the purpose of these variables.
# <JAS 2019-03-28>
global CATALOGUNIT
global CONDUNIT
CATALOGUNIT = 12
CONDUNIT = 13
# Set path of heclib_x64.dll.
module_dir = os.path.dirname(__file__)
dll_path = os.path.join(module_dir, r'bin\heclib_x64.dll')
dsslib = C.cdll.LoadLibrary(dll_path)


# %% Define functions.
def open_dss(fpi):
    r"""
    Summary
    -------
    No documentation as of 2019-03-29.

    """
    zopen = getattr(dsslib, 'ZOPEN_')
    zopen.argparse = [C.c_long*600, C.c_char*800, C.c_long, C.c_long]
    zopen.restype = None
    ifltab = (C.c_long * 600)()
    iostat = C.c_long()
    #fp = (C.c_char*800)(*fpi)
    #print(fp)
    zopen(C.byref(ifltab), fpi.encode('ascii'), C.byref(iostat), len(fpi.encode('ascii')))
    return([ifltab, iostat.value])


def close_dss(ifltab):
    r"""
    Summary
    -------
    No documentation as of 2019-03-29.

    """
    zclose = getattr(dsslib, 'ZCLOSE_')
    zclose.argparse = [C.c_long*600]
    stdout = zclose(ifltab)
    return(stdout)


def dss_exists(fp):
    r"""
    Summary
    -------
    No documentation as of 2019-03-29.

    """
    zfname = getattr(dsslib, 'ZFNAME_')
    zfname.argparse = [C.c_char*800, C.c_char*800, C.c_long, C.c_bool,
                       C.c_long, C.c_long]
    # arguments: path in, path with extension if exists[out], length of returned name,
    #            file exists logical, length of path in, length of path out
    zfname.restype = None
    plen = C.c_long()
    lexist = C.c_bool()
    new_fp = (C.c_char*800)()

    zfname(fp.encode('ascii'), new_fp, C.byref(plen),
           C.byref(lexist), len(fp.strip()), len(new_fp))

    #print(lexist.value)
    #print(new_fp)
    #print(plen)
    return([lexist.value, new_fp[0:].decode().strip()])


def read_regts(ifltab, cpath, cdate, ctime, nvalsi):
    r"""
    Summary
    -------
    No documentation as of 2019-03-29.

    """
    zrrts = getattr(dsslib, 'ZRRTS_')
# Args: (IFLTAB(600)  | CPATH*80   | CDATE   | CTIME   | NVALS   | VALUES   |  CUNITS   | CTYPE   | IOFSET   | ISTAT)
# DTYPE:(INTEGER(600) | CHAR*80    | CHAR*20 | CHAR*4  | INTEGER | REAL(NVALS)| CHAR*8  | CHAR*8  | INTEGER  | INTEGER)
#        INPUT        | INPUT      | INPUT   | INPUT   | INPUT   | OUTPUT   | OUTPUT    | OUTPUT | OUTPUT    | OUTPUT)
    zrrts.argparse = [C.c_long*600, C.c_char*800, C.c_char*20, C.c_char*4,
                      C.c_long, C.c_float, C.c_char*8, C.c_char*8, C.c_long,
                      C.c_long, C.c_long, C.c_long, C.c_long, C.c_long,
                      C.c_long]  # lengths of: cpath, cdate, ctime, cunits, ctype
    zrrts.restype = None

    nvals = C.c_long(nvalsi)
    vals = (C.c_float*nvalsi)()
    cunits = (C.c_char*8)()
    ctype = (C.c_char*8)()
    iofset = C.c_long()
    istat = C.c_long()
    zrrts(ifltab, cpath.encode('ascii'), cdate.encode('ascii'),
          ctime.encode('ascii'), C.byref(nvals), C.byref(vals), cunits, ctype,
          C.byref(iofset), C.byref(istat), len(cpath), len(cdate), len(ctime),
          len(cunits), len(ctype))

    return([nvals.value, np.float32(vals), cunits[0:].decode('utf-8').strip(),
            ctype[0:].decode('utf-8').strip(), iofset.value, istat.value])


def read_regtsd(ifltab, cpath, cdate, ctime, nvalsi, lgetdob_in=True):
    r"""
    Summary
    -------
    No summary as of 2019-03-29.

    Parameters
    ----------
    No documentation of Parameters as of 2019-03-29.

    Returns
    -------
    - 0 : nvals.value
    - 1 : dVals
    - 2 : cunits[0:].decode('utf-8').strip()
    - 3 : ctype[0:].decode('utf-8').strip()
    - 4 : iofset.value
    - 5 : istat.value
    - 6 : csupp[0:].decode('utf-8')
    - 7 : coords_info
    - 8 : ctzone[0:].decode('utf-8')
    - 9 : jqual
    - 10 : lfildob
    - 11 : itzone

    Notes
    -----
    # Args: (IFLTAB(600)  | CPATH*80   | CDATE   | CTIME   | KVALS   | NVALS   |
    #        INPUT        | INPUT      | INPUT   | INPUT   | INPUT   | OUTPUT   | OUTPUT    | OUTPUT | OUTPUT    | OUTPUT)
    # DTYPE:(INTEGER(600) | CHAR*80    | CHAR*20 | CHAR*4  | INTEGER | INTEGER |

    # Args: LGETDOB | LFILDOB | sVALUES    | dVALUES       | JQUAL | LQUAL | LQREAD |
    #       LOGICAL | LOGICAL | REAL(NVALS)| Double(NVALS) |INTEGER|LOGICAL| LOGICAL|

    #Args: CUNITS   | CTYPE   | CSUPP   | IOFSET   | JCOMP | ITZONE |CTZONE| COORDS |
    #       CHAR*8  | CHAR*8  | CHAR*80 | INTEGER  |INTEGER|INTEGER |CHAR*30 | DOUBLE |

    #Args: ICDESC |LCOORDS | ISTAT | L_CPATH | L_CDATE | L_CTIME | L_CUNITS | L_CTYPE | L_CSUPP | L_CTZONE)
    #      INTEGER|INTEGER |INTEGER| INTEGER | INTEGER | INTEGER | INTEGER  | INTEGER | INTEGER | INTEGER

    """
    zrrtsc = getattr(dsslib, 'ZRRTSC_')
    zrrtsc.argparse = [
            C.c_long*600, C.c_char*800, C.c_char*20, C.c_char*4, C.c_long, C.c_long,  # IFLTAB, CPATH, CDATE, CTIME, KVALS, NVALS
            C.c_bool, C.c_bool, C.c_float, C.c_double, C.c_long, C.c_bool, C.c_bool,  # LGETDOB, LFILDOB, sVALUES, dVALUES, JQUAL, LQUAL, LQREAD
            C.c_char*8, C.c_char*8, C.c_char*80, C.c_long, C.c_long, C.c_long, C.c_char*30, C.c_double, # CUNITS, CTYPE, CSUPP, IOFSET, JCOMP, ITZONE, CTZONE, COORDS
            C.c_long, C.c_long, C.c_long, C.c_long, C.c_long, C.c_long, C.c_long, C.c_long, C.c_long, C.c_long ] # ICDESC, LCOORDS< ISTAT, L_CPATH, L_CDATE, L_CTIME, L_CUNITS, L_CTYPE, L_CSUPP L_CTZONE

    zrrtsc.restype = None

    nvals = kvals = C.c_long(nvalsi)
    sVals = (C.c_float*nvalsi)()
    dVals = (C.c_double*nvalsi)()
    lgetdob = C.c_bool(lgetdob_in)
    lfildob = C.c_bool()
    jqual = (C.c_long*nvalsi)()
    lqual = C.c_bool()
    lqread = C.c_bool()
    csupp = (C.c_char*80)()
    cunits = (C.c_char*8)()
    ctype = (C.c_char*8)()
    iofset = C.c_long()
    istat = C.c_long()
    jcomp = C.c_long()
    itzone = C.c_long()
    ctzone = (C.c_char*30)()
    coords = (C.c_double*3)()
    icdesc = (C.c_long*6)()
    lcoords = C.c_bool()

    ecpath = cpath.encode('ascii')
    ecdate = cdate.encode('ascii')
    ectime = ctime.encode('ascii')

    l_cpath = len(ecpath)
    l_cdate = len(ecdate)
    l_ctime = len(ectime)
    l_cunits = len(cunits)
    l_ctype = len(ctype)
    l_csupp = len(csupp)
    l_ctzone = len(ctzone)

    zrrtsc(ifltab, ecpath, ecdate, ectime, C.byref(kvals), C.byref(nvals), C.byref(lgetdob), C.byref(lfildob),   #8
           C.byref(sVals),  C.byref(dVals), C.byref(jqual), C.byref(lqual), C.byref(lqread),    #+5 = 13
           cunits, ctype, csupp, C.byref(iofset), C.byref(jcomp), C.byref(itzone),            #+6 = 19
           ctzone, C.byref(coords), C.byref(icdesc), C.byref(lcoords), C.byref(istat),        #+5 = 24
           l_cpath, l_cdate, l_ctime, l_cunits, l_ctype, l_csupp, l_ctzone)                   #7 = 31

    # icdesc is a 6-item list with information on the coordinates
    # icdesc[0] = Coordiante system: (following values follow the DSSVue convention - as far as I can tell though, you can use whatever scheme you want)
    #               0 = No coordinates set
    #               1 = Lat/Long
    #               2 = State Plane/FIPS
    #               3 = State Plane/ADS
    #               4 = UTM
    #               5 = Local/other
    # icdesc[1] = coordinate ID (some integer)
    # icdesc[2] = Datum Units
    #             0 = not specified
    #             1 = English (ft or miles i guess?)
    #             2 = SI   (m or km?)
    #             3 = Decimal degrees
    #             4 = degrees, minutes, seconds
    # icdesc[3] = horizontal datum
    #          0 = unset
    #          1 = NAD83
    #          2 = NAD27
    #          3 = WAS84
    #          4 = WAS72
    #          5 = local/other
    #

    # compile coordinate information into a single dictionary for clarity
    coords_info = icdesc_to_dict(coords, icdesc)

#            ByVal strCDate As String, ByVal strCTime As String, _
#            lngKVals As Long, lngNVals As Long, lngLGETDOB As Long, _
#            lngLFILDOB As Long, sngValues As Single, dblValues As Double, _
#            lngJQUAL As Long, lngLQUAL As Long, lngLQREAD As Long, _
#            ByVal strCUnits As String, ByVal strCType As String, _
#            ByVal strCSUPP As String, lngIOFSET As Long, _
#            lngJCOMP As Long, lngITZONE As Long, ByVal strCTZONE As String, _
#            dblCOORDS As Double, lngICDESC As Long, lngLCoords As Long, _
#            lngISTAT As Long, ByVal lngL_CPath As Long, ByVal lngL_CDate As Long, _
#            ByVal lngL_CTime As Long, ByVal lngL_CUnits As Long, _
#            ByVal lngL_CType As Long, ByVal lngL_CSUPP As Long, ByVal lngL_CTZONE As Long

    return([nvals.value, dVals, cunits[0:].decode('utf-8').strip(), ctype[0:].decode('utf-8').strip(),
            iofset.value, istat.value, csupp[0:].decode('utf-8'), coords_info, ctzone[0:].decode('utf-8'), jqual, lfildob, itzone])


def icdesc_to_dict(coords, icdesc):
    r"""
    Summary
    -------
    Converts the ICDESC list as returned from a DSS file to a dictionary with
    keys as identifiers, according to DSSVue scheme.

    No additional documentation as of 2019-03-29.

    """
    coords_info = {}

    coords_info['X_Long'] = coords[0]
    coords_info['Y_Lat'] = coords[1]

    coordsys = icdesc[0]
    if coordsys == 0:
        coords_info['CoordSys'] = ['No coordinates set', 0]
    elif coordsys == 1:
        coords_info['CoordSys'] = ['LatLong', 1]
    elif coordsys == 2:
        coords_info['CoordSys'] = ['State Plane/FIPS', 2]
    elif coordsys == 3:
        coords_info['CoordSys'] = ['State Plane/ADS', 3]
    elif coordsys == 4:
        coords_info['CoordSys'] = ['UTM', 4]
    elif coordsys == 5:
        coords_info['CoordSys'] = ['Local/other', 5]
    else:
        coords_info['CoordSys'] = ['Hard to say, really', coordsys]

    coordID = icdesc[1]
    coords_info['CoordID'] = coordID

    datumUnits = icdesc[2]
    if datumUnits == 0:
        coords_info['DatumUnit'] = ['Not specified', 0]
    elif datumUnits == 1:
        coords_info['DatumUnit'] = ['English', 1]
    elif datumUnits == 2:
        coords_info['DatumUnit'] = ['SI', 2]
    elif datumUnits == 3:
        coords_info['DatumUnit'] = ['Decimal Degrees', 3]
    elif datumUnits == 4:
        coords_info['DatumUnit'] = ['degrees, minutes, seconds', 4]
    else:
        coords_info['DatumUnit'] = ['Unknown', None]

    datum = icdesc[3]
    #print("Datum from file is: %s" %datum)
    if datum == 0:
        coords_info['Datum'] = ['Not specified', 0]
    elif datum == 1:
        coords_info['Datum'] = ['NAD83', 1]
    elif datum == 2:
        coords_info['Datum'] = ['NAD27', 2]
    elif datum == 3:
        coords_info['Datum'] = ['WAS84', 3]
    elif datum == 4:
        coords_info['Datum'] = ['WAS72', 4]
    elif datum == 5:
        coords_info['Datum'] = ['local/other', 5]
    else:
        coords_info['Datum'] = ['Unknown', None]

    return(coords_info)


def open_catalog(fp, icunitin, lgenca_in=True, lgencd_in=False):
    r"""
    Summary
    -------
    No documentation as of 2019-03-29.

    """
    # need to open catalog before you can read it - this gets called from 'get_catalog'
    zopnca = getattr(dsslib, 'ZOPNCA_')
    # Args:   CDSSFI  | ICUNIT   | LGENCA   |LOPNCA   | LCATLG   | ICDUNT   | LGENCD  |LOPNCD   | LCATCD  |NRECS  | *Len(CDSSFI)
    #         CHAR*64 | INT      | LOGICAL  | LOGICAL  | LOGICAL | INT      | LOGICAL |LOGICAL  |LOGICAL  | INT   | INT
    #          INPUT  | INPUT    |INPUT     | OUT      | OUT     | INPUT    | INPUT  | OUTPUT   | OUTPUT   | OUTPUT  | IN
    zopnca.argparse = [C.c_char,  C.c_long, C.c_bool, C.c_bool, C.c_bool, C.c_long, C.c_bool, C.c_bool, C.c_bool, C.c_long, C.c_long ]
    zopnca.restype = None
    icunit = C.c_long(icunitin)
    lgenca = C.c_bool(lgenca_in)
    lopnca = C.c_bool()
    lcatlg = C.c_bool()
    icdunt = C.c_long()  #@jmg 2017.11.30 - set this explicitly
    lgencd = C.c_bool(lgencd_in)
    lopncd = C.c_bool()
    lcatcd = C.c_bool()
    nrecs = C.c_long()
    #C.byref(icunit)
    zopnca(fp.encode('ascii'), C.byref(icunit), C.byref(lgenca),
           C.byref(lopnca), C.byref(lcatlg), C.byref(icdunt), C.byref(lgencd),
           C.byref(lopncd), C.byref(lcatcd), C.byref(nrecs),
           len(fp.encode('ascii')))
    return([lgenca.value, lopnca.value, lcatlg.value, lgencd.value,
            lopncd.value, lcatcd.value, nrecs.value])


def read_catalog(lopnca, icunitin=12):
    r"""
    Summary
    -------
    No documentation as of 2019-03-29.

    """
    # zrdcat - read pathnames from catalog
    # first - check to see if catalog is open - if not, prompt to open it
    if not lopnca:
        print("Catalog is not open - do this first using open_catalog function")
        return([None, None])

    zrdcat = getattr(dsslib, 'ZRDCAT_')
    # ARGS:  ICUNIT   |  LALL | IOUNIT  |CTAGS        | NDIM  | CPATHS | NPATHS  | NFOUND | *Len(ctags) | *Len(cpath)
    #DTYPES:  INT    | LOGICAL| INT     |CHAR(NDIM)*8 | INT   |
    zrdcat.argparse = [C.c_long, C.c_bool, C.c_long, C.c_char, C.c_long,
                       C.c_char, C.c_long, C.c_long, C.c_long, C.c_long]
    icunit = C.c_long(icunitin)
    lall = C.c_bool(True)
    iounit = C.c_long(0)
    ndim = C.c_long(1) #nrecs  #C.c_long()
    ctags = (C.c_char*8)()
    cpaths = (C.c_char*80)()
    npaths = C.c_long()
    nfound = C.c_long(1)
    path_list = []

    while nfound.value == 1:
        zrdcat(C.byref(icunit), C.byref(lall), C.byref(iounit), ctags,
               C.byref(ndim), cpaths, C.byref(npaths), C.byref(nfound),
               len(ctags), len(cpaths))
        #print(nfound.value)
        #print(cpaths.value)
        #print(ctags.value)
        nfoundv = nfound.value
        if nfoundv == 1:
            path_list.append(cpaths.value.decode('utf-8').strip())
            #print(cpaths.value)
    catstatus  = fortran_close_file(icunitin)  # close the catalog file
    if catstatus == 0:
        lopnca= False
    else:
        lopnca= True
    return([path_list, lopnca])


def fortran_close_file(iunit_in):
    fortclose = getattr(dsslib,'FORTRANCLOSE_')
    # ARGES: IUNIT
    # DTYPES: INT
    fortclose.argparse = [C.c_long]
    iunit = C.c_long(iunit_in)
    status = fortclose(C.byref(iunit))
    return(status)


def write_regts(ifltab, cpath, cdate, ctime, nvalsi, valsi, cunits, ctype,
                _IPLAN=0):
    r"""
    Summary
    ------
    This function writes a regular time series, stored as a Python list, to a
    DSS file. The target DSS file to which the data should be written needs to
    be opened before it can be written to. Opening the file (using `open_dss`
    function) will create a new file if the file doesn't already exist.

    Parameters
    ----------
    ifltab :
        should be a long integer retrieved when opening a file
    cpath :
        string of a properly formatted DSS pathname, should have 7 forward
        slashes /, but no D-part (start date - this comes from the cdate
        argument)
    cdate :
        a string representing the start date - should be of the form
        "DDmonYYYY", e.g. "31Oct1999" for Oct 31, 1999
    ctime :
        a string representing start time in hours since midnight on `cdate`,
        common practice is to use "2400" for monthly data
    nvalsi :
        number of values to be written, as an integer; the end date is defined
        by the start date, number of values, and the period length given in the
        E-part of the variable path
    valsi :
        a python list of length `nvalsi` to be written to the DSS file; missing
        values should be filled with -901.0
    cunits :
        a string indicating units of the time series values
    ctype :
        a string indicating the time period summary type that the values
        represent, commonly "PER-AVG" for period average
    _IPLAN : (default = 0)
        flag to indicate whether to write over existing data
            0 = always write over existing data
            1 = only replace missing data flags in the record, as indicated by
                (-901)
            4 = if an input value is missing (-901), do not allow it to replace
                a non-missing value

    Returns
    -------
    istat :
        flag indicating if data was stored successfully (0), if all data was
        missing (-901), or other "fatal" errors(values >10, see HECLIB
        documentation for details for now...someday I'll include code to catch
        the errors and tell you something useful)

    """
    zsrts = getattr(dsslib, 'ZSRTS_')
# Args: (IFLTAB(600)  | CPATH*80   | CDATE   | CTIME   | NVALS   | VALUES   |  CUNITS   | CTYPE   | IPLAN   | ISTAT)
# DTYPE:(INTEGER(600) | CHAR*80    | CHAR*20 | CHAR*4  | INTEGER | REAL(NVALS)| CHAR*8  | CHAR*8  | INTEGER  | INTEGER)
#        INPUT        | INPUT      | INPUT   | INPUT   | INPUT   | OUTPUT   | OUTPUT    | OUTPUT | OUTPUT    | OUTPUT)
    zsrts.argparse = []

    zsrts.argparse = [C.c_long*600, C.c_char*800, C.c_char*20, C.c_char*4,
                      C.c_long, (C.c_float*nvalsi)(), C.c_char*8, C.c_char*8,
                      C.c_long, C.c_long,
                      C.c_long, C.c_long, C.c_long, C.c_long, C.c_long]  # lengths of: cpath, cdate, ctime, cunits, ctype
    zsrts.restype = None
    valsi = np.float32(valsi) # enforce 32-bit single precision on the values

    nvals = C.c_long(nvalsi)
    vals = (C.c_float*nvalsi)(*valsi)
    #print(vals)  # for debugging
    cunits = cunits  #(C.c_char*8)(cunits)
    ctype = ctype  #(C.c_char*8)(ctype)
    IPLAN = C.c_long(_IPLAN)
    istat = C.c_long()
    zsrts(C.byref(ifltab), cpath.encode('ascii'), cdate.encode('ascii'),
          ctime.encode('ascii'), C.byref(nvals), C.byref(vals),
          cunits.encode('ascii'), ctype.encode('ascii'), C.byref(IPLAN),
          C.byref(istat), len(cpath), len(cdate), len(ctime), len(cunits),
          len(ctype))

    return([istat])


def write_regtsd(ifltab, cpath, cdate, ctime, vals, cunits, ctype,
                 coords=[0.0, 0.0, 0.0], icdesc=[0, 0, 0, 0, 0, 0],
                 csupp='', ctzone='', iplan=0):
    r"""
    Summary
    -------
    Python function to extract regular time series and supplemental information
    from DSS file. This function was written with guidance from HECLIB
    documentation of 'ZSRTSX' [1]_ and source code in DSSExcel.xlam [2]_.

    Parameters
    ----------
    ifltab : c_long_Array_600
        Integer long returned when opening a DSS file with open_dss().
    cpath : str
        Pathname of regular time series to be stored.
    cdate : str
        Starting date of regular time series in "DDMMMYYYY" format. If the year
        is in the 1900s, date can be in "DDMMMYY" format.
    ctime : str
        Time since midnight of starting date.
    vals : list of float
        Sequential list of regular time series.
    cunits : str
        Units of time series data.
    ctype : str
        Type of time series data (e.g. 'PER-AVER', 'PER-CUM', 'INST-VAL')
    coords : list of float, default [0.0, 0.0, 0.0], optional
        Coordinates of time series data.

        - coords[0] : X-coordinate
        - coords[1] : Y-coordinate
        - coords[2] : Unknown; default set to 0.

    icdesc : list of int, default [0, 0, 0, 0, 0, 0], optional
        Metadata for coordinates of time series data.

        +------------+-----------------------------------------------------+
        | List Entry | Description                                         |
        +============+=====================================================+
        | icdesc[0]  | One of the following options for Coordinate System. |
        |            |                                                     |
        |            | - 0 = No Coordinates Set                            |
        |            | - 1 = Latitude/Longitude                            |
        |            | - 2 = State Plane/FIPS                              |
        |            | - 3 = State Plane/ADS                               |
        |            | - 4 = UTM                                           |
        |            | - 5 = Local/other                                   |
        +------------+-----------------------------------------------------+
        | icdesc[1]  | Coordinate ID Number.                               |
        +------------+-----------------------------------------------------+
        | icdesc[2]  | Horizontal Datum Units.                             |
        |            |                                                     |
        |            | - 0 = Not Specified                                 |
        |            | - 1 = English (ft or miles i guess?)                |
        |            | - 2 = SI   (m or km?)                               |
        |            | - 3 = Decimal degrees                               |
        |            | - 4 = degrees, minutes, seconds                     |
        +------------+-----------------------------------------------------+
        | icdesc[3]  | Horizontal Datum.                                   |
        |            |                                                     |
        |            | - 0 = Unset                                         |
        |            | - 1 = NAD83                                         |
        |            | - 2 = NAD27                                         |
        |            | - 3 = WAS84                                         |
        |            | - 4 = WAS72                                         |
        |            | - 5 = local/other                                   |
        +------------+-----------------------------------------------------+
        | icdesc[4]  | Unknown; default set to 0                           |
        +------------+-----------------------------------------------------+
        | icdesc[5]  | Unknown; default set to 0.                          |
        +------------+-----------------------------------------------------+

    csupp : str, default '', optional
        Description of time series data.
    ctzone : str, default '', optional
        Time Zone Identification.
    iplan : int, default 0, optional
        Argument for writing over existing data according to the following
        table from HECLIB documentation of 'ZSRTSX' [1]_.

        +-------+-------------------------------------------------------+
        | iplan | Description                                           |
        +=======+=======================================================+
        | 0     | Always write over existing data.                      |
        +-------+-------------------------------------------------------+
        | 1     | Only replace missing data flags in the record (-901). |
        +-------+-------------------------------------------------------+
        | 4     | If an input value is missing (-901), do not allow it  |
        |       | to replace a non-missing value.                       |
        +-------+-------------------------------------------------------+

    Returns
    -------
    istat : None
        Status indicator of writing operation with the following possible
        returns as stated in the HECLIB documentation of 'ZSRTSX' [1]_.

        +-------+---------------------------------------------------------+
        | istat | Description                                             |
        +=======+=========================================================+
        | 0     | The data was successfully stored.                       |
        +-------+---------------------------------------------------------+
        | 4     | All of the input data provided were missing data flags  |
        |       | (-901).                                                 |
        +-------+---------------------------------------------------------+
        | >10   | A "fatal" error occurred."                              |
        +-------+---------------------------------------------------------+
        | 11    | The number of values to store (nvals) is less than one. |
        +-------+---------------------------------------------------------+
        | 12    | Unrecognized time interval (E part).                    |
        +-------+---------------------------------------------------------+
        | 15    | The starting date or time is invalid.                   |
        +-------+---------------------------------------------------------+
        | 24    | The pathname given does not meet the regular-interval   |
        |       | time series conventions.                                |
        +-------+---------------------------------------------------------+
        | 51    | Unrecognized data compression scheme.                   |
        +-------+---------------------------------------------------------+
        | 53    | Invalid precision exponent specified for the delta      |
        |       | compression method. The precision exponent range is     |
        |       | -6 to +6.                                               |
        +-------+---------------------------------------------------------+

    Notes
    -----
    Variables without the ``_input`` suffix are used directly into the
    heclib_x64.dll ``'ZSRTSC_'`` function. Variables with ``_input`` suffix are
    transformed before input into the heclib_x64.dll ``ZSRTSC_`` function.

    Set missing values to -901.0.

    The following table summarizes arguments for ``'ZSRTSC_'`` based on
    information from HECLIB documentation of 'ZSRTSX' [1_] and source code in
    DSSExcel.xlam [2]_.

    +----------+----------+--------------+--------------+---------------------+
    | Sequence | Argument | Data Type    | Input/Output | Description         |
    +==========+==========+==============+==============+=====================+
    | 1        | IFLTAB   | INTEGER(600) | INPUT        | The DSS work space  |
    |          |          |              |              | used to manage the  |
    |          |          |              |              | DSS file.           |
    +----------+----------+--------------+--------------+---------------------+
    | 2        | CPATH    | CHARACTER*80 | INPUT        | The pathname of the |
    |          |          |              |              | data to store.      |
    +----------+----------+--------------+--------------+---------------------+
    | 3        | CDATE    | CHARACTER*20 | INPUT        | The beginning date  |
    |          |          |              |              | of the time window. |
    +----------+----------+--------------+--------------+---------------------+
    | 4        | CTIME    | CHARACTER*4  | INPUT        | The beginning time  |
    |          |          |              |              | of the time window. |
    +----------+----------+--------------+--------------+---------------------+
    | 5        | NVALS    | INTEGER      | INPUT        | The number of values|
    |          |          |              |              | to store for        |
    |          |          |              |              | SVALUES, defining   |
    |          |          |              |              | the end of the time |
    |          |          |              |              | window.             |
    +----------+----------+--------------+--------------+---------------------+
    | 6        | DOUBLE   | INTEGER      | INPUT        | The number of values|
    |          |          |              |              | to store for        |
    |          |          |              |              | DVALUES, defining   |
    |          |          |              |              | the end of the time |
    |          |          |              |              | window.             |
    +----------+----------+--------------+--------------+---------------------+
    | 7        | SVALUES  | REAL(NVALS)  | INPUT        | List of time series |
    |          |          |              |              | values in sequential|
    |          |          |              |              | order to be stored  |
    |          |          |              |              | in single precision.|
    +----------+----------+--------------+--------------+---------------------+
    | 8        | DVALUES  | REAL(DOUBLE) | INPUT        | List of time series |
    |          |          |              |              | values in sequential|
    |          |          |              |              | order to be stored  |
    |          |          |              |              | in double precision.|
    +----------+----------+--------------+--------------+---------------------+
    | 9        | JQUAL    | INTEGER      | INPUT        | Unknown, but likely |
    |          |          |              |              | an array containing |
    |          |          |              |              | thirty-two bit      |
    |          |          |              |              | flags. Not stored if|
    |          |          |              |              | LQUAL is false.     |
    +----------+----------+--------------+--------------+---------------------+
    | 10       | LQUAL    | LOGICAL      | INPUT        | Variable indicating |
    |          |          |              |              | whether or not to   |
    |          |          |              |              | store JQUAL.        |
    +----------+----------+--------------+--------------+---------------------+
    | 11       | CUNITS   | CHARACTER*8  | INPUT        | The units of the    |
    |          |          |              |              | data (e.g., 'FEET').|
    +----------+----------+--------------+--------------+---------------------+
    | 12       | CTYPE    | CHARACTER*8  | INPUT        | The type of the data|
    |          |          |              |              | (e.g., 'PER-AVER'). |
    +----------+----------+--------------+--------------+---------------------+
    | 13       | COORDS   | REAL         | INPUT        | Coordinates of the  |
    |          |          |              |              | time series data.   |
    +----------+----------+--------------+--------------+---------------------+
    | 14       | NCOORDS  | INTEGER      | INPUT        | Length of COORDS.   |
    +----------+----------+--------------+--------------+---------------------+
    | 15       | ICDESC   | INTEGER      | INPUT        | Metadata for COORDS.|
    +----------+----------+--------------+--------------+---------------------+
    | 16       | NCDESC   | INTEGER      | INPUT        | Length of ICDESC.   |
    +----------+----------+--------------+--------------+---------------------+
    | 17       | CSUPP    | CHARACTER*80 | INPUT        | Description of the  |
    |          |          |              |              | time series data.   |
    +----------+----------+--------------+--------------+---------------------+
    | 18       | ITZONE   | INTEGER      | INPUT        | Time offset in      |
    |          |          |              |              | minutes from UTC.   |
    +----------+----------+--------------+--------------+---------------------+
    | 19       | CTZONE   | CHARACTER*30 | INPUT        | Time Zone ID.       |
    +----------+----------+--------------+--------------+---------------------+
    | 20       | IPLAN    | INTEGER      | INPUT        | Variable indicating |
    |          |          |              |              | whether or not to   |
    |          |          |              |              | write over existing |
    |          |          |              |              | data.               |
    +----------+----------+--------------+--------------+---------------------+
    | 21       | JCOMP    | INTEGER      | INPUT        | Data compression    |
    |          |          |              |              | method.             |
    +----------+----------+--------------+--------------+---------------------+
    | 22       | BASEV    | REAL         | INPUT        | Baseline value for  |
    |          |          |              |              | data compression.   |
    +----------+----------+--------------+--------------+---------------------+
    | 23       | LBASEV   | LOGICAL      | INPUT        | Variable indicating |
    |          |          |              |              | whether or not to   |
    |          |          |              |              | store BASEV.        |
    +----------+----------+--------------+--------------+---------------------+
    | 24       | LDHIGH   | LOGICAL      | INPUT        | Setting for         |
    |          |          |              |              | preallocating       |
    |          |          |              |              | for data            |
    |          |          |              |              | compression.        |
    +----------+----------+--------------+--------------+---------------------+
    | 25       | NPREC    | INTEGER      | INPUT        | Precision exponent  |
    |          |          |              |              | for data            |
    |          |          |              |              | compression.        |
    +----------+----------+--------------+--------------+---------------------+
    | 26       | ISTAT    | INTEGER      | OUTPUT       | Status parameter    |
    |          |          |              |              | indicating success  |
    |          |          |              |              | of storage.         |
    +----------+----------+--------------+--------------+---------------------+
    | 27       | L_Cpath  | INTEGER      | INPUT        | Length of CPATH.    |
    +----------+----------+--------------+--------------+---------------------+
    | 28       | L_CDate  | INTEGER      | INPUT        | Length of CDATE.    |
    +----------+----------+--------------+--------------+---------------------+
    | 29       | L_CTime  | INTEGER      | INPUT        | Length of CTIME.    |
    +----------+----------+--------------+--------------+---------------------+
    | 30       | L_CUnits | INTEGER      | INPUT        | Length of CUNITS.   |
    +----------+----------+--------------+--------------+---------------------+
    | 31       | L_CType  | INTEGER      | INPUT        | Length of CTYPE.    |
    +----------+----------+--------------+--------------+---------------------+
    | 32       | L_CSUPP  | INTEGER      | INPUT        | Length of CSUPP.    |
    +----------+----------+--------------+--------------+---------------------+
    | 33       | L_CTZONE | INTEGER      | INPUT        | Length of CTZONE.   |
    +----------+----------+--------------+--------------+---------------------+

    Stored data sets are without flags. For more information on flagging, see
    HECLIB documentation, Appendix C [1]_.

    Stored data is not compressed. For more information on data compression,
    see HECLIB documentation, Chapter 10 [1]_.

    References
    ----------

    The references below are formatted according to Chicago Manual of Style,
    16th Edition.

    .. [1] CEWRC-IWR-HEC. CPD-57.pdf. PDF. Davis, CA: US Army Corps of
       Engineers Institute for Water Resources Hydrologic Engineering Center,
       May 1991.
       Title: "HECLIB Volume 2: HECDSS Subroutines, Programmer's Manual"
       URL: http://www.hec.usace.army.mil/publications/ComputerProgramDocumentation/CPD-57.pdf
       Accessed: 2018-10-04

    .. [2] Steissberg, Todd. DSSExcel.xlam. Microsoft Excel XLAM. Davis, CA: US
       Army Corps of Engineers Institute for Water Resources Hydrologic
       Engineering Center, February 11, 2016. Title: "HEC-DSS MS Excel Data
       Exchange"

    """
    # Get data length from 'vals' list.
    nvals = len(vals)
    # Get DLL function for Storing Regular-Interval Time Series Data
    # (Extended Version).
    zsrtsc = getattr(dsslib, 'ZSRTSC_')
    # Initialize input declarations for function.
    # ???: Is initialization required?
    # <JAS 2018-10-02>
    zsrtsc.argparse = []
    # Set input declarations for function, mapped to DLL function variables.
    # ???: Why is CPATH multiplied by 800 instead of 80?
    # <JAS 2018-10-03>
    zsrtsc.argparse = [C.c_long*600,  # IFLTAB
                       C.c_char*800,  # CPATH
                       C.c_char*20,   # CDATE
                       C.c_char*4,    # CTIME
                       C.c_long,      # NVALS
                       C.c_long,      # DOUBLE
                       (C.c_float*nvals)(),   # SVALUES
                       (C.c_double*nvals)(),  # DVALUES
                       C.c_long,      # JQUAL
                       C.c_bool,      # LQUAL
                       C.c_char*8,    # CUNITS
                       C.c_char*8,    # CTYPE
                       C.c_double,    # COORDS
                       C.c_long,      # NCOORDS
                       C.c_long,      # ICDESC
                       C.c_long,      # NCDESC
                       C.c_char*80,   # CSUPP
                       C.c_long,      # ITZONE
                       C.c_char*30,   # CTZONE
                       C.c_long,      # IPLAN
                       C.c_long,      # JCOMP
                       C.c_float,     # BASEV
                       C.c_bool,      # LBASEV
                       C.c_bool,      # LDHIGH
                       C.c_long,      # NPREC
                       C.c_long,      # ISTAT
                       C.c_long,      # L_CPATH  25
                       C.c_long,      # L_CDATE  26
                       C.c_long,      # L_CTIME  27
                       C.c_long,      # L_CUNITS 28
                       C.c_long,      # L_CTYPE  29
                       C.c_long,      # L_CSUPP  30
                       C.c_long]      # L_CTZONE 31
    # Set return type of function to NoneType.
    zsrtsc.restype = None
    # Indicate that there are no flags in the data set.
    jqual = [0 for i in range(nvals)]
    lqual = False
    # Set coordinates.
    ncoords = len(coords)
    # Set length of Coordinate System Info.
    ncdesc = len(icdesc)
    # Set Time Zone Offset.
    # NOTE: Offset seems to be in minutes whereas milliseconds in HEC-DSS Vue.
    # TODO: Set itzone given ctzone.
    # <JAS 2018-10-05>
    if ctzone == 'PST':
        itzone = -420
    else:
        itzone = 0
    # Set parameters for no data compression.
    jcomp = 0
    basev = 0
    lbasev = False
    ldhigh = False
    nprec = 0
    # Prepare function inputs. Distiguish function inputs with "_" suffix.
    ifltab_ = ifltab
    cpath_ = cpath.encode('ascii')
    cdate_ = cdate.encode('ascii')
    ctime_ = ctime.encode('ascii')
    nvals_ = C.c_long(nvals)
    double_ = C.c_long(nvals)
    svalues_ = (C.c_float*nvals)(*vals)
    dvalues_ = (C.c_double*nvals)(*vals)
    jqual_ = (C.c_long*nvals)(*jqual)
    lqual_ = C.c_bool(lqual)
    cunits_ = cunits.encode('ascii')
    ctype_ = ctype.encode('ascii')
    coords_ = (C.c_double*3)(*coords)
    ncoords_ = C.c_long(ncoords)
    icdesc_ = (C.c_long*6)(*icdesc)
    ncdesc_ = C.c_long(ncdesc)
    csupp_ = csupp.encode('ascii')  # FUNCTION INPUT
    itzone_ = C.c_long(itzone)
    ctzone_ = ctzone.encode('ascii')
    iplan_ = C.c_long(iplan)
    jcomp_ = C.c_long(jcomp)
    basev_ = C.c_float(basev)
    lbasev_ = C.c_bool(lbasev)
    ldhigh_ = C.c_bool(ldhigh)
    nprec_ = C.c_long(nprec)
    istat_ = C.c_long()
    l_cpath_ = len(cpath_)
    l_cdate_ = len(cdate_)
    l_ctime_ = len(ctime_)
    l_cunits_ = len(cunits_)
    l_ctype_ = len(ctype_)
    l_csupp_ = len(csupp_)
    l_ctzone_ = len(ctzone_)
    # Pass variables to DLL function.
    zsrtsc(C.byref(ifltab_),
           cpath_,
           cdate_,
           ctime_,
           C.byref(nvals_),
           C.byref(double_),
           C.byref(svalues_),
           C.byref(dvalues_),
           C.byref(jqual_),
           C.byref(lqual_),
           cunits_,
           ctype_,
           C.byref(coords_),
           C.byref(ncoords_),
           C.byref(icdesc_),
           C.byref(ncdesc_),
           csupp_,
           C.byref(itzone_),
           ctzone_,
           C.byref(iplan_),
           C.byref(jcomp_),
           C.byref(basev_),
           C.byref(lbasev_),
           C.byref(ldhigh_),
           C.byref(nprec_),
           C.byref(istat_),
           l_cpath_,
           l_cdate_,
           l_ctime_,
           l_cunits_,
           l_ctype_,
           l_csupp_,
           l_ctzone_)
    return istat_.value


def create_catalog(ifltab, icunit, icdunt, inunit, cinstr, labrev, lsort): #, lcdcat, nrecs ):
    r"""
    Summary
    -------
    No documentation as of 2019-03-29.

    """
    zcat = getattr(dsslib, 'ZCAT_')
    zcat.argparse = []

    zcat.argparse = [C.c_long*600, C.c_long, C.c_long, C.c_long, C.c_char*800,
                     C.c_bool, C.c_bool, C.c_bool,
                     C.c_long]  # length of cinstr input string

    zcat.restype = None

    icunit = C.c_long(icunit)
    icdunt = C.c_long(icdunt)
    inunit = C.c_long(inunit)
    cinstr = (C.c_char*800)(*cinstr.encode('ascii'))
    labrev = C.c_bool(labrev)
    lsort = C.c_bool(lsort)
    lcdcat = C.c_bool()
    nrecs = C.c_long()

    #C.byref(icunit)
    zcat(C.byref(ifltab), C.byref(icunit), C.byref(icdunt), C.byref(inunit),
         cinstr, C.byref(labrev), C.byref(lsort), C.byref(lcdcat),
         C.byref(nrecs), len(cinstr))
    return([lcdcat.value, nrecs.value])


def get_catalog(fp):
    r"""
    Summary
    -------
    No documentation as of 2019-03-29.

    """
    # the DLL used in the Excel Add-in, which is the one we reference here for these
    # functions, creates a 'lock' file (*.dsk) when catalog operations are made;
    # once this lock file is created in a session, any *outside* modifications
    # to the *.dsc or *.dsd files will result in calls to open, create, or read
    # catalogs not working properly; this happens with the catalog functions in
    # the excel add-in as well;
    # the 'get_catalog' function can be called (and should work correctly) multiple
    # times within a session *as long as the *.dsc and *.dsd files are not modified
    # by another process*
    # as a precaution, we will try to delete any existing *.dsk files associated
    # with 'fp' upon calling this function - if the *.dsk file cannot be deleted,
    # a message will be returned indicating this condition

    dskfp = os.path.join(os.path.dirname(fp), os.path.basename(fp)[:-1]+'k')
    try:
        os.remove(dskfp)
    except:
        print("\n-------------------------------------------------------------\n"+ \
              "Could not remove the *dsk file - it is locked\nfor use by a previous " +\
              "call of this function \nor another system process (Excel add-in, I'm looking at you...)\n\n" +\
              "---DO NOT MODIFY THE CATALOG FILES (*.dsc, *.dsd) \n---WHILE RUNNING THIS PYTHON SESSION!!!\n" +\
              "================================================================\n")

    [ifltab, iostat] = open_dss(fp)  #open the DSS file to start
    if iostat != 0:
        print("couldn't open the dss file - exiting now...")
        return([None, None, None])

    dscfp = os.path.join(os.path.dirname(fp), os.path.basename(fp)[:-1]+'c')
    if os.path.exists(dscfp):
        if os.stat(dscfp).st_size==0:
            fortran_close_file(CATALOGUNIT)
            os.remove(dscfp)
    [lgenca, lopnca, lcatlg, lgencd, lopncd, lcatcd, nrecs] = open_catalog(fp, CATALOGUNIT, lgenca_in=True, lgencd_in=True)

    if not lopnca:
        print("Couldn't open catalog - something is not right...exiting..")
        return([None,None, None])
    else:
        NStdCatRecs = NCatalogedRecs = nrecs
        ValidCat = lcatlg

    # add some checks for number of records - ZINQR vs what's returned with open_catalog

    # check if a valid catalog was found
    if ValidCat:
        proceed=True
        print("Valid catalog found: %s \nNumber of records: %s" %(ValidCat, nrecs))
    else:
        print("Creating a catalog...")

        [lcatcd, nrecs] = create_catalog(ifltab, CATALOGUNIT, CONDUNIT, 0,'',False,True)
        #        Call vbaZCAT(CATALOGFILENUMBER, 0)
        #                          12,           0
        #                     TargetFile,  SourceFile
        #                     CatUnit,     SourceUnit
        print("Created condensed catalog: %s \nNumber of records: %s" %(lcatcd, nrecs))

        if nrecs>0:
            proceed=True
        else:
            proceed=False

    if proceed:

        if nrecs == 0:
            print("no records found in file...something is wrong again...exiting")
            fortran_close_file(CATALOGUNIT)
            fortran_close_file(CONDUNIT)
            close_dss(ifltab)
            return([None, None, None])
        else:
            [pathlist, lopnca] = read_catalog(lopnca, icunitin=CATALOGUNIT)
            close_dss(ifltab)
            fortran_close_file(CATALOGUNIT)
            fortran_close_file(CONDUNIT)
            return([pathlist, nrecs, lopnca])
    else:
        print("Couldn't create a catalog for some reason...exiting")
        close_dss(ifltab)
        fortran_close_file(CATALOGUNIT)
        fortran_close_file(CONDUNIT)
        return([None, None, None])


# %% Establish Code Body.
if __name__ == '__main__':
    msg = ('This module is intended to be imported for use into another '
           + 'module. It is not intended to be run as a __main__ file.')
    print(msg)
