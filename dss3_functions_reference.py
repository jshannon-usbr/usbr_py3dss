# -*- coding: utf-8 -*-
"""
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
"""

import ctypes as C
import os
import numpy as np

global CATALOGUNIT
global CONDUNIT

CATALOGUNIT=12
CONDUNIT = 13


dll_path = r'C:\Users\jmgilbert\01_Programs' +\
r'\Excel_DSS_Add-in_V3.3_Distribution-2015-01-28' +\
r'\Excel DSS Add-in V3.3.jmg\libraries\64-bit\heclib_x64.dll'
dsslib = C.cdll.LoadLibrary(dll_path)

def open_dss(fpi):
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
    zclose = getattr(dsslib, 'ZCLOSE_')
    zclose.argparse = [C.c_long*600]
    stdout = zclose(ifltab)
    return(stdout)
    
def dss_exists(fp):
    zfname = getattr(dsslib, 'ZFNAME_')
    zfname.argparse = [C.c_char*800, C.c_char*800, C.c_long, C.c_bool, C.c_long, C.c_long]
    # arguments: path in, path with extension if exists[out], length of returned name, 
    #            file exists logical, length of path in, length of path out
    zfname.restype = None
    plen = C.c_long()
    lexist = C.c_bool()
    new_fp = (C.c_char*800)()
    
    zfname(fp.encode('ascii'), new_fp, C.byref(plen),C.byref(lexist),len(fp.strip()),len(new_fp))
    
    #print(lexist.value)
    #print(new_fp)
    #print(plen)
    return([lexist.value, new_fp[0:].decode().strip()])
    
    
def read_regts(ifltab, cpath, cdate, ctime, nvalsi):
    zrrts = getattr(dsslib, 'ZRRTS_')
# Args: (IFLTAB(600)  | CPATH*80   | CDATE   | CTIME   | NVALS   | VALUES   |  CUNITS   | CTYPE   | IOFSET   | ISTAT)
# DTYPE:(INTEGER(600) | CHAR*80    | CHAR*20 | CHAR*4  | INTEGER | REAL(NVALS)| CHAR*8  | CHAR*8  | INTEGER  | INTEGER)
#        INPUT        | INPUT      | INPUT   | INPUT   | INPUT   | OUTPUT   | OUTPUT    | OUTPUT | OUTPUT    | OUTPUT)
    zrrts.argparse = [C.c_long*600, C.c_char*800, C.c_char*20, C.c_char*4, C.c_long, C.c_float, C.c_char*8, C.c_char*8, C.c_long, C.c_long, \
                      C.c_long, C.c_long, C.c_long, C.c_long, C.c_long]  # lengths of: cpath, cdate, ctime, cunits, ctype
    zrrts.restype = None
    
    nvals = C.c_long(nvalsi)  
    vals = (C.c_float*nvalsi)()
    cunits = (C.c_char*8)()
    ctype = (C.c_char*8)()
    iofset = C.c_long()
    istat = C.c_long()
    zrrts(ifltab, cpath.encode('ascii'), cdate.encode('ascii'), ctime.encode('ascii'), C.byref(nvals), C.byref(vals), cunits, ctype,C.byref(iofset), C.byref(istat), \
      len(cpath), len(cdate), len(ctime), len(cunits), len(ctype))
    
    return([nvals.value, np.float32(vals), cunits[0:].decode('utf-8').strip(), ctype[0:].decode('utf-8').strip(), iofset.value, istat.value])

def read_regtsd(ifltab, cpath, cdate, ctime, nvalsi, lgetdob_in=True):
    '''
    # Args: (IFLTAB(600)  | CPATH*80   | CDATE   | CTIME   | KVALS   | NVALS   | 
    #        INPUT        | INPUT      | INPUT   | INPUT   | INPUT   | OUTPUT   | OUTPUT    | OUTPUT | OUTPUT    | OUTPUT)
    # DTYPE:(INTEGER(600) | CHAR*80    | CHAR*20 | CHAR*4  | INTEGER | INTEGER | 
    
    # Args: LGETDOB | LFILDOB | sVALUES    | dVALUES       | JQUAL | LQUAL | LQREAD |  
    #       LOGICAL | LOGICAL | REAL(NVALS)| Double(NVALS) |INTEGER|LOGICAL| LOGICAL|

    #Args: CUNITS   | CTYPE   | CSUPP   | IOFSET   | JCOMP | ITZONE |CTZONE| COORDS |
    #       CHAR*8  | CHAR*8  | CHAR*80 | INTEGER  |INTEGER|INTEGER |CHAR*30 | DOUBLE |
    
    #Args: ICDESC |LCOORDS | ISTAT | L_CPATH | L_CDATE | L_CTIME | L_CUNITS | L_CTYPE | L_CSUPP | L_CTZONE)
    #      INTEGER|INTEGER |INTEGER| INTEGER | INTEGER | INTEGER | INTEGER  | INTEGER | INTEGER | INTEGER 
    '''
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
           C.byref(sVals), C.byref(dVals), C.byref(jqual), C.byref(lqual), C.byref(lqread),    #+5 = 13
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
    
    return([nvals.value, dVals, cunits[0:].decode('utf-8').strip(), ctype[0:].decode('utf-8').strip(),
            iofset.value, istat.value, csupp[0:].decode('utf-8'), coords_info, ctzone[0:].decode('utf-8'), jqual, lfildob, itzone])
    
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

def icdesc_to_dict(coords, icdesc):
    '''
    Converts the ICDESC list as returned from a DSS file to a dictionary with 
    keys as identifiers, according to DSSVue scheme
    '''
    
    coords_info = {}
    
    coords_info['X_Long'] = coords[0]
    coords_info['Y_Lat'] = coords[1]
    
    coordsys = icdesc[0]
    if coordsys==0:
        coords_info['CoordSys'] = ['No coordinates set', 0]
    elif coordsys==1:
        coords_info['CoordSys'] = ['LatLong', 1]
    elif coordsys==2:
        coords_info['CoordSys'] = ['State Plane/FIPS', 2]
    elif coordsys==3:
        coords_info['CoordSys'] = ['State Plane/ADS', 3]
    elif coordsys==4:
        coords_info['CoordSys'] = ['UTM', 4]
    elif coordsys==5:
        coords_info['CoordSys'] = ['Local/other', 5]
    else:
        coords_info['CoordSys'] = ['Hard to say, really', coordsys]
    
    coordID = icdesc[1]
    coords_info['CoordID'] = coordID
    
    datumUnits = icdesc[2]
    if datumUnits==0:
        coords_info['DatumUnit'] = ['Not specified', 0]
    elif datumUnits==1:
        coords_info['DatumUnit'] = ['English', 1]
    elif datumUnits==2:
        coords_info['DatumUnit'] = ['SI', 2]
    elif datumUnits==3:
        coords_info['DatumUnit'] = ['Decimal Degrees', 3]
    elif datumUnits==4:
        coords_info['DatumUnit'] = ['degrees, minutes, seconds', 4]
    else:
        coords_info['DatumUnit'] = ['Unknown', None]
    
    datum = icdesc[3]
    print("Datum from file is: %s" %datum)
    if datum==0:
        coords_info['Datum'] = ['Not specified', 0]
    elif datum==1:
        coords_info['Datum'] = ['NAD83', 1]
    elif datum==2:
        coords_info['Datum'] = ['NAD27', 2]
    elif datum==3:
        coords_info['Datum'] = ['WAS84', 3]
    elif datum==4:
        coords_info['Datum'] = ['WAS72', 4]
    elif datum==5:
        coords_info['Datum'] = ['local/other', 5]
    else:
        coords_info['Datum'] = ['Unknown', None]
    
    
    
    return(coords_info)
        
        
    
def open_catalog(fp, icunitin, lgenca_in=True, lgencd_in=False):   
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
    zopnca(fp.encode('ascii'), C.byref(icunit) , C.byref(lgenca), C.byref(lopnca), C.byref(lcatlg), \
           C.byref(icdunt), C.byref(lgencd), C.byref(lopncd), C.byref(lcatcd), C.byref(nrecs), len(fp.encode('ascii')))
    return([lgenca.value, lopnca.value, lcatlg.value, lgencd.value,lopncd.value, lcatcd.value, nrecs.value])
   
def read_catalog(lopnca, icunitin=12):
    # zrdcat - read pathnames from catalog
    # first - check to see if catalog is open - if not, prompt to open it
    if not lopnca:
        print("Catalog is not open - do this first using open_catalog function")
        return([None, None])
    
    zrdcat = getattr(dsslib, 'ZRDCAT_')
    # ARGS:  ICUNIT   |  LALL | IOUNIT  |CTAGS        | NDIM  | CPATHS | NPATHS  | NFOUND | *Len(ctags) | *Len(cpath)
    #DTYPES:  INT    | LOGICAL| INT     |CHAR(NDIM)*8 | INT   | 
    zrdcat.argparse= [C.c_long, C.c_bool, C.c_long, C.c_char, C.c_long, C.c_char, C.c_long, C.c_long, C.c_long, C.c_long ]
    icunit = C.c_long(icunitin)
    lall = C.c_bool(True)
    iounit = C.c_long(0)
    ndim = C.c_long(1) #nrecs  #C.c_long()
    ctags = (C.c_char*8)()
    cpaths = (C.c_char*80)()
    npaths = C.c_long()
    nfound = C.c_long(1)
    path_list = []
    

    while nfound.value ==1:
        zrdcat(C.byref(icunit), C.byref(lall), C.byref(iounit), ctags, C.byref(ndim),cpaths, C.byref(npaths), C.byref(nfound), len(ctags), len(cpaths) )
        #print(nfound.value)
        #print(cpaths.value)
        #print(ctags.value)
        nfoundv = nfound.value
        if nfoundv==1:
            path_list.append(cpaths.value.decode('utf-8').strip())
            #print(cpaths.value) 
    catstatus  = fortran_close_file(icunitin)  # close the catalog file
    if catstatus==0:
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
    
  
def write_regts(ifltab, cpath, cdate, ctime, nvalsi, valsi, cunits, ctype,_IPLAN=0):
    '''
    This function writes a regular time series, stored as a Python list, to a DSS file.
    The target DSS file to which the data should be written needs to be opened before 
    it can be written to. Opening the file (using `open_dss` function) will create
    a new file if the file doesn't already exist.
    INPUTS:
    -------
    ifltab: should be a long integer retrieved when opening a file
    cpath:  string of a properly formatted DSS pathname, should have 7 forward slashes /, 
            but no D-part (start date - this comes from the cdate argument)
    cdate:  a string representing the start date - should be of the form 
            "DDmonYYYY", e.g. "31Oct1999" for Oct 31, 1999
    ctime:  a string representing start time in hours since midnight on `cdate`,
            common practice is to use "2400" for monthly data
    nvalsi: number of values to be written, as an integer; the end date is defined
            by the start date, number of values, and the period length given in 
            the E-part of the variable path
    valsi:  a python list of length `nvalsi` to be written to the DSS file; missing values
            should be filled with -901.0 
    cunits: a string indicating units of the time series values
    ctype:  a string indicating the time period summary type that the values represent,
            commonly "PER-AVG" for period average
    _IPLAN:  (default = 0) - flag to indicate whether to write over existing data
            0 =  always write over existing data
            1 = only replace missing data flags in the record, as indicated by (-901)
            4 = if an input value is missing (-901), do not allow it to replace a
                non-missing value
    OUTPUTS:
    --------
    istat:   flag indicating if data was stored successfully (0), if all data was missing (-901),
             or other "fatal" errors(values >10, see HECLIB documentation for details for now...
             someday I'll include code to catch the errors and tell you something useful)
    ''' 
    zsrts = getattr(dsslib, 'ZSRTS_')
# Args: (IFLTAB(600)  | CPATH*80   | CDATE   | CTIME   | NVALS   | VALUES   |  CUNITS   | CTYPE   | IPLAN   | ISTAT)
# DTYPE:(INTEGER(600) | CHAR*80    | CHAR*20 | CHAR*4  | INTEGER | REAL(NVALS)| CHAR*8  | CHAR*8  | INTEGER  | INTEGER)
#        INPUT        | INPUT      | INPUT   | INPUT   | INPUT   | OUTPUT   | OUTPUT    | OUTPUT | OUTPUT    | OUTPUT)
    zsrts.argparse = []


    zsrts.argparse = [C.c_long*600, C.c_char*800, C.c_char*20, C.c_char*4, C.c_long, (C.c_float*nvalsi)(), C.c_char*8, C.c_char*8, C.c_long, C.c_long, \
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
    zsrts(C.byref(ifltab), cpath.encode('ascii'), cdate.encode('ascii'), ctime.encode('ascii'),
          C.byref(nvals), C.byref(vals), cunits.encode('ascii'), ctype.encode('ascii'),C.byref(IPLAN), C.byref(istat), \
          len(cpath), len(cdate), len(ctime), len(cunits), len(ctype))
    
    return([istat])
    
    
def create_catalog(ifltab, icunit, icdunt, inunit, cinstr, labrev, lsort): #, lcdcat, nrecs ):
    zcat = getattr(dsslib, 'ZCAT_')
    zcat.argparse = []

    zcat.argparse = [C.c_long*600,C.c_long, C.c_long, C.c_long, C.c_char*800, C.c_bool, C.c_bool, C.c_bool,
                     C.c_long]  # length of cinstr input string

    zcat.restype = None
    
    icunit = C.c_long(icunit)
    icdunt = C.c_long(icdunt)
    inunit= C.c_long(inunit)
    cinstr = (C.c_char*800)(*cinstr.encode('ascii'))
    labrev = C.c_bool(labrev)
    lsort = C.c_bool(lsort)
    lcdcat = C.c_bool()
    nrecs = C.c_long()
    

    #C.byref(icunit)
    zcat(C.byref(ifltab), C.byref(icunit) , C.byref(icdunt), C.byref(inunit), cinstr, \
           C.byref(labrev), C.byref(lsort), C.byref(lcdcat), C.byref(nrecs), len(cinstr))
    return([lcdcat.value, nrecs.value])
    
def get_catalog(fp):
    
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
    if iostat!=0:
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
        
        if nrecs==0:
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
    
    