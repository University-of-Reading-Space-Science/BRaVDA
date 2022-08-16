"""
Created on Thu Jun 11 12:50:49 2020

@author: mathewjowens
Editted by Matthew Lang
"""
import httplib2
import urllib.request
import os
from pyhdf.SD import SD, SDC
import numpy as np


# Get MAS data from MHDweb
def getMASboundaryconditions(
        cr=np.NaN, dirLoc='', observatory='',
        runtype='', runnumber=''
):
    """
    A function to grab the  Vr and Br boundary conditions from MHDweb. An order
    of preference for observatories is given in the function. Checks first if
    the data already exists in the BRaVDA MAS ensemble folder

    Parameters
    ----------
    cr : INT
        Carrington rotation number
    dirLoc : str
        Directory to which MAS conditions should be written
    observatory : STRING
        Name of preferred observatory (e.g., 'hmi','mdi','solis',
        'gong','mwo','wso','kpo').
        Empty if no preference and automatically selected
    runtype : STRING
        Name of preferred MAS run type (e.g., 'mas','mast','masp').
        Empty if no preference and automatically selected
    runnumber : STRING
        Name of preferred MAS run number (e.g., '0101','0201').
        Empty if no preference and automatically selected

    Returns
    -------
    flag : INT
        1 = successful download. 0 = files exist, -1 = no file found.

    """

    assert (not np.isnan(cr))

    # the order of preference for different MAS run results
    overwrite = False
    if not observatory:
        observatories_order = [
            'hmi', 'mdi', 'solis', 'gong', 'mwo', 'wso', 'kpo'
        ]
    else:
        observatories_order = [str(observatory)]
        # if the user wants a specific observatory,
        # overwrite what's already downloaded
        overwrite = True

    if not runtype:
        runtype_order = ['masp', 'mas', 'mast']
    else:
        runtype_order = [str(runtype)]
        overwrite = True

    if not runnumber:
        runnumber_order = ['0201', '0101']
    else:
        runnumber_order = [str(runnumber)]
        overwrite = True

    # get the BRaVDA MAS boundary condition directory
    _boundary_dir_ = dirLoc
    # dirs['boundary_conditions']

    # example URL:
    #   'http://www.predsci.com/data/runs/cr2010-medium/
    #    mdi_mas_mas_std_0101/helio/vr_r0.hdf'
    heliomas_url_front = 'http://www.predsci.com/data/runs/cr'
    heliomas_url_end = '_r0.hdf'

    vrfilename = f'HelioMAS_CR{int(cr)}_vr{heliomas_url_end}'

    if ((
            not os.path.exists(
                os.path.join(_boundary_dir_, vrfilename)
            )
    ) or overwrite):
        # check if the files already exist
        # Search MHDweb for a HelioMAS run, in order of preference
        h = httplib2.Http()
        foundfile = False
        for masob in observatories_order:
            for masrun in runtype_order:
                for masnum in runnumber_order:
                    urlbase = (
                        f'{heliomas_url_front}{int(cr)}-medium/'
                        f'{masob}_{masrun}_mas_std_{masnum}/helio/'
                    )
                    url = f'{urlbase}vr{heliomas_url_end}'
                    # print(url)

                    # see if this vr file exists
                    resp = h.request(url, 'HEAD')
                    if int(resp[0]['status']) < 400:
                        foundfile = True
                        # print(url)

                    # exit all the loops - clumsy, but works
                    if foundfile:
                        break
                if foundfile:
                    break
            if foundfile:
                break

        if not foundfile:
            print('No data available for given CR and observatory'
                  ' preferences')
            return -1

        # download the vr files
        print(f'Downloading from: {urlbase}')
        print(f' to {os.path.join(_boundary_dir_, vrfilename)}')
        urllib.request.urlretrieve(
            f'{urlbase}vr{heliomas_url_end}', os.path.join(_boundary_dir_, vrfilename)
        )

        return 1
    else:
        print(f'Files already exist for CR{int(cr)}')
        return 0


def readMASvr(cr, dirLoc):
    """
    A function to read in the MAS boundary conditions for a given CR

    Parameters
    ----------
    cr : INT
        Carrington rotation number
    dirLoc : str
        Directory to which MAS conditions should be written

    Returns
    -------
    MAS_vr : NP ARRAY (NDIM = 2)
        Solar wind speed at 30rS, in km/s
    MAS_vr_Xa : NP ARRAY (NDIM = 1)
        Carrington longitude of Vr map, in rad
    MAS_vr_Xm : NP ARRAY (NDIM = 1)
        Latitude of Vr as angle down from N pole, in rad
    """

    # get the boundary condition directory
    _boundary_dir_ = dirLoc
    # dirs['boundary_conditions']

    # create the filenames
    heliomas_url_end = '_r0.hdf'
    vrfilename = f'HelioMAS_CR{int(cr)}_vr{heliomas_url_end}'

    filepath = os.path.join(_boundary_dir_, vrfilename)
    assert os.path.exists(filepath)
    # print(os.path.exists(filepath))

    file = SD(filepath, SDC.READ)
    # print(file.info())
    # datasets_dic = file.datasets()
    # for idx,sds in enumerate(datasets_dic.keys()):
    #     print(idx,sds)

    sds_obj = file.select('fakeDim0')  # select sds
    MAS_vr_Xa = sds_obj.get()  # get sds data
    sds_obj = file.select('fakeDim1')  # select sds
    MAS_vr_Xm = sds_obj.get()  # get sds data
    sds_obj = file.select('Data-Set-2')  # select sds
    MAS_vr = (sds_obj.get()).transpose()  # get sds data

    # convert from model to physicsal units
    MAS_vr = MAS_vr * 481.0  # * u.km/u.s
    MAS_vr_Xa = MAS_vr_Xa  # * u.rad
    MAS_vr_Xm = MAS_vr_Xm  # * u.rad

    return MAS_vr, MAS_vr_Xa, MAS_vr_Xm


def get_MAS_long_profile(cr, dirLoc, nLon=128, lat=0.0):
    """
    a function to download, read and process MAS output
    to provide HUXt boundary conditions at a given latitude

    Parameters
    ----------
    cr : INT
        Carrington rotation number
    dirLoc : str
        Directory to which MAS conditions should be written
    nLon: Integer
        Number of longitude gridspaces
    lat : FLOAT
        Latitude at which to extract the longitudinal profile,
        measure up from equator

    Returns
    -------
    vr_in : NP ARRAY (nLon)
        Solar wind speed as a function of Carrington longitude
        at solar equator.
        Interpolated to HUXt longitudinal resolution. In km/s
    """

    assert ((not np.isnan(cr)) and (cr > 0))
    assert (lat >= -90.0)
    assert (lat <= 90.0)

    # convert angle from equator to angle down from N pole
    ang_from_N_pole = (np.pi / 2.0) - (lat * ((2 * np.pi) / 180.0))

    # check the data exist, if not, download them
    flag = getMASboundaryconditions(cr, dirLoc)
    # getMASboundaryconditions(cr,observatory='mdi')
    assert (flag > -1)

    # read the HelioMAS data
    MAS_vr, MAS_vr_Xa, MAS_vr_Xm = readMASvr(cr, dirLoc)

    # extract the value at the given latitude
    vr = np.ones(len(MAS_vr_Xa))
    for i in range(0, len(MAS_vr_Xa)):
        vr[i] = np.interp(
            ang_from_N_pole, MAS_vr_Xm, MAS_vr[i][:]
        )

    # now interpolate on to the HUXt longitudinal grid
    dLon = (2 * np.pi) / nLon
    halfDLon = dLon / 2.0
    lon = np.arange(halfDLon, 2 * np.pi - halfDLon, dLon)

    vr_in = np.interp(
        lon, MAS_vr_Xa, vr
    )  # *u.km/u.s

    return vr_in


def get_MAS_maps(cr, dirLoc):
    """
    a function to download, read and process MAS output to
    provide HUXt boundary
    conditions as lat-long maps, along with angle from equator for the maps
    maps returned in native resolution, not HUXt resolution

    Parameters
    ----------
    cr : INT
        Carrington rotation number
    dirLoc : str
        Directory to which MAS conditions should be written


    Returns
    -------
    vr_map : NP ARRAY
        Solar wind speed as a Carrington longitude-latitude map.
        In km/s
    vr_lats :
        The latitudes for the Vr map, in radians from trhe equator
    vr_longs :
        The Carrington longitudes for the Vr map, in radians
    """

    assert ((not np.isnan(cr)) and (cr > 0))

    # check the data exist, if not, download them
    flag = getMASboundaryconditions(cr, dirLoc)
    # getMASboundaryconditions(cr,observatory='mdi')
    assert (flag > -1)

    # read the HelioMAS data
    MAS_vr, MAS_vr_Xa, MAS_vr_Xm = readMASvr(
        cr, dirLoc)

    vr_map = MAS_vr

    # convert the lat angles from N-pole to equator centred
    vr_lats = (np.pi / 2) - MAS_vr_Xm

    # flip lats, so they're increasing in value
    vr_lats = np.flipud(vr_lats)
    vr_map = np.fliplr(vr_map)

    vr_longs = MAS_vr_Xa

    return vr_map, vr_lats, vr_longs
