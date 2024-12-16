

# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 08:54:19 2023

@author: mathewjowens

A function to run BRaVDA on a given WSA map and generate HUXt input 
(i.e. speed as a function of Carr lon)

"""

import astropy.units as u
import numpy as np
import os as os
import glob
from astropy.io import fits
from scipy import interpolate
from sunpy.coordinates import sun

#from HUXt
import huxt_inputs as Hin
import huxt as H
import huxt_analysis as HA

#from HUXt_tools
import huxt_ensembles as Hens

#from BRaVDA
import startBravda



#load in WSA data
def _zerototwopi_(angles):
    """
    Function to constrain angles to the 0 - 2pi domain.

    Args:
        angles: a numpy array of angles.
    Returns:
            angles_out: a numpy array of angles constrained to 0 - 2pi domain.
    """
    twopi = 2.0 * np.pi
    angles_out = angles
    a = -np.floor_divide(angles_out, twopi)
    angles_out = angles_out + (a * twopi)
    return angles_out

#delete the output_dir so that teh correct posterior and prior are read in
def remove_non_empty_dir(dir_path):
    for root, dirs, files in os.walk(dir_path, topdown=False):
        for name in files:
            file_path = os.path.join(root, name)
            os.remove(file_path)
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(dir_path)
    
    
def readWSA(wsafilepath):
    """
    

    Parameters
    ----------
    wsafilepath : STRING
        Full filepath to WSA map, in FITS format

    Returns
    -------
    wsa_vr_map : 2D ARRAY, units of u.km/u.s
        Carrington map of WSA solar wind speeds
    vr_longs : 1D ARRAY, units of u.rad
        Carrington longitudes of the was_vr_map
    vr_lats : 1D ARRAY, units of u.rad
        Heliolatitudes of the was_vr_map
    wsa_map_date : ASTROPY TIME object
        date of the WSA map

    """
    
    
    #load the WSA data
    assert os.path.exists(wsafilepath)
    hdul = fits.open(wsafilepath)
    
    keys = hdul[0].header
    assert 'CARROT' in keys
    cr_num = hdul[0].header['CARROT']
    
    # different versions of WSA data have different keywords?
    if 'GRID' in keys:
        dgrid = hdul[0].header['GRID'] * np.pi / 180
    else:
        assert 'LONSTEP' in keys
        dgrid = hdul[0].header['LONSTEP'] * np.pi / 180
    
    # The map edge longitude is given by the CARRLONG variable.
    # This is 60 degrees from Central meridian (i.e. Earth Carrington longitude)
    carrlong = _zerototwopi_((hdul[0].header['CARRLONG']) * np.pi / 180)
    
    data = hdul[0].data
    vr_map_fits = data[1, :, :]
    
    hdul.close()
    
    # compute the Carrington map grids
    vr_long_edges = np.arange(0, 2 * np.pi + 0.00001, dgrid)
    vr_long_centres = (vr_long_edges[1:] + vr_long_edges[:-1]) / 2
    
    vr_lat_edges = np.arange(-np.pi / 2, np.pi / 2 + 0.00001, dgrid)
    vr_lat_centres = (vr_lat_edges[1:] + vr_lat_edges[:-1]) / 2
    
    vr_longs = vr_long_centres * u.rad
    vr_lats = vr_lat_centres * u.rad
    
    
    # rotate the maps so they are in the Carrington frame
    vr_map = np.empty(vr_map_fits.shape)
    
    for nlat in range(0, len(vr_lat_centres)):
        interp = interpolate.interp1d(_zerototwopi_(vr_long_centres + carrlong),
                                      vr_map_fits[nlat, :], kind="nearest",
                                      fill_value="extrapolate")
        vr_map[nlat, :] = interp(vr_long_centres)
    
    wsa_vr_map = vr_map * u.km / u.s
    cr_frac = cr_num + (2*np.pi - carrlong)/(2*np.pi)
    wsa_map_date = sun.carrington_rotation_time(cr_frac)
    
    return wsa_vr_map, vr_longs, vr_lats, wsa_map_date
    

def generate_HUXt_inputs_from_WSA_BRaVDA(wsafilepath, 
             r_min = 21.5*u.solRad, assim_end_date = None, obs = ['OMNI'],
             deacc = True, huxt_bravda_scale = 0.95,  
             smooth_posterior = True, smooth_val = 2,
             bravda_dir = None, runname = 'temp',
             Nens = 500, 
             lat_rot_sigma = 10*np.pi/180 *u.rad, 
             lat_dev_sigma = 0*np.pi/180 *u.rad,  
             long_dev_sigma = 10*np.pi/180 *u.rad
             ):
    
    """
    A function to load a WSA map, produce the ensemble required by BRaVDA,
    run BRaVDA and process the output for use by HUXt
    
    Args:
        wsafilepath : STRING. The full filepath to the WSA map
        r_min : FLOAT, units of u.solRad. The heliocentric distance of the WSA map, typically 21.5
        assim_end_date : DATETIME. Time of the last in situ obs to assimilate. If none, 
                            defaults to the WSA map date.
        obs : LIST. Which in situ data to assimilate, typically one or more of ['OMNI', 'STERA', 'STERB']
        deacc : BOOL. whether or not to deaccelerate the WSA speeds from 1 au to r_min
        huxt_bravda_scale : FLOAT. how much to rescale the BRaVDA output for use in HUXT (1 = no change)
        smooth_posterior: BOOL. WHether to smooth the BRaVDA output for use with HUXt
        smooth_val: FLOAT. Level of Gaussian smoothing
        bravda_dir: STRING. Location of BRaVDA install. Set to None for this file location
        runname: STRING. Name for new run dir. Will delete whatever's already in there.
        Nens: INT. Number of ensemble members for BRaVDA covariance matrix
        lat_rot_sigma : FLOAT, units of u.rad. Latitudinal rotational perturbation for ensemble generations
        lat_dev_sigma : FLOAT, units of u.rad. Latitudinal DC perturbation for ensemble generations
        long_rot_sigma : FLOAT, units of u.rad. Longitudinal rotational perturbation for ensemble generations  
        
    Returns:
       prior_wsalongs : FLOAT ARRAY, units of u.km/u.s. The BRaVDA prior, in Carrington longitude
       post_wsalongs : FLOAT ARRAY, units of u.km/u.s. The BRaVDA posterior, with processing, in Carrington longitude
       cr: INT. The carrington rotation number
       cr_lon_init: FLOAT, units of u.rad. Carrington longitude of Earth
       r_min: FLOAT, units of u.solRad. The heliocentric distance of the WSA map, typically 21.5
       E_lat: FLOAT, units of u.rad. Heliolaitude of Earth.
    """

    # assume the working dir is the same as teh WSA file
    # Extract the filename without extension
    data_dir =  os.path.dirname(wsafilepath)
    
    #bravda ensemble directory
    if bravda_dir is None:
        bravda_dir = os.path.abspath(os.path.dirname(__file__))
    bravda_ens_dir = os.path.join(bravda_dir,'masEns')
    output_dir = os.path.join(data_dir,runname)
    


    # Delete the directory and all its contents
    try:
        remove_non_empty_dir(output_dir)
        print(f"Directory {output_dir} and all its contents have been deleted successfully.")
    except OSError as e:
        print(f"Error: {e.strerror} - {e.filename}")
    
    #check the output directory exists
    if not os.path.exists(output_dir):
        # Create the directory
        os.makedirs(output_dir)
    
    # read in the WSA data
    wsa_vr_map, vr_longs, vr_lats, wsa_map_date = readWSA(wsafilepath)

    
    
    if deacc:
        #deaccelerate the WSA map from 1-AU calibrated speeds to expected 21.5 rS values
        for nlat in range (1, len(vr_lats)):
            wsa_vr_map[nlat,:], lon_temp = Hin.map_v_inwards(wsa_vr_map[nlat,:], 215*u.solRad, 
                                                     vr_longs, r_min)
            
      
    # check the end assimilation date. if none is provided, use the WSA map date
    if assim_end_date is None:
        assim_end_date = wsa_map_date.to_datetime()
            
    
    #get the huxt params for the start time
    cr, cr_lon_init = Hin.datetime2huxtinputs(assim_end_date)
    
    
    
    #Use the HUXt ephemeris data to get Earth lat over the CR
    #========================================================
    dummymodel = H.HUXt(v_boundary=np.ones((128))*400* (u.km/u.s), simtime=27.27*u.day, 
                       cr_num= cr, cr_lon_init = cr_lon_init, 
                       lon_out=0.0*u.deg,
                       r_min=r_min)
    
    #retrieve a bodies position at each model timestep:
    earth = dummymodel.get_observer('earth')   
    reflats = np.interp(vr_longs,earth.lon_c,earth.lat_c)
    
    
    # generate WSA ensemble for the DA and put it in the BRaVDA dir 
    phi, theta = np.meshgrid(vr_longs, vr_lats, indexing = 'xy')
    
    vr_ensemble = Hens.generate_input_ensemble(phi, theta, wsa_vr_map, 
                                reflats, Nens = Nens, 
                                lat_rot_sigma = lat_rot_sigma, 
                                lat_dev_sigma = lat_dev_sigma,
                                long_dev_sigma = long_dev_sigma)
    
    #resample the ensemble to 128 longitude bins
    vr128_ensemble = np.ones((Nens,128))  
    dphi = 2*np.pi/128
    phi128 = np.linspace(dphi/2, 2*np.pi - dphi/2, 128)
    for i in range(0, Nens):
        vr128_ensemble[i,:] = np.interp(phi128,
                      vr_longs.value,vr_ensemble[i,:])
    
    
    #also save vr128 to a .dat file for use in BRaVDA
    #outEnsTxtFile = open(f'{bravda_ens_dir}customensemble.dat', 'w')
    outEnsTxtFile = os.path.join(bravda_ens_dir, 'customensemble.dat')
    np.savetxt(outEnsTxtFile, vr128_ensemble)
    #outEnsTxtFile.close()


    #==============================================================================
    #==============================================================================
    #Run startBravda
    #==============================================================================
    #==============================================================================
    
    
    startBravda.bravdafunction(assim_end_date, obsToAssim = obs, usecustomens = True,
                               runoutputdir = output_dir, plottimeseries = False,
                               corona = 'WSA')
    
    #fmjd = int(Time(assim_end_date).mjd)
    
    #read in the posterior solution to blend with the WSA map   
    #BRaVDA files are labelled with teh start time, 27 or 28 days previous. helpful
    post_file = glob.glob(os.path.join(output_dir,'posterior','posterior_MJDstart*'))[0]
    prior_file = glob.glob(os.path.join(output_dir,'prior','prior_MJDstart*'))[0]
    
    posterior = np.loadtxt(post_file)
    prior = np.loadtxt(prior_file)
    
    #the inner boundary value is given as a function of time from the inition time
    prior_vt = prior[0,:]
    post_vt = posterior[0,:]
    
    #smooth the posterior for use in HUXt
    post_vt_unsmooth = post_vt
    if smooth_posterior:
       from scipy.ndimage import gaussian_filter1d
       post_vt = gaussian_filter1d(post_vt_unsmooth, sigma=smooth_val, mode='wrap') 
       # plt.figure()
       # plt.plot(post_vt_unsmooth, label='posterior unsmoothed')
       # plt.plot(post_vt, label='posterior smoothed, sigma = ' + str(smooth_val))
       # plt.legend()
       
    #scale the speeds for use in HUXt
    post_vt = post_vt * huxt_bravda_scale
    
    #convert to function of longitude
    post_vlong = np.flipud(post_vt)
    prior_vlong = np.flipud(prior_vt)
    #post_vlong = post_inner
    #find the associated carr_longs
    cr, cr_lon_init = Hin.datetime2huxtinputs(assim_end_date)
    post_carrlongs = _zerototwopi_(phi128 + cr_lon_init.value)
    prior_carrlongs = _zerototwopi_(phi128 + cr_lon_init.value)
    
    
    
    
    #interpolate the posterior at the inner boundary to the WSA CarrLong grid
    interp = interpolate.interp1d(post_carrlongs,
                                  post_vlong, kind="nearest",
                                  fill_value="extrapolate")
    post_wsalongs = interp(vr_longs.value) * u.km/u.s
    
    interp = interpolate.interp1d(prior_carrlongs,
                                  prior_vlong, kind="nearest",
                                  fill_value="extrapolate")
    prior_wsalongs = interp(vr_longs.value) * u.km/u.s
    
    E_lat =  np.nanmean(reflats)
    
    
    
    # Delete the directory and all its contents - clean up
    try:
        remove_non_empty_dir(output_dir)
        print(f"Directory {output_dir} and all its contents have been deleted successfully.")
    except OSError as e:
        print(f"Error: {e.strerror} - {e.filename}")

    return prior_wsalongs, post_wsalongs, cr, cr_lon_init, r_min, E_lat


# <codecell> Test with HUXt

runnow = True

if runnow:
    obs = ['OMNI', 'STERA', 'STERB']
    wsafilepath = os.path.join(os.environ['DBOX'], 'python_repos', 'HUXt_tools', 
                               'data', 'models%2Fenlil%2F2022%2F11%2F9%2F0%2Fwsa.gong.fits')
    
    prior_wsalongs, post_wsalongs, cr, cr_lon_init, r_min, E_lat = generate_HUXt_inputs_from_WSA_BRaVDA(wsafilepath,
                                    obs = obs)
    
    
    
    #to plot the assimilation window, use cr_num = cr-1. To plot a forecast, use cr_num = cr
    model = H.HUXt(v_boundary=prior_wsalongs, 
                   cr_num=cr-1, cr_lon_init=cr_lon_init, latitude = E_lat,
                   simtime=27.27*u.day, dt_scale=4, r_min = r_min, lon_out=0*u.rad)
    model.solve([])
    
    #plot it
    fig, ax = HA.plot_earth_timeseries(model, plot_omni = True)
    
    
    model = H.HUXt(v_boundary=post_wsalongs, 
                   cr_num=cr, cr_lon_init=cr_lon_init, latitude = E_lat,
                   simtime=27.27*u.day, dt_scale=4, r_min = r_min, lon_out=0*u.rad)
    model.solve([])
    
    #plot it
    fig, ax = HA.plot_earth_timeseries(model, plot_omni = True)
        
