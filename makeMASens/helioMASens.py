# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 16:23:15 2021

@author: mathewjowens
Ediited by Matthew Lang
"""

import os
import sys
import numpy as np
import astropy.time as ast
import makeMASens.downMASfc as dmf
import makeMASens.bravdaEns as bEns
import pandas as pd
import h5py

def observerReadH5(ephFile, body, times):
    # Function to extract the locations of a celestrial body/spacecraft specified in the
    # ephemeris HDF5 file at specific times.
    # Input:
    # ephFile: File containing ephemeris data to be extracted
    # body: Celestial object/spacecraft data to be extracted
    # times: Times (astropy.Time object) at which the data is to be extracted
    #
    # Output:
    # df: Directory containing the data at the required times. If there is no data at the times specified,
    #     the data is replaced with NaNs

    fileName = ephFile
    startDate = times.min()
    endDate = times.max()

    # Check body specified is a viable option
    possibleBodies = ['MERCURY', 'VENUS', 'EARTH', 'MARS',
                      'PSP', 'STERA', 'STERB', 'SO']

    if body.upper() in possibleBodies:
        body = body.upper()
    else:
        print('Warning, body requested was not recognised.')
        print(f'Only {possibleBodies} are valid.')
        print('Defaulting to Earth')
        body = ['EARTH']


    with h5py.File(fileName, "r") as f:
        # Extract variables

        # Generate time object to contain all dates in fileName
        allTime = ast.Time(f[body]['date'], format='mjd')
        allTime.format = 'datetime'

        # Define mask to restrict to only the relevant dates
        dtMask = (
                (allTime >= startDate) & (allTime <= endDate)
        )

        # Extract the variables at the required times
        bodyDates = ast.Time(allTime[dtMask])
        bodyRad = f[body]['radAU'][dtMask]
        bodyLat = f[body]['hgiLatDeg'][dtMask]
        bodyLon = f[body]['hgiLonDeg'][dtMask]

    # Interpolate the data to the required timeseries
    # Check that body dates are within times requested
    if len(bodyDates) > 0:
        bodyRadInterp = np.interp(times.jd, bodyDates.jd, bodyRad)
        bodyLatInterp = np.interp(times.jd, bodyDates.jd, bodyLat)
        bodyLonInterp = np.interp(times.jd, bodyDates.jd, bodyLon)

        # Fill position arrays with NaNs for [times < (bodyDates.min() - 2 days)] and
        #  for [times > (bodyDates.max() + 2 days)]
        twoDaySec = ast.TimeDelta(2 * 24 * 60 * 60, format='sec')
        lowMask = times.jd < (bodyDates.min() - twoDaySec).jd
        hiMask = times.jd > (bodyDates.max() + twoDaySec).jd
        interpMask = (lowMask | hiMask)

        if len(interpMask) > 0:
            bodyRadInterp[interpMask] = np.nan
            bodyLatInterp[interpMask] = np.nan
            bodyLonInterp[interpMask] = np.nan

            # Create dictionary to store required data
            df = pd.DataFrame(
                {'radAU': bodyRadInterp, 'hgiLatDeg': bodyLatInterp, 'hgiLonDeg': bodyLonInterp},
                index=times
            )
    else:
        df = pd.DataFrame(
            {'radAU': np.nan, 'hgiLatDeg': np.nan, 'hgiLonDeg': np.nan},
            index=times
        )

    return df


def makeMASens(
        crMJDFile, cr_start, cr_end, nMASens, nLon, ephemFile, masMapsDir, ensSaveDir,
        lat_rot_sigma=5 * np.pi / 180, lat_dev_sigma=2 * np.pi / 180,
        long_dev_sigma=2 * np.pi / 180, r_in=30
):
    # Function to generate MAS ensemble members and download MAS model runs if they don't exist
    # crMJDFile: File containing the start times of all Carrington Rotations
    # cr_start: Initial carrington rotation required
    # cr_end: Final Carrington Rotation required
    # nMASens: Number of MAS ensemble members required
    # nLon: Number of longitude points in each ens. member
    # ephemFile: Location of ephemeris file containing Earth's positional data
    # masMapsDir: Directory containing MAS maps (if they exist, if it doesn't, new directory will be made)
    # ensSaveDir: Directory to save ensembles to (if doesn't exist, new directory will be created)
    # lat_rot_sigma: The standard deviation of the Gaussian from which the rotational perturbation is drawn
    # lat_dev_sigma: The standard deviation of the Gaussian from which the linear
    #               latitudinal perturbation is drawn
    # lon_dev_sigma: The standard deviation of the Gaussian from which the linear
    #                longitudinal perturbation is drawn
    # r_in: Radial distance of speed map required (in rS, default is 30rS)

    # CR to be extracted
    # crMJDFile = 'C:/Users/mslan/Desktop/HUXtDA/CR_MJD/2039_2300CRMJDstart.csv'
    # cr_start = 2050
    # cr_end = 2230

    # If masMapsDir does not exist, create it
    if not os.path.isdir(masMapsDir):
        os.makedirs(masMapsDir)

    # If ensSaveDir does not exist, create it
    if not os.path.isdir(ensSaveDir):
        os.makedirs(ensSaveDir)

    # Read CR_MJD file to extract MJD start times for each Carrington Rotation
    crLines = []
    mjdLines = []
    with open(crMJDFile) as f:
        crMJDLines = f.readlines()

    for ln in crMJDLines:
        splitLn = (ln.strip()).split(',')

        crLines.append(float(splitLn[0]))
        mjdLines.append(float(splitLn[1]))

    # Generate MAS ensemble for required Carrington Rotation
    for nCR in range(cr_start, cr_end):
        cr = nCR
        print(f'CR = {cr}')

        reqIndex = next(
            x for x, val in enumerate(crLines) if val > cr
        )

        winCR = int(crLines[reqIndex - 1])

        if winCR != cr:
            print(f'winCR ({winCR}) != cr ({cr})')
            print('Exiting now')
            sys.exit()

        winCRStart = ast.Time(mjdLines[reqIndex - 1], format='mjd')
        winCREnd = ast.Time(mjdLines[reqIndex], format='mjd')

        # Make time-series of all times that need to be extracted from each MAS map
        deltaT = ast.TimeDelta(winCREnd - winCRStart, format='sec') / nLon

        timeGrid = ast.Time(
            [winCRStart + (deltaT * i)
             for i in range(nLon)]
        )
        crLonGrid = np.arange(360, 0, -360.0 / nLon)

        # Get the MAS maps
        vr_map, vr_lats, vr_longs = dmf.get_MAS_maps(cr, masMapsDir)

        # vr_map will be a single int if there is no data
        if not isinstance(vr_map, int):

            # Use the HUXt ephemeris data to get Earth lat over the CR
            earth = observerReadH5(ephemFile, 'Earth', timeGrid)

            # get Earth lat as a function of longitude (not time)
            E_lat = np.interp(vr_longs * 180 / np.pi, np.flipud(crLonGrid),
                              np.flipud(earth['hgiLatDeg']))

            # Convert E_lat to radians
            E_lat = E_lat * np.pi / 180.0

            # ==============================================================================
            # generate the input ensemble
            # ==============================================================================
            # generate the meshed grid
            phi, theta = np.meshgrid(vr_longs, vr_lats)

            vr_ensemble = bEns.generate_input_ensemble(
                phi, theta, vr_map, reflats=E_lat, Nens=nMASens,
                lat_rot_sigma=lat_rot_sigma, lat_dev_sigma=lat_dev_sigma, long_dev_sigma=long_dev_sigma
            )

            # resample the ensemble to 128 longitude bins
            vr128_ensemble = np.ones((nMASens, nLon))
            dphi = 2 * np.pi / nLon
            phi128 = np.linspace(dphi / 2, 2 * np.pi - dphi / 2, nLon)

            for i in range(0, nMASens):
                vr128_ensemble[i, :] = np.interp(
                    phi128, vr_longs, vr_ensemble[i, :]
                )

            # ==============================================================================
            # save the ensemble for use in  BRaVDA
            # ==============================================================================
            h5FileName = os.path.join(ensSaveDir, f'HelioMAS_CR{cr}_vin_ensemble.h5')
            h5f = h5py.File(h5FileName, 'w')
            h5f.create_dataset('Vin_ensemble', data=vr128_ensemble)

            outEnsTxtFile = open(f'{ensSaveDir}/vin_ensemble_CR{cr}.dat', 'w')
            np.savetxt(outEnsTxtFile, vr128_ensemble)
            outEnsTxtFile.close()

            h5f.attrs['lat_rot_sigma'] = lat_rot_sigma
            h5f.attrs['lat_dev_sigma'] = lat_dev_sigma
            h5f.attrs['long_dev_sigma'] = long_dev_sigma

            # this is used only to identify the source files.
            filepath = f'dmf.get_MAS_maps({cr}, {masMapsDir})'
            h5f.attrs['source_file'] = filepath
            h5f.attrs['r_in_rS'] = r_in
            h5f.attrs['Carrington_rotation'] = cr
            h5f.close()
