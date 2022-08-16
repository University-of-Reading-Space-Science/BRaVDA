# Import necessary packages
import numpy as np
import time
from scipy import optimize
import os
import shutil
import sys

import bravdaMethodsExport as bme
from makeMASens import helioMASens


def runBravDA(configFile, huxVarFile, outputDir, obsToUse, setupOfR,
              initDate, noOfWindows, nMASens, locRad, gTol, makePlots):
    # Initialise timer
    start_time = time.time()

    ##################################################
    # Find current directory
    currentDir = os.path.dirname(os.path.abspath(__file__))

    #########################################################
    # Read in file containing all constants/directory names
    # ((should be kept in same directory as this script))
    #########################################################
    # Open config file that contains all the directories and file locations
    print(f'Config file: {configFile}')

    with open(configFile, 'r') as fConfig:
        configLines = fConfig.readlines()

    # Set up variables containing directories/file names
    fileCRstartMJD = os.path.join(currentDir, configLines[1].strip())
    downMASdir = os.path.join(currentDir, configLines[3].strip(), '')
    dirMASens = os.path.join(currentDir, configLines[5].strip(), '')
    steraLocFile = os.path.join(currentDir, configLines[7].strip())
    sterbLocFile = os.path.join(currentDir, configLines[9].strip())
    earthLocFile = os.path.join(currentDir, configLines[11].strip())
    fileObsA = os.path.join(currentDir, configLines[13].strip())
    fileObsB = os.path.join(currentDir, configLines[15].strip())
    fileObsACE = os.path.join(currentDir, configLines[17].strip())

    #########################################################
    # Read in file containing required inputs for the HUX solar wind
    # propagation model
    # ((should be kept in same directory as this script))
    #########################################################
    # Open parameter file that contains inputs for
    print(f'HUX var file: {huxVarFile}')

    with open(huxVarFile, 'r') as fhv:
        huxVarLines = fhv.readlines()

    noOfLonPoints = int(huxVarLines[1].strip())
    innerRad = float(huxVarLines[3].strip())
    outerRad = float(huxVarLines[5].strip())
    deltaR = float(huxVarLines[7].strip())
    alpha = float(huxVarLines[9])
    rH = float(huxVarLines[11])

    print(f'Output directory: {outputDir}')

    # Make output directory if it doesn't already exist
    if not os.path.isdir(outputDir):
        os.makedirs(outputDir)

    # Copy config and huxVar files used for this run to the
    # output directory so the user can see what parameters
    # were used to generate the most *recent* output in that output dir
    shutil.copy2(configFile, outputDir)
    shutil.copy2(huxVarFile, outputDir)

    # Write file containing all additional inputs to the output directory
    inWriteFile = os.path.join(outputDir, 'inputs.dat')
    with open(inWriteFile, 'w') as iwf:
        iwf.write(f'outputDir = {outputDir}\n')
        iwf.write(f'obsToUse = {obsToUse}\n')
        iwf.write(f'setupOfR = {setupOfR}\n')
        iwf.write(f'initDate = {initDate}\n')
        iwf.write(f'noOfConsecWindows = {noOfWindows}\n')
        iwf.write(f'noOfMASens = {nMASens}\n')
        iwf.write(f'locRad = {locRad}\n')
        iwf.write(f'gTol = {gTol}\n')

    # Make all necessary directories in output directory to save files
    # List containing all folders to make
    dirsToMake = ['meanMAS', 'prior', 'posterior', 'ACE', 'STERA', 'STERB']
    for dirName in dirsToMake:
        dirPath = os.path.join(outputDir, dirName, '')

        # Check if dir exists, if not, make it
        if not os.path.isdir(dirPath):
            os.makedirs(dirPath)

    # Make all necessary directories in output directory to save plots
    if makePlots:
        dirsToMake = ['swOver27d']
        plotDir = os.path.join(outputDir, 'plots', '')

        # Make directory to hold plots
        if not os.path.isdir(plotDir):
            os.makedirs(plotDir)

        # Make subdirectories in plots directory
        for dirName in dirsToMake:
            dirPath = os.path.join(plotDir, dirName, '')

            # Check if dir exists, if not, make it
            if not os.path.isdir(dirPath):
                os.makedirs(dirPath)

    #########################################
    # Define constants
    #########################################
    # Set Initial/Final Carrington Rotations to be run
    initMJD = bme.dateToMJD(initDate)
    finalMJD = initMJD + (27 * noOfWindows)

    print(f'Initial MJD: {initMJD}')
    print(f'Final MJD  : {finalMJD}')
    print(f'No. Windows: {noOfWindows}')

    # Longitudinal constants
    deltaPhiDeg = float(360.0 / noOfLonPoints)
    deltaPhi = deltaPhiDeg * (np.pi / 180.0)
    totalLonPoints = noOfLonPoints * noOfWindows

    # Radial constants (rS = solar radius)
    rS = 700000  # 695700

    # Convert innerRad, outerRad and deltaR to km
    innerRadRs = innerRad * rS
    outerRadRs = outerRad * rS
    deltaRrs = deltaR * rS
    r = np.arange(innerRadRs, outerRadRs, deltaRrs)
    noOfRadPoints = len(r)

    # Constants for solar wind propagation model
    solarRotFreq = (2 * np.pi) / (25.38 * 24 * 3600.0)

    # Convert rH to km from rS
    rH = rH * rS

    #######################################
    # Read in which obs. to assimilate
    #######################################
    # Observations assimilated
    # Input a string of A, B and C where
    # A=STERA, B=STERB and C='ACE' eg. for
    # all assim, input obsToUse='ABC'
    # obsToUse = str(vSWDALines[33].strip()).upper()
    # Order the characters alphabetically and make upper-case
    obsToUse = (obsToUse.strip()).upper()
    obsToUse = ''.join(sorted(obsToUse))
    print(f'Observations to be assimilated: {obsToUse}')

    ########################################################
    # Initialise arrays for use later in script
    ########################################################
    # Observation vectors
    yAllA = 9999 * np.ones(totalLonPoints)
    yAllB = 9999 * np.ones(totalLonPoints)
    yAllC = 9999 * np.ones(totalLonPoints)

    # Vectors to hold solar wind speeds
    vPrior = np.zeros((noOfWindows, noOfRadPoints, noOfLonPoints))
    vPosterior = np.zeros((noOfWindows, noOfRadPoints, noOfLonPoints))
    vMASMean = np.zeros((noOfWindows, noOfRadPoints, noOfLonPoints))

    # If user wishes to plot all windows consecutively, this array holds the solar wind speeds sequentially
    vPriorAll = np.zeros((noOfRadPoints, totalLonPoints))
    vPosteriorAll = np.zeros((noOfRadPoints, totalLonPoints))

    # Initialise arrays to store RMSEs of the solar wind speeds at obs. locations
    RMSESTERAPrior = np.zeros(noOfWindows)
    RMSESTERAPosterior = np.zeros(noOfWindows)
    RMSESTERAMASMean = np.zeros(noOfWindows)

    RMSESTERBPrior = np.zeros(noOfWindows)
    RMSESTERBPosterior = np.zeros(noOfWindows)
    RMSESTERBMASMean = np.zeros(noOfWindows)

    RMSEACEPrior = np.zeros(noOfWindows)
    RMSEACEPosterior = np.zeros(noOfWindows)
    RMSEACEMASMean = np.zeros(noOfWindows)

    ##################################################################################
    # Make MJD and initial ensemble files and store them in an array to be read in later
    ##################################################################################
    fileMJD = []  # Variable to hold MJD file names
    initMJDCR = []  # Variable to hold CR of initial MJD, so we can extract from the correct MAS ens file
    MJDLines = []
    CRLines = []

    currMJD = []  # Variables to hold each window's MJD and CR starting points
    currCR = []
    midMJD = []  # Variable to hold mid-point MJD of each 27-day period

    # Read in table of MJD values and Carrington rotations
    with open(fileCRstartMJD, 'r') as fMJDstart:
        fMJDLines = fMJDstart.readlines()

    # Read in MJD values and CR values and store them in an array
    for ln in range(len(fMJDLines)):
        splitLine = (fMJDLines[ln].strip()).split(',')

        CRLines.append(int(splitLine[0]))
        MJDLines.append(float(splitLine[1]))

    ###########################################################################################
    # Find MJD starting/mid-points for all windows, find which Carrington Rotation they're in
    # and make MJD files for each longitude point in each 27 day period
    ###########################################################################################
    for w in range(noOfWindows):
        # Calculate current window's initial MJD
        currMJD.append(initMJD + (27 * w))

        # Calculate current window's mid-point MJD
        midMJD.append(currMJD[-1] + 13.5)

        # Calculate which CR currMJD is in
        # Find index where currMJD first exceeds MJD lines, currMJD will be in the previous line's CR
        reqIndex = next(x for x, val in enumerate(MJDLines) if val > currMJD[-1])
        currCR.append(CRLines[reqIndex - 1])
        initMJDCR.append(MJDLines[reqIndex - 1])

        # Store name of required files for creating/reading later
        dirMJDpath = os.path.join(outputDir, 'MJDfiles')
        fileMJDpath = os.path.join(dirMJDpath, f'MJD_{int(currMJD[-1])}.dat')
        fileMJD.append(fileMJDpath)

        # Check if MJD file exists
        # If file does not exist, make it
        if not os.path.isfile(fileMJDpath):
            # Check if directory exists, if not make it
            if not os.path.isdir(dirMJDpath):
                os.makedirs(dirMJDpath)

            # Make MJD file, downloading new data if necessary
            bme.makeMJDfile(currMJD[-1], noOfLonPoints, fileMJDpath, daySolarRot=27)

    #####################################################################################
    # Check existence of and make initial ensemble files if required for each window
    #####################################################################################
    # Directory to hold solar wind ens files specific to this run
    filesVrEns = []  # Variable to hold required initial solar wind speed ensemble file names for this run
    for w in range(noOfWindows):
        filesVrEns.append(f'{dirMASens}vin_ensemble_CR{currCR[w]}.dat')

        if not os.path.isfile(filesVrEns[w]):
            # If MAS ensemble file does not exist, make it
            print(f'Generating MAS ensembles for CR {currCR[w]}...')

            # Generate path to ephemerisMSL.hdf5 file (should be in makeMASens, don't move)
            ephemFile = os.path.join('makeMASens', 'ephemerisMSL.hdf5')
            helioMASens.makeMASens(
                fileCRstartMJD, currCR[w], currCR[w] + 1, nMASens, noOfLonPoints,
                ephemFile, downMASdir, dirMASens,
                lat_rot_sigma=5 * np.pi / 180, lat_dev_sigma=2 * np.pi / 180,
                long_dev_sigma=2 * np.pi / 180, r_in=innerRad
            )
        else:
            print(f'MAS ensembles exist for CR {currCR[w]}')

    ############################################################################
    # Check if mid-points file exists for these windows, if not, make it
    ############################################################################
    mjdCRFile = os.path.join(
        outputDir, 'MJDFiles', f'MJDMidPoints_{int(currMJD[0])}_{int(currMJD[-1])}.dat'
    )
    if not os.path.isfile(mjdCRFile):
        with open(mjdCRFile, 'w') as mjdMidFile:
            for i in range(noOfWindows):
                mjdMidFile.write(f'{i}\t{midMJD[i]}\n')

    ##################################################################################
    # Read in mid-points and locations, then extract the radii/longitudes of
    # STEREO A and B that are relevant for this run
    ##################################################################################
    # Read in observation positions
    readSterAPos = bme.readObsLocFileWindows(steraLocFile, mjdCRFile, noOfWindows)
    readSterBPos = bme.readObsLocFileWindows(sterbLocFile, mjdCRFile, noOfWindows)
    readEarthPos = bme.readObsLocFileWindows(earthLocFile, mjdCRFile, noOfWindows)

    # Extract difference between Earth and observation locations
    sterAEarthPosDiff = readEarthPos - readSterAPos
    sterBEarthPosDiff = readEarthPos - readSterBPos
    aceEarthPosDiff = readEarthPos - readEarthPos

    ##########################################
    # Extract radial positions
    ##########################################
    sterARadAU = readSterAPos[1, :]
    sterBRadAU = readSterBPos[1, :]
    earthRadAU = readEarthPos[1, :]
    # At present, ACE is assumed to be at the same place as Earth. Change if necessary
    aceRadAU = readEarthPos[1, :]

    # Convert radial positions from AU to solar radii
    sterARadRs = sterARadAU * 215
    sterBRadRs = sterBRadAU * 215
    earthRadRs = earthRadAU * 215
    aceRadRs = aceRadAU * 215

    sterARadKm = sterARadRs * rS
    sterBRadKm = sterBRadRs * rS
    aceRadKm = aceRadRs * rS
    earthRadKm = earthRadRs * rS

    # Get radial coordinates of obs.
    sterARadCoord = ((sterARadKm - innerRadRs) / deltaRrs - 1).round().astype(int)
    sterBRadCoord = ((sterBRadKm - innerRadRs) / deltaRrs - 1).round().astype(int)
    aceRadCoord = ((aceRadKm - innerRadRs) / deltaRrs - 1).round().astype(int)
    earthRadCoord = ((earthRadKm - innerRadRs) / deltaRrs - 1).round().astype(int)

    # Extract longitude positions
    sterALon = np.mod(360 * np.ones(noOfWindows) - sterAEarthPosDiff[0, :], 360)
    sterBLon = np.mod(360 * np.ones(noOfWindows) - sterBEarthPosDiff[0, :], 360)
    earthLon = np.zeros(noOfWindows)
    aceLon = np.mod(360 * np.ones(noOfWindows) - aceEarthPosDiff[0, :], 360)

    # Set up model coordinates for the longitudes
    # STEREO satellite locations
    sterALonCoord = (sterALon / deltaPhiDeg).round().astype(int)
    sterBLonCoord = (sterBLon / deltaPhiDeg).round().astype(int)
    earthLonCoord = (earthLon / deltaPhiDeg).round().astype(int)
    aceLonCoord = (aceLon / deltaPhiDeg).round().astype(int)

    # Output time it has taken to get to start of windows loop
    print('\n---------------------------------------------------------------------------------')
    print(f'Time Taken to reach start of windows loop = {int((time.time() - start_time) / 60)} minutes '
          f'{np.mod(time.time() - start_time, 60)} seconds')
    print('---------------------------------------------------------------------------------\n')

    ###############################################################################
    # Read in file names of observation files
    ###############################################################################
    # Open output file to preserve printed output
    outFileName = os.path.join(outputDir, 'output.txt')
    outFile = open(outFileName, 'w')

    # Loop over all the required windows
    for w in range(noOfWindows):
        startWin_time = time.time()
        print(f'Window No.: {w}/{noOfWindows}')
        print(f'Carr. Rot.: {currCR[w]}')

        #########################################
        # Make mean and B matrix
        #########################################
        # Generate prior mean and prior covariance matrix
        unpertEnsMem, meanPrior, B = bme.createUnpertAndBMatrix(
            deltaPhiDeg, locRad, noOfLonPoints, filesVrEns[w],
            initMJDCR[w], currMJD[w], nMASens
        )

        ################################################
        # Initialise ensemble for this CR
        ################################################
        # Initialise Prior and MAS ensemble vectors
        forwardStatePrior = np.zeros((noOfRadPoints, noOfLonPoints))
        forwardStateMASMean = np.zeros((noOfRadPoints, noOfLonPoints))

        ##############################################
        # Generate prior and MAS Mean state
        ##############################################
        # MAS Mean is just mean of distribution
        forwardStateMASMean[0, :] = np.copy(meanPrior[:])
        forwardStatePrior[0, :] = unpertEnsMem.copy()

        # Reshape for use in forward model
        vPrior[w, 0, :] = forwardStatePrior[0, :].copy()
        vMASMean[w, 0, :] = forwardStateMASMean[0, :].copy()

        ##########################################################################
        # Run solar wind propagation model to get estimates of the solar wind throughout domain
        ##########################################################################
        for rIndex in range(1, noOfRadPoints):
            # Run prior state forward
            forwardStatePrior[rIndex, :] = bme.forwardRadModelNoLat(
                vPrior[w, rIndex - 1, :], vPrior[w, 0, :],
                r[rIndex - 1], rIndex - 1, deltaRrs, deltaPhi,
                solarRotFreq, alpha, rH, noOfLonPoints
            )
            vPrior[w, rIndex, :] = np.copy(forwardStatePrior[rIndex, :])

            # Run MAS Mean forward
            forwardStateMASMean[rIndex, :] = bme.forwardRadModelNoLat(
                vMASMean[w, rIndex - 1, :], vMASMean[w, 0, :],
                r[rIndex - 1], rIndex - 1, deltaRrs, deltaPhi, solarRotFreq,
                alpha, rH, noOfLonPoints
            )
            vMASMean[w, rIndex, :] = np.copy(forwardStateMASMean[rIndex, :])

        #############################################################################
        # Generate all observations
        #############################################################################
        # Initialise variables
        noOfObsA = noOfLonPoints
        noOfObsB = noOfLonPoints
        noOfObsC = noOfLonPoints

        obsToBeTakenA = range(noOfLonPoints)
        obsToBeTakenB = range(noOfLonPoints)
        obsToBeTakenC = range(noOfLonPoints)

        obsNotTakenA = []
        obsNotTakenB = []
        obsNotTakenC = []

        ###################################
        # Read observation files
        ###################################
        yA = bme.readObsFileAvg(fileObsA, fileMJD[w], currMJD[w], currMJD[w] + 27)
        yB = bme.readObsFileAvg(fileObsB, fileMJD[w], currMJD[w], currMJD[w] + 27)
        yC = bme.readObsFileAvg(fileObsACE, fileMJD[w], currMJD[w], currMJD[w] + 27)

        # Reverse order to get obs in longitude order
        yA = np.copy(yA[::-1])
        yB = np.copy(yB[::-1])
        yC = np.copy(yC[::-1])

        # Extract all obs for plotting later
        yAPlot = np.copy(yA)  # [::-1])
        yBPlot = np.copy(yB)  # [::-1])
        yCPlot = np.copy(yC)  # [::-1])

        ################################################################################
        # Filter out observations that are unrealistic/were not recorded
        # (negative SW speeds and extremely large SW speeds)
        #################################################################################
        #############
        # STEREO A
        #############
        # Initialise temporary variables
        noOfObsATemp = np.copy(noOfObsA)
        yATemp = list(np.copy(yA))
        obsToBeTakenATemp = list(np.copy(obsToBeTakenA))

        # Check if obs. speed negative or greater than 5000km/s remove observation
        for i in range(noOfObsA):
            if (abs(yA[i]) > 5000) or (yA[i] < 0):
                # Remove observation from list
                yATemp.remove(yA[i])

                # Update plotting variable as NaN
                yAPlot[i] = np.nan

                # Reduce number of obs. to be taken and record which
                # longitude obs is removed from (and which obs. are to be taken)
                noOfObsATemp = noOfObsATemp - 1
                obsToBeTakenATemp.remove(obsToBeTakenA[i])
                obsNotTakenA.append(obsToBeTakenA[i])

        ################
        # STEREO B
        ################
        # Initialise temporary variables
        noOfObsBTemp = np.copy(noOfObsB)
        yBTemp = list(np.copy(yB))
        obsToBeTakenBTemp = list(np.copy(obsToBeTakenB))

        # Check if obs. speed negative or greater than 2000km/s remove observation
        for i in range(noOfObsB):
            if (abs(yB[i]) > 5000) or (yB[i] < 0):
                # Remove observation from list
                yBTemp.remove(yB[i])

                # Update plotting variable as NaN
                yBPlot[i] = np.nan

                # Reduce number of obs. to be taken and record which
                # longitude obs is removed from (and which obs. are to be taken)
                noOfObsBTemp = noOfObsBTemp - 1
                obsToBeTakenBTemp.remove(obsToBeTakenB[i])
                obsNotTakenB.append(obsToBeTakenB[i])

        ###############
        # ACE
        ###############
        # Initialise temporary variables
        noOfObsCTemp = np.copy(noOfObsC)
        yCTemp = list(np.copy(yC[:]))
        obsToBeTakenCTemp = list(np.copy(obsToBeTakenC))

        # Check if obs. speed negative or greater than 2000km/s remove observation
        for i in range(noOfObsC):
            if (abs(yC[i]) > 5000) or (yC[i] < 0):
                # Remove observation from list
                yCTemp.remove(yC[i])

                # Update plotting variable as NaN
                yCPlot[i] = np.nan

                # Reduce number of obs. to be taken and record which
                # longitude obs is removed from (and which obs. are to be taken)
                noOfObsCTemp = noOfObsCTemp - 1
                obsToBeTakenCTemp.remove(obsToBeTakenC[i])
                obsNotTakenC.append(obsToBeTakenC[i])

        ######################################
        # Update obs. vectors
        ######################################
        # Update all observations vector
        yAllA[w * noOfLonPoints:(w + 1) * noOfLonPoints] = np.copy(yAPlot)
        yAllB[w * noOfLonPoints:(w + 1) * noOfLonPoints] = np.copy(yBPlot)
        yAllC[w * noOfLonPoints:(w + 1) * noOfLonPoints] = np.copy(yCPlot)

        # Update no. of obs. at each satellite location
        noOfObsA = np.copy(noOfObsATemp)
        noOfObsB = np.copy(noOfObsBTemp)
        noOfObsC = np.copy(noOfObsCTemp)

        # Update number of obs. to be taken
        obsToBeTakenA = np.copy(obsToBeTakenATemp)
        obsToBeTakenB = np.copy(obsToBeTakenBTemp)
        obsToBeTakenC = np.copy(obsToBeTakenCTemp)

        # Update the observation vector
        yA = np.copy(yATemp)
        yB = np.copy(yBTemp)
        yC = np.copy(yCTemp)

        #########################################################
        # Depending on what observations the user has specified to be assimilated
        # Update variables according to the radius at which the STEREO A/B satellites are
        # and the number of observations that occur at each radius
        #########################################################
        if obsToUse == 'A':
            radObs = np.zeros(1)
            nRadObs = np.zeros((len(radObs)))

            # Input radius of observations
            radObs[0] = sterARadCoord[w]

            # Input number of observations at each radius
            nRadObs[0] = noOfObsA

        elif obsToUse == 'B':
            radObs = np.zeros(1)
            nRadObs = np.zeros((len(radObs)))

            # Input radius of observations
            radObs[0] = sterBRadCoord[w]

            # Input number of observations at each radius
            nRadObs[0] = noOfObsB

        elif obsToUse == 'C':
            radObs = np.zeros(1)
            nRadObs = np.zeros((len(radObs)))

            # Input radius of observations
            radObs[0] = earthRadCoord[w]

            # Input number of observations at each radius
            nRadObs[0] = noOfObsC

        elif obsToUse == 'AB':
            radObs = np.zeros(2)
            nRadObs = np.zeros((len(radObs)))

            # Input radius of observations
            radObs[0] = sterARadCoord[w]
            radObs[1] = sterBRadCoord[w]

            # Input number of observations at each radius
            nRadObs[0] = noOfObsB
            nRadObs[1] = noOfObsA + noOfObsB

        elif obsToUse == 'AC':
            radObs = np.zeros(2)
            nRadObs = np.zeros((len(radObs)))

            # Input radius of observations
            radObs[0] = sterARadCoord[w]
            radObs[1] = earthRadCoord[w]

            # Input number of observations at each radius
            nRadObs[0] = noOfObsA
            nRadObs[1] = noOfObsA + noOfObsC
        elif obsToUse == 'BC':
            radObs = np.zeros(2)
            nRadObs = np.zeros((len(radObs)))

            # Input radius of observations
            radObs[0] = sterBRadCoord[w]
            radObs[1] = earthRadCoord[w]

            # Input number of observations at each radius
            nRadObs[0] = noOfObsB
            nRadObs[1] = noOfObsB + noOfObsC
        elif obsToUse == 'ABC':
            radObs = np.zeros(3)
            nRadObs = np.zeros((len(radObs)))

            # Input radius of observations
            radObs[0] = sterARadCoord[w]
            radObs[1] = sterBRadCoord[w]
            radObs[2] = earthRadCoord[w]

            # Input number of observations at each radius
            nRadObs[0] = noOfObsA
            nRadObs[1] = noOfObsA + noOfObsB
            nRadObs[2] = noOfObsA + noOfObsB + noOfObsC
        else:
            print(
                'obsToUse must equal "A", "B", "C", "AB", "AC", "BC" or "ABC" '
                '(or some permutation of these) without spaces'
            )
            print(
                'where A corresponds to STEREO A, B corresponds to STEREO B and'
                'C corresponds to ACE data being assimilated'
            )
            print('Please update [33] accordingly. System will now exit...')
            sys.exit()

        noOfObsTotal = int(nRadObs[-1])

        # Print number of observations
        print(f'len(yA) = {len(yA)}')
        print(f'len(yB) = {len(yB)}')
        print(f'len(yC) = {len(yC)}')
        print(f'Total number of observations: {noOfObsTotal}')

        #####################################
        # Initialise observation variables
        #####################################
        # Initialise obs.
        y = np.zeros(noOfObsTotal)

        # Obs. error covar. matrices for STERA, STERB and ACE
        RA = np.zeros((noOfObsA, noOfObsA))
        RB = np.zeros((noOfObsB, noOfObsB))
        RC = np.zeros((noOfObsC, noOfObsC))

        R = np.zeros((noOfObsTotal, noOfObsTotal))  # Obs. error covar. for all obs. to be assimilated

        ###################################################
        # Generate observation operators
        ###################################################
        H = np.zeros((noOfObsTotal, noOfLonPoints))
        HA = bme.obsOp(noOfLonPoints, obsToBeTakenA, [sterALonCoord[w]])
        HB = bme.obsOp(noOfLonPoints, obsToBeTakenB, [sterBLonCoord[w]])
        HC = bme.obsOp(noOfLonPoints, obsToBeTakenC, [aceLonCoord[w]])

        ######################################################################
        # Make observation error covariance matrices (assumed diagonal,
        # i.e. all observations are uncorrelated)
        ######################################################################
        # Read in whether a constant obs. covariance matrix is required or if
        # the uncertainty should be a percentage of the prior solar wind speed at the observation radius
        splitLine = setupOfR.split(' ')
        print(splitLine)

        # Check that the correct number of variables have been specified
        if len(splitLine) != 4:
            print(
                'Number of arguments in setupOfR should be four,'
                'each separated with spaces, with structure of:'
            )
            print('[B or C] <float> <float> <float>.')
            print('System will now exit')
            sys.exit()

        if str(splitLine[0]) == 'B':
            # Generate observation errors proportional to mean of solar wind speed at obs. radius
            obsUncA = (float(splitLine[1]) * vPrior[w, sterARadCoord[w], :].mean()) * np.ones(noOfObsA)
            obsUncB = (float(splitLine[2]) * vPrior[w, sterBRadCoord[w], :].mean()) * np.ones(noOfObsB)
            obsUncC = (float(splitLine[3]) * vPrior[w, earthRadCoord[w], :].mean()) * np.ones(noOfObsC)
        elif str(splitLine[0]) == 'C':
            # Generate observation errors as a constant supplied by the user
            obsUncA = float(splitLine[1]) * np.ones(noOfObsA)
            obsUncB = float(splitLine[2]) * np.ones(noOfObsB)
            obsUncC = float(splitLine[3]) * np.ones(noOfObsC)
        else:
            print('First character should be a "B" or "C"')
            print(('where B corresponds to a observation error standard deviation '
                   'proportional to mean prior solar wind at obs. radius'))
            print('and C corresponds to a constant observation error standard deviation being used.')
            print('Please update setupOfR accordingly. System will now exit...')
            sys.exit()

        # Assume observations at different satellites are not correlated
        for i in range(noOfObsA):
            RA[i, i] = obsUncA[i] * obsUncA[i]
        for i in range(noOfObsB):
            RB[i, i] = obsUncB[i] * obsUncB[i]
        for i in range(noOfObsC):
            RC[i, i] = obsUncC[i] * obsUncC[i]

        #########################################################################################
        # Update the generic DA variables depending upon what obs. the user wishes to assimilate
        #########################################################################################
        # Input appropriate parts into y, R and H to be used in DA
        if obsToUse == 'A':
            # Make Observations
            y[:noOfObsA] = np.copy(yA)

            # Input appropriate components into full observation covariance matrix
            R[:noOfObsA, :noOfObsA] = np.copy(RA)

            # Input appropriate components into full observation covariance matrix
            H[:noOfObsA, :] = np.copy(HA)

        elif obsToUse == 'B':
            # Make Observations
            y[:noOfObsB] = np.copy(yB)

            # Input appropriate components into full observation covariance matrix
            R[:noOfObsB, :noOfObsB] = np.copy(RB)

            # Input appropriate components into full observation covariance matrix
            H[:noOfObsB, :] = np.copy(HB)
        elif obsToUse == 'C':
            # Make Observations
            y[:noOfObsC] = np.copy(yC)

            # Input appropriate components into full observation covariance matrix
            R[:noOfObsC, :noOfObsC] = np.copy(RC)

            # Input appropriate components into full observation covariance matrix
            H[:noOfObsC, :] = np.copy(HC)

        elif obsToUse == 'AB':
            # Make Observations
            y[:noOfObsA] = np.copy(yA)
            y[noOfObsA:(noOfObsA + noOfObsB)] = np.copy(yB)

            # Input appropriate components into full observation covariance matrix
            R[:noOfObsA, :noOfObsA] = np.copy(RA)
            R[noOfObsA:(noOfObsA + noOfObsB), noOfObsA:(noOfObsA + noOfObsB)] = np.copy(RB)

            # Input appropriate components into full observation covariance matrix
            H[:noOfObsA, :] = np.copy(HA)
            H[noOfObsA:(noOfObsA + noOfObsB), :] = np.copy(HB)
        elif obsToUse == 'AC':
            # Make Observations
            y[:noOfObsA] = np.copy(yA)
            y[noOfObsA:] = np.copy(yC)

            # Input appropriate components into full observation covariance matrix
            R[:noOfObsA, :noOfObsA] = np.copy(RA)
            R[noOfObsA:, noOfObsA:] = np.copy(RC)

            # Input appropriate components into full observation covariance matrix
            H[:noOfObsA, :] = np.copy(HA)
            H[noOfObsA:, :] = np.copy(HC)
        elif obsToUse == 'BC':
            # Make Observations
            y[:noOfObsB] = np.copy(yB)
            y[noOfObsB:] = np.copy(yC)

            # Input appropriate components into full observation covariance matrix
            R[:noOfObsB, :noOfObsB] = np.copy(RB)
            R[noOfObsB:, noOfObsB:] = np.copy(RC)

            # Input appropriate components into full observation covariance matrix
            H[:noOfObsB, :] = np.copy(HB)
            H[noOfObsB:, :] = np.copy(HC)
        elif obsToUse == 'ABC':
            # Make Observations
            y[:noOfObsA] = np.copy(yA)
            y[noOfObsA:(noOfObsA + noOfObsB)] = np.copy(yB)
            y[(noOfObsA + noOfObsB):] = np.copy(yC)

            # Input appropriate components into full observation covariance matrix
            R[:noOfObsA, :noOfObsA] = np.copy(RA)
            R[noOfObsA:(noOfObsA + noOfObsB), noOfObsA:(noOfObsA + noOfObsB)] = np.copy(RB)
            R[(noOfObsA + noOfObsB):, (noOfObsA + noOfObsB):] = np.copy(RC)

            # Input appropriate components into full observation covariance matrix
            H[:noOfObsA, :] = np.copy(HA)
            H[noOfObsA:(noOfObsA + noOfObsB), :] = np.copy(HB)
            H[(noOfObsA + noOfObsB):, :] = np.copy(HC)
        else:
            print(
                'obsToUse must equal "A", "B", "C", "AB", "AC", "BC" or "ABC"'
                '(or some permutation of these) without spaces,\n'
                'where A corresponds to STEREO A, B corresponds to STEREO B and C'
                'corresponds to ACE data being assimilated\n'
                'Please update obsToUse accordingly. System will now exit...'
            )
            sys.exit()

        #############################################################################
        # Data Assimilation
        #############################################################################
        # Initialise variables to hold state before and after DA
        vIter = np.zeros((2, noOfRadPoints, noOfLonPoints))
        forwardStateIter = np.zeros((2, noOfRadPoints, noOfLonPoints))
        costFuncVal = np.zeros(2)

        # Set the initial solar wind speed as equal to the prior solar wind speed
        vIter[0, 0, :] = np.copy(vPrior[w, 0, :])
        forwardStateIter[0, 0, :] = np.copy(vIter[0, 0, :])

        # Run initial solar wind speed out into the heliosphere (from 30rS -> 215rS)
        # !Is this for loop necessary (vIter[0, :,:] = vPrior[w, :, :])?
        for rIndex in range(1, noOfRadPoints):
            forwardStateIter[0, rIndex, :] = bme.forwardRadModelNoLat(
                vIter[0, rIndex - 1, :], vIter[0, 0, :],
                r[rIndex - 1], rIndex - 1, deltaRrs, deltaPhi,
                solarRotFreq, alpha, rH, noOfLonPoints
            )
            vIter[0, rIndex, :] = np.copy(forwardStateIter[0, rIndex, :])

        # Initialise state vector variables at inner radius
        xb = forwardStateIter[0, 0, :]

        # Calculate cost function at initial iteration
        costFuncVal[0] = bme.calcCostFuncNoLat(
            B, R, H, forwardStateIter[0, 0, :], xb, vIter[0, :, :], y,
            radObs, nRadObs
        )

        ###########################################################################################
        # Minimise the cost function
        ###########################################################################################
        # Perform minimisation of cost function to obtain analysis state
        resOpt = optimize.minimize(
            fun=bme.calcCostFuncForCGNoLat, x0=xb,
            args=(
                B, R, H, xb, y, radObs, nRadObs,
                r, rH, deltaRrs, deltaPhi, alpha, solarRotFreq, noOfRadPoints, noOfLonPoints
            ), method='BFGS', jac=bme.makeGradCGNoLat, options={'gtol': gTol, 'disp': True}
        )

        ###########################################################################
        # Extract the  analysis solar wind speed in both speed matrix and state vector form
        ###########################################################################
        vIter[1, 0, :] = np.copy(resOpt.x)
        forwardStateIter[1, 0, :] = np.copy(resOpt.x)

        #####################################################################
        # Run analysed DA state out into the heliosphere using the numerical model
        #####################################################################
        for rIndex in range(1, noOfRadPoints):
            forwardStateIter[1, rIndex, :] = bme.forwardRadModelNoLat(
                vIter[1, rIndex - 1, :], vIter[1, 0, :], r[rIndex - 1], rIndex - 1,
                deltaRrs, deltaPhi, solarRotFreq, alpha, rH, noOfLonPoints
            )
            vIter[1, rIndex, :] = np.copy(forwardStateIter[1, rIndex, :])

        ###################################################################################
        # Calculate cost function after DA analysis and store in cost function variable
        ###################################################################################
        costFuncVal[1] = bme.calcCostFuncNoLat(
            B, R, H, forwardStateIter[1, 0, :], xb, vIter[1, :, :], y,
            radObs, nRadObs
        )

        print('\nMinimisation complete for this window\n')

        ################################################
        # Extract variables into larger arrays
        ################################################
        # Store analysis for current window in file
        vPosterior[w, :, :] = vIter[1, :, :]
        vPriorAll[:, (w * noOfLonPoints):((w + 1) * noOfLonPoints)] = vIter[0, :, :]
        vPosteriorAll[:, (w * noOfLonPoints):((w + 1) * noOfLonPoints)] = vIter[1, :, :]

        #############################################################
        # Extract variables for plotting at each observation location
        #############################################################
        # Initialise variables
        vPlotTempA = np.zeros((3, noOfLonPoints))
        vPlotTempB = np.zeros((3, noOfLonPoints))
        vPlotTempC = np.zeros((3, noOfLonPoints))

        # Initialise plotting variables for each obs. location
        # Column 0=Obs. values, column 1=MAS Mean speed,
        # column 2=Prior speed, column 3= Posterior speed
        vSterAPlot = np.ones((4, noOfLonPoints)) * np.nan
        vSterBPlot = np.ones((4, noOfLonPoints)) * np.nan
        vACEPlot = np.ones((4, noOfLonPoints)) * np.nan

        ######################################################
        # Extract variables for ACE satellite and order variables
        # as ascending in time (-Carr. Rot.)
        ######################################################
        # Extract velocities at STERA radii (and reorder)
        vPlotTempA[0, :] = np.copy(vMASMean[w, sterARadCoord[w], ::-1])
        vPlotTempA[1, :] = np.copy(vIter[0, sterARadCoord[w], ::-1])
        vPlotTempA[2, :] = np.copy(vIter[-1, sterARadCoord[w], ::-1])

        # Extract velocities at STERA radii (and reorder)
        vPlotTempB[0, :] = np.copy(vMASMean[w, sterBRadCoord[w], ::-1])
        vPlotTempB[1, :] = np.copy(vIter[0, sterBRadCoord[w], ::-1])
        vPlotTempB[2, :] = np.copy(vIter[-1, sterBRadCoord[w], ::-1])

        # Extract velocities at ACE radius (and reorder)
        vPlotTempC[0, :] = np.copy(vMASMean[w, aceRadCoord[w], ::-1])
        vPlotTempC[1, :] = np.copy(vIter[0, aceRadCoord[w], ::-1])
        vPlotTempC[2, :] = np.copy(vIter[-1, aceRadCoord[w], ::-1])

        ######################################################
        # Copy observation data (and reorder)
        ######################################################
        vSterAPlot[0, :] = np.copy(yAPlot[::-1])
        vSterBPlot[0, :] = np.copy(yBPlot[::-1])
        vACEPlot[0, :] = np.copy(yCPlot[::-1])

        for i in range(noOfLonPoints):
            # Extract relevant solar wind velocities at STEREO A location
            sterATrans = int(np.mod(i + sterALonCoord[w], noOfLonPoints))
            vSterAPlot[1, sterATrans] = np.copy(vPlotTempA[0, i])
            vSterAPlot[2, sterATrans] = np.copy(vPlotTempA[1, i])
            vSterAPlot[3, sterATrans] = np.copy(vPlotTempA[2, i])

            # Extract relevant solar wind velocities at STEREO B location
            sterBTrans = int(np.mod(i + sterBLonCoord[w], noOfLonPoints))
            vSterBPlot[1, sterBTrans] = np.copy(vPlotTempB[0, i])
            vSterBPlot[2, sterBTrans] = np.copy(vPlotTempB[1, i])
            vSterBPlot[3, sterBTrans] = np.copy(vPlotTempB[2, i])

            # Extract relevant solar wind velocities at ACE location
            aceTrans = int(np.mod(i + aceLonCoord[w], noOfLonPoints))
            vACEPlot[1, aceTrans] = np.copy(vPlotTempC[0, i])
            vACEPlot[2, aceTrans] = np.copy(vPlotTempC[1, i])
            vACEPlot[3, aceTrans] = np.copy(vPlotTempC[2, i])

        ########################################################################
        # Calculate RMSEs of Prior, MASMean and Posterior compared to obs. data
        ########################################################################
        # Calculate RMSEs
        RMSESTERAMASMean[w] = bme.calcStateObsRMSE(vSterAPlot[1, :], vSterAPlot[0, :])
        RMSESTERAPrior[w] = bme.calcStateObsRMSE(vSterAPlot[2, :], vSterAPlot[0, :])
        RMSESTERAPosterior[w] = bme.calcStateObsRMSE(vSterAPlot[3, :], vSterAPlot[0, :])

        # STEREO B
        RMSESTERBMASMean[w] = bme.calcStateObsRMSE(vSterBPlot[1, :], vSterBPlot[0, :])
        RMSESTERBPrior[w] = bme.calcStateObsRMSE(vSterBPlot[2, :], vSterBPlot[0, :])
        RMSESTERBPosterior[w] = bme.calcStateObsRMSE(vSterBPlot[3, :], vSterBPlot[0, :])

        # ACE
        RMSEACEMASMean[w] = bme.calcStateObsRMSE(vACEPlot[1, :], vACEPlot[0, :])
        RMSEACEPrior[w] = bme.calcStateObsRMSE(vACEPlot[2, :], vACEPlot[0, :])
        RMSEACEPosterior[w] = bme.calcStateObsRMSE(vACEPlot[3, :], vACEPlot[0, :])

        ##############################################################################
        # Plot solar wind speeds at STEREO A, STEREO B and ACE
        ##############################################################################
        swFileDir = os.path.join(outputDir, 'Plots', 'swOver27d')
        fig1, ax1 = bme.plotSWspeed(
            deltaPhiDeg, vSterAPlot, swFileDir, 'STA', currMJD[w], fontSize=18, lWid=2.0
        )
        # plt.show()

        fig2, ax2 = bme.plotSWspeed(
            deltaPhiDeg, vSterBPlot, swFileDir, 'STB', currMJD[w], fontSize=18, lWid=2.0
        )
        # plt.show()

        fig3, ax3 = bme.plotSWspeed(
            deltaPhiDeg, vACEPlot, swFileDir, 'ACE', currMJD[w], fontSize=18, lWid=2.0
        )
        # plt.show()

        ###########################################################################
        # Output RMSEs for user
        ###########################################################################
        print('------------------------------------------------------')
        print('Window no. ' + str(w))
        print(' ')
        print('RMSE STEREO A MASMean = ' + str(RMSESTERAMASMean[w]))
        print('RMSE STEREO A Prior = ' + str(RMSESTERAPrior[w]))
        print('RMSE STEREO A Posterior= ' + str(RMSESTERAPosterior[w]))
        print(' ')
        print('RMSE STEREO B MASMean = ' + str(RMSESTERBMASMean[w]))
        print('RMSE STEREO B Prior = ' + str(RMSESTERBPrior[w]))
        print('RMSE STEREO B Posterior= ' + str(RMSESTERBPosterior[w]))
        print(' ')
        print('RMSE ACE MASMean = ' + str(RMSEACEMASMean[w]))
        print('RMSE ACE Prior = ' + str(RMSEACEPrior[w]))
        print('RMSE ACE Posterior= ' + str(RMSEACEPosterior[w]))

        ###########################################################################
        # Write RMSEs into output file for further use
        ###########################################################################
        outFile.write('------------------------------------------------------\n')
        outFile.write(f'Window no. {w} \n\n')
        outFile.write(f'RMSE STEREO A MASMean = {RMSESTERAMASMean[w]} \n')
        outFile.write(f'RMSE STEREO A Prior = {RMSESTERAPrior[w]} \n')
        outFile.write(f'RMSE STEREO A Posterior = {RMSESTERAPosterior[w]} \n\n')

        outFile.write(f'RMSE STEREO B MASMean = {RMSESTERBMASMean[w]}\n')
        outFile.write(f'RMSE STEREO B Prior = {RMSESTERBPrior[w]}\n')
        outFile.write(f'RMSE STEREO B Posterior= {RMSESTERBPosterior[w]}\n\n')

        outFile.write(f'RMSE ACE MASMean = {RMSEACEMASMean[w]} \n')
        outFile.write(f'RMSE ACE Prior = {RMSEACEPrior[w]}\n')
        outFile.write(f'RMSE ACE Posterior= {RMSEACEPosterior[w]}\n')
        outFile.write('------------------------------------------------------\n\n')

        #############################################################################
        # Write MASMean, prior and posterior arrays to file for later use
        # Data output is ASCENDING IN TIME
        #############################################################################
        # MAS Mean
        outWindowMASMeanFile = os.path.join(
            outputDir, 'meanMAS', f'meanMAS_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outWindowMASMeanFile, 'w') as fOutMAS:
            np.savetxt(fOutMAS, vMASMean[w, :, ::-1])

        # Prior
        outWindowPriorFile = os.path.join(
            outputDir, 'prior', f'prior_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outWindowPriorFile, 'w') as fOutPrior:
            np.savetxt(fOutPrior, vPrior[w, :, ::-1])

        # Posterior
        outWindowPostFile = os.path.join(
            outputDir, 'posterior', f'posterior_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outWindowPostFile, 'w') as fOutPost:
            np.savetxt(fOutPost, vPosterior[w, :, ::-1])

        ######################################################################################
        # Write ACE, STEREOA and STEREOB obs and prior, post, MAS ens in obs loc.
        # during solar rotation to file for later use
        # Data output is ASCENDING IN TIME
        ######################################################################################
        ###################
        # STEREO-A
        ###################
        outSTERAObsFile = os.path.join(
            outputDir, 'STERA', f'obs_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outSTERAObsFile, 'w') as fOutSTA:
            np.savetxt(fOutSTA, yAPlot[::-1])

        outSTERAmasFile = os.path.join(
            outputDir, 'STERA', f'masMean_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outSTERAmasFile, 'w') as fOutSTA:
            np.savetxt(fOutSTA, vSterAPlot[1, :])

        outSTERApriorFile = os.path.join(
            outputDir, 'STERA', f'prior_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outSTERApriorFile, 'w') as fOutSTA:
            np.savetxt(fOutSTA, vSterAPlot[2, :])

        outSTERApostFile = os.path.join(
            outputDir, 'STERA', f'post_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outSTERApostFile, 'w') as fOutSTA:
            np.savetxt(fOutSTA, vSterAPlot[3, :])

        ###############################
        # STEREO-B
        ###############################
        outSTERBObsFile = os.path.join(
            outputDir, 'STERB', f'obs_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outSTERBObsFile, 'w') as fOutSTB:
            np.savetxt(fOutSTB, yBPlot[::-1])

        outSTERBmasFile = os.path.join(
            outputDir, 'STERB', f'masMean_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outSTERBmasFile, 'w') as fOutSTB:
            np.savetxt(fOutSTB, vSterBPlot[1, :])

        outSTERBpriorFile = os.path.join(
            outputDir, 'STERB', f'prior_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outSTERBpriorFile, 'w') as fOutSTB:
            np.savetxt(fOutSTB, vSterBPlot[2, :])

        outSTERBpostFile = os.path.join(
            outputDir, 'STERB', f'post_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outSTERBpostFile, 'w') as fOutSTB:
            np.savetxt(fOutSTB, vSterBPlot[3, :])

        ##################################
        # ACE
        ##################################
        outACEObsFile = os.path.join(
            outputDir, 'ACE', f'obs_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outACEObsFile, 'w') as fOutACE:
            np.savetxt(fOutACE, yAPlot[::-1])

        outACEmasFile = os.path.join(
            outputDir, 'ACE', f'masMean_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outACEmasFile, 'w') as fOutACE:
            np.savetxt(fOutACE, vACEPlot[1, :])

        outACEpriorFile = os.path.join(
            outputDir, 'ACE', f'prior_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outACEpriorFile, 'w') as fOutACE:
            np.savetxt(fOutACE, vACEPlot[2, :])

        outACEpostFile = os.path.join(
            outputDir, 'ACE', f'post_MJDstart{int(currMJD[w])}.txt'
        )
        with open(outACEpostFile, 'w') as fOutACE:
            np.savetxt(fOutACE, vACEPlot[3, :])

        print('\n---------------------------------------------------------------------------------')
        print(f'------ Time Taken to run window = {int((time.time() - startWin_time) / 60)} minutes '
              f'{np.mod(time.time() - startWin_time, 60)} seconds ------')
        print('---------------------------------------------------------------------------------\n')

    # Close outFile after all windows have been processed
    outFile.close()

    print('---------------------------------------------------------------------------------')
    print(f'------ Total time Taken = {int((time.time() - start_time) / 60)} minutes '
          f'{np.mod(time.time() - start_time, 60)} seconds ------')
    print('---------------------------------------------------------------------------------\n')

    return None
