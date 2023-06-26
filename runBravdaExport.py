# Import necessary packages
import numpy as np
import pandas as pd
import time
from scipy import optimize
import os
import shutil
import sys
import matplotlib.pyplot as plt

import bravdaMethodsExport as bme
from makeMASens import helioMASens


def runBravDA(configFile, huxVarFile, outputDir, obsToAssim, setupOfR,
              initDate, noOfWindows, nMASens, locRad, gTol, makePlots,
              usecustomens = False):
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
    fileObs = os.path.join(currentDir, configLines[7].strip())
    earthLocFile = os.path.join(currentDir, configLines[9].strip())

    # Read obs. location file
    with open(fileObs, 'r') as fObs:
        obsLines = fObs.readlines()

    # Calculate the number of observation streams contained within fileObs
    nObsStreams = len(obsLines) - 1  # Subtract 1 to remove the header from nObsStreams

    # Make sure obsLines contains more than just the header
    if nObsStreams < 1:
        print(f"{fileObs} should contain more than just a header.")
        print("Add observation name, observation relative filepath, observation location relative filepath,")
        print(" observation error covariance type (B = proportional to prior mean or C = Constant)")
        print(" and observation error covariance value")
        print("In the following format....")
        print("obsName obsFile obsLocFile obsErrCovType obsErrCov")
        print("<name> <rel. obsFilePath> <rel. obsLocFilePath> <B or C> <obsErrCov>")
        print("All quantities should be separated by a <space>")
        print("Don't delete the header!")

    # Extract observation file information from obsLines
    obsFileDf = pd.DataFrame(columns=["obsFilePath", "obsLocFilePath", "obsErrCovType", "obsErrCov"])
    for i, ln in enumerate(obsLines):
        # Ignore the header
        if i > 0:
            lnSplit = ln.split()

            obsName = (lnSplit[0].strip()).upper()
            obsFilePath = os.path.join(currentDir, lnSplit[1].strip())
            obsLocFilePath = os.path.join(currentDir, lnSplit[2].strip())
            obsErrCovType = lnSplit[3].strip()
            obsErrCov = float(lnSplit[4].strip())
            obsFileDf.loc[obsName] = [obsFilePath, obsLocFilePath, obsErrCovType, obsErrCov]

    # Ensure all values in obsToAssim are upper-case
    print(f'Observations to be assimilated: {obsToAssim}')
    obsToAssim = [ota.upper() for ota in obsToAssim]

    # Check all obsToAssim are named in the obsFile.dat
    for ota in obsToAssim:
        if ota not in obsFileDf.index:
            print(f"{ota} is not in {fileObs}. Please update obsToAssim or {fileObs}")
            print(" such that the required observation is in both.")
            print("System will now exit...")
            sys.exit()

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
        iwf.write(f'obsToAssim = {obsToAssim}\n')
        iwf.write(f'setupOfR = {setupOfR}\n')
        iwf.write(f'initDate = {initDate}\n')
        iwf.write(f'noOfConsecWindows = {noOfWindows}\n')
        iwf.write(f'noOfMASens = {nMASens}\n')
        iwf.write(f'locRad = {locRad}\n')
        iwf.write(f'gTol = {gTol}\n')

    # Make all necessary directories in output directory to save files
    # List containing all folders to make
    dirsToMake = ['meanMAS', 'prior', 'posterior']
    for dirName in dirsToMake:
        dirPath = os.path.join(outputDir, dirName, '')

        # Check if dir exists, if not, make it
        if not os.path.isdir(dirPath):
            os.makedirs(dirPath)

    #Make directories to hold all observation files
    for dirName in obsFileDf.index:
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

    ########################################################
    # Initialise arrays for use later in script
    ########################################################
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
            bme.makeMJDfile(currMJD[-1], noOfLonPoints, fileMJDpath, daySolarRot= 27.2753)

    #####################################################################################
    # Check existence of and make initial ensemble files if required for each window
    #####################################################################################
    # Directory to hold solar wind ens files specific to this run
    filesVrEns = []  # Variable to hold required initial solar wind speed ensemble file names for this run
    if usecustomens:
        filesVrEns.append(f'{dirMASens}customensemble.dat')
    else:
        for w in range(noOfWindows):
            filesVrEns.append(f'{dirMASens}vin_ensemble_CR{currCR[w]}.dat')
           
    
            if not os.path.isfile(filesVrEns[w]):
                # If MAS ensemble file does not exist, make it
                print(f'Generating MAS ensembles for CR {currCR[w]}...')
    
                # Generate path to ephemerisMSL.hdf5 file (should be in makeMASens, don't move)
                ephemFile = os.path.join(currentDir,'makeMASens', 'ephemerisMSL.hdf5')
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
        outputDir, 'MJDfiles', f'MJDMidPoints_{int(currMJD[0])}_{int(currMJD[-1])}.dat'
    )
    if not os.path.isfile(mjdCRFile):
        with open(mjdCRFile, 'w') as mjdMidFile:
            for i in range(noOfWindows):
                mjdMidFile.write(f'{i}\t{midMJD[i]}\n')

    ##################################################################################
    # Read in mid-points and locations, then extract the radii/longitudes of
    # STEREO A and B that are relevant for this run
    ##################################################################################
    obsPosDf = pd.DataFrame(columns=["radAU", "radRs", "radKm", "radCoord", "lon", "lonCoord"])
    for obsName in obsFileDf.index:
        fObsLoc = obsFileDf.loc[obsName]["obsLocFilePath"]
        radAU, radRs, radKm, radCoord, lon, lonCoord = bme.getObsPos(
            fObsLoc, earthLocFile, mjdCRFile, noOfWindows,
            rS, innerRadRs, deltaRrs, deltaPhiDeg
        )
        obsPosDf.loc[obsName] = [radAU, radRs, radKm, radCoord, lon, lonCoord]

    print(obsPosDf)

    # Output time it has taken to get to start of windows loop
    print('\n---------------------------------------------------------------------------------')
    print(f'Time Taken to reach start of windows loop = {int((time.time() - start_time) / 60)} minutes '
          f'{np.mod(time.time() - start_time, 60)} seconds')
    print('---------------------------------------------------------------------------------\n')

    ###############################################################################
    # Open output file to preserve printed output
    #################################################
    outFileName = os.path.join(outputDir, 'output.txt')
    outFile = open(outFileName, 'w')

    # Loop over all the required windows
    for w in range(noOfWindows):
        startWin_time = time.time()
        print(f'Window No.: {w}/{noOfWindows}')
        print(f'Carr. Rot.: {currCR[w]}')

        B, forwardStatePrior, forwardStateMASMean = bme.makePriors(
            filesVrEns[w], initMJDCR[w], currMJD[w], locRad, nMASens,
            r, rH, deltaRrs, deltaPhi, alpha, solarRotFreq,
            noOfRadPoints, noOfLonPoints, incMASMean=True
        )

        vPrior[w, :, :] = np.copy(forwardStatePrior)
        vMASMean[w, :, :] = np.copy(forwardStateMASMean)

        #############################################################################
        # Generate all observations
        #############################################################################
        obsCompDf = pd.DataFrame(columns=["y", "yPlot", "obsToBeTaken", "noOfObs"])

        for obsName in obsFileDf.index:
            fObs = obsFileDf.loc[obsName]["obsFilePath"]
            yTemp, yTempPlot, obsToBeTakenTemp, noOfObsTemp = bme.makeObs(
                fObs, fileMJD[w], currMJD[w], noOfLonPoints,
                lowerCutOff=0, upperCutOff=5000
            )
            obsCompDf.loc[obsName] = [yTemp, yTempPlot, obsToBeTakenTemp, noOfObsTemp]

        # Place all obs into a dictionary
        obsCompDict = obsCompDf.to_dict()["y"]

        # Set up dictionaries containing the number of obs and which are to be taken
        nObsDict = obsCompDf.to_dict()["noOfObs"]
        obsToBeTakenDict = obsCompDf.to_dict()["obsToBeTaken"]

        # Set up dictionaries containing the obsservation's positional data
        radCoordDict = {}
        lonCoordDict = {}
        for index, row in obsPosDf.iterrows():
            radCoordDict[index] = row["radCoord"][w]
            lonCoordDict[index] = row["lonCoord"][w]

        # Extract the observation error uncertainty
        obsUncDf = obsFileDf[["obsErrCovType", "obsErrCov"]]

        #####################################
        # Initialise observation variables
        #####################################
        # Initialise obs.
        y = []
        H = [] # Obs. operator
        R = []  # Obs. error covar. for all obs. to be assimilated

        #Initialise radObs and nRadObs
        radObs = []
        nRadObs = []

        # Extract radial positions from each obs. source
        for obsName in obsToAssim:
            radObs, nRadObs = bme.extractRadObs(obsName, radCoordDict, nObsDict, radObs, nRadObs)

            y, H, R = bme.makeObsForDA(
                y, H, R, obsCompDict, obsName, radCoordDict, lonCoordDict,
                obsUncDf, obsToBeTakenDict, nObsDict, vPrior[w, :, :], noOfLonPoints
            )

        # Extract total number of observations
        noOfObsTotal = int(nRadObs[-1])

        # Print number of observations
        for obsName in obsCompDf.index:
            print(f'len(y_{obsName}) = {len(obsCompDf.loc[obsName]["y"])}')
        print(f'Total number of observations: {noOfObsTotal}')

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
        #print(f'y={y}')
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
        vPlotTempA[0, :] = np.copy(vMASMean[w, obsPosDf.loc["A"]["radCoord"][w], ::-1])
        vPlotTempA[1, :] = np.copy(vIter[0, obsPosDf.loc["A"]["radCoord"][w], ::-1])
        vPlotTempA[2, :] = np.copy(vIter[-1, obsPosDf.loc["A"]["radCoord"][w], ::-1])

        # Extract velocities at STERA radii (and reorder)
        vPlotTempB[0, :] = np.copy(vMASMean[w, obsPosDf.loc["B"]["radCoord"][w], ::-1])
        vPlotTempB[1, :] = np.copy(vIter[0, obsPosDf.loc["B"]["radCoord"][w], ::-1])
        vPlotTempB[2, :] = np.copy(vIter[-1, obsPosDf.loc["B"]["radCoord"][w], ::-1])

        # Extract velocities at ACE radius (and reorder)
        vPlotTempC[0, :] = np.copy(vMASMean[w, obsPosDf.loc["C"]["radCoord"][w], ::-1])
        vPlotTempC[1, :] = np.copy(vIter[0, obsPosDf.loc["C"]["radCoord"][w], ::-1])
        vPlotTempC[2, :] = np.copy(vIter[-1, obsPosDf.loc["C"]["radCoord"][w], ::-1])

        ######################################################
        # Copy observation data (and reorder)
        ######################################################
        vSterAPlot[0, :] = np.copy(yAPlot[::-1])
        vSterBPlot[0, :] = np.copy(yBPlot[::-1])
        vACEPlot[0, :] = np.copy(yCPlot[::-1])

        for i in range(noOfLonPoints):
            # Extract relevant solar wind velocities at STEREO A location
            sterATrans = int(np.mod(i + obsPosDf.loc["A"]["lonCoord"][w], noOfLonPoints))
            vSterAPlot[1, sterATrans] = np.copy(vPlotTempA[0, i])
            vSterAPlot[2, sterATrans] = np.copy(vPlotTempA[1, i])
            vSterAPlot[3, sterATrans] = np.copy(vPlotTempA[2, i])

            # Extract relevant solar wind velocities at STEREO B location
            sterBTrans = int(np.mod(i + obsPosDf.loc["B"]["lonCoord"][w], noOfLonPoints))
            vSterBPlot[1, sterBTrans] = np.copy(vPlotTempB[0, i])
            vSterBPlot[2, sterBTrans] = np.copy(vPlotTempB[1, i])
            vSterBPlot[3, sterBTrans] = np.copy(vPlotTempB[2, i])

            # Extract relevant solar wind velocities at ACE location
            aceTrans = int(np.mod(i + obsPosDf.loc["C"]["lonCoord"][w], noOfLonPoints))
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
        if makePlots:
            swFileDir = os.path.join(outputDir, 'plots', 'swOver27d')
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
            np.savetxt(fOutSTA, vSterAPlot[0, :])

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
            np.savetxt(fOutSTB, vSterBPlot[0, :])

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
            np.savetxt(fOutACE, vACEPlot[0, :])

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
