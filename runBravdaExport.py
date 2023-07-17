# Import necessary packages
import numpy as np
import pandas as pd
import time
import scipy
import os
import shutil
import sys
import matplotlib.pyplot as plt

import bravdaMethodsExport as bme
from makeMASens import helioMASens
import lineSearchMin as lsm


def runBravDA(configFile, huxVarFile, outputDir, obsToAssim, setupOfR,
              initDate, noOfWindows, nMASens, locRad, gTol, makePlots,
              usecustomens = False, useLogTrans=False, precondState=False):
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

    # Make directories to hold all observation files
    for dirName in obsFileDf.index:
        dirPath = os.path.join(outputDir, dirName, '')
        print(dirPath)
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

    if useLogTrans:
        vLogPrior = np.zeros((noOfWindows, noOfRadPoints, noOfLonPoints))

    # If user wishes to plot all windows consecutively, this array holds the solar wind speeds sequentially
    vPriorAll = np.zeros((noOfRadPoints, totalLonPoints))
    vPosteriorAll = np.zeros((noOfRadPoints, totalLonPoints))

    # Create dataframe to hold all RMSE values over all windows
    rmseAllWinDf = pd.DataFrame(columns=["ensMean", "prior", "posterior"])

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

        if precondState:
            if useLogTrans:
                Bhalf, forwardStatePrior, forwardStateMASMean, BlogHalf, forwardLogPrior = bme.makePriors(
                    filesVrEns[w], initMJDCR[w], currMJD[w], locRad, nMASens, r, rH, deltaRrs, deltaPhi, alpha,
                    solarRotFreq, noOfRadPoints, noOfLonPoints, useBlogTrans=useLogTrans, precondState=True
                )

                vLogPrior[w, :, :] = np.copy(forwardLogPrior)
            else:
                Bhalf, forwardStatePrior, forwardStateMASMean = bme.makePriors(
                    filesVrEns[w], initMJDCR[w], currMJD[w], locRad, nMASens, r, rH, deltaRrs, deltaPhi, alpha,
                    solarRotFreq, noOfRadPoints, noOfLonPoints, useBlogTrans=useLogTrans, precondState=True
                )
        else:
            B, forwardStatePrior, forwardStateMASMean = bme.makePriors(
                filesVrEns[w], initMJDCR[w], currMJD[w], locRad, nMASens, r, rH, deltaRrs, deltaPhi, alpha,
                solarRotFreq, noOfRadPoints, noOfLonPoints, precondState=False
            )

        vPrior[w, :, :] = np.copy(forwardStatePrior)
        vMASMean[w, :, :] = np.copy(forwardStateMASMean)

        #############################################################################
        # Generate all observations
        #############################################################################
        # if useLogTrans:
        #     obsCompDf = pd.DataFrame(columns=["y", "yLog", "yPlot", "obsToBeTaken", "noOfObs"])
        # else:
        obsCompDf = pd.DataFrame(columns=["y", "yPlot", "obsToBeTaken", "noOfObs"])

        for obsName in obsFileDf.index:
            # if useLogTrans:
            #     fObs = obsFileDf.loc[obsName]["obsFilePath"]
            #     yTemp, yLogTemp, yTempPlot, obsToBeTakenTemp, noOfObsTemp = bme.makeObs(
            #         fObs, fileMJD[w], currMJD[w], noOfLonPoints,
            #         lowerCutOff=0, upperCutOff=5000, useLogTrans=useLogTrans
            #     )
            #     obsCompDf.loc[obsName] = [yTemp, yLogTemp, yTempPlot, obsToBeTakenTemp, noOfObsTemp]
            # else:
            fObs = obsFileDf.loc[obsName]["obsFilePath"]
            yTemp, yTempPlot, obsToBeTakenTemp, noOfObsTemp = bme.makeObs(
                fObs, fileMJD[w], currMJD[w], noOfLonPoints,
                lowerCutOff=0, upperCutOff=5000
            )
            obsCompDf.loc[obsName] = [yTemp, yTempPlot, obsToBeTakenTemp, noOfObsTemp]

        # Place all obs into a dictionary
        obsCompDict = obsCompDf.to_dict()["y"]

        # if useLogTrans:
        #     obsLogCompDict = obsCompDf.to_dict()["yLog"]

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

        # if useLogTrans:
        #     yLog = []
        #     Rlog = []

        #Initialise radObs and nRadObs
        radObs = []
        nRadObs = []

        # Extract radial positions from each obs. source
        for obsName in obsToAssim:
            radObs, nRadObs = bme.extractRadObs(obsName, radCoordDict, nObsDict, radObs, nRadObs)

            # if useLogTrans:
            #     y, H, R, yLog, Rlog = bme.makeObsForDA(
            #         y, H, R, obsCompDict, obsName, radCoordDict, lonCoordDict,
            #         obsUncDf, obsToBeTakenDict, nObsDict, vPrior[w, :, :], noOfLonPoints,
            #         useLogTrans=useLogTrans, yLogDict=obsLogCompDict, yLog=yLog, Rlog=Rlog
            #     )
            # else:
            y, H, R = bme.makeObsForDA(
                y, H, R, obsCompDict, obsName, radCoordDict, lonCoordDict,
                obsUncDf, obsToBeTakenDict, nObsDict, vPrior[w, :, :], noOfLonPoints
            )

        # Extract total number of observations
        noOfObsTotal = int(nRadObs[-1])
        lenRadObs = len(radObs)

        # Print number of observations
        for obsName in obsCompDf.index:
            print(f'len(y_{obsName}) = {len(obsCompDf.loc[obsName]["y"])}')
        print(f'Total number of observations: {noOfObsTotal}')
        preCondState2 = True
        #############################################################################
        # Data Assimilation
        #############################################################################
        if preCondState2:
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
                np.transpose(Bhalf).dot(Bhalf), R, H, forwardStateIter[0, 0, :], xb, vIter[0, :, :], y,
                radObs, nRadObs
            )

            # Precondition the state
            chi = np.zeros((noOfLonPoints))

            ###########################################################################################
            # Minimise the cost function
            ###########################################################################################
            # Perform minimisation of cost function to obtain analysis state
            resOpt = scipy.optimize.minimize(
                fun=bme.calcCostFuncPrecond2, x0=chi,
                args=(
                    xb, Bhalf, R, H, y, radObs, nRadObs,
                    r, rH, deltaRrs, deltaPhi, alpha, solarRotFreq,
                    noOfRadPoints, noOfLonPoints, False, False
                ), method='BFGS', jac=bme.makeGradCGPrecond2, options={'gtol': gTol, 'disp': True}
            )

            ###########################################################################
            # Extract the  analysis solar wind speed in both speed matrix and state vector form
            ###########################################################################
            vIter[1, 0, :] = xb + Bhalf.dot(resOpt.x)
            forwardStateIter[1, 0, :] = xb + Bhalf.dot(resOpt.x)

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
                np.transpose(Bhalf).dot(Bhalf), R, H, forwardStateIter[1, 0, :], xb, vIter[1, :, :], y,
                radObs, nRadObs
            )
        elif precondState:
            print("In precondState loop")
            outerLoopIter = 5
            innerLoopIter = 1

            # Initialise variables to hold state before and after DA
            vIter = np.zeros((2, noOfRadPoints, noOfLonPoints))
            vOuterIter = np.zeros((outerLoopIter + 1, noOfRadPoints, noOfLonPoints))
            totalVarOuter = np.zeros(outerLoopIter)
            innov = np.zeros((outerLoopIter, noOfObsTotal))
            chiIter = np.zeros((outerLoopIter + 1, noOfLonPoints))
            chiB = np.zeros((outerLoopIter + 1, noOfLonPoints))
            costFuncInnerVal = np.zeros(innerLoopIter + 1)
            costFuncOuterVal = np.zeros(outerLoopIter + 1)

            # Set the initial solar wind speed as equal to the prior solar wind speed
            vIter[0, :, :] = np.copy(vPrior[w, :, :])
            vOuterIter[0, 0, :] = np.copy(vPrior[w, 0, :])

            if useLogTrans:
                vLogIter = np.zeros((outerLoopIter + 1, noOfRadPoints, noOfLonPoints))
                vLogIter[0, 0, :] = np.copy(np.log(vOuterIter[0, 0, :] - 200))
            ##########################TRANSFORM TO LOG SPACE!!!!!!!!!
            # vPriorLogSpace = np.log(vPrior[w, :, :] - 200)
            ###########################################################

            for outer in range(outerLoopIter):
                printCost = False
                ###################################################################
                # Run outer loop
                ###################################################################
                # vOuterIter is the state we linearise around
                if outer == 0:
                    for rIndex in range(1, noOfRadPoints):
                        vOuterIter[outer, rIndex, :] = bme.forwardRadModelNoLat(
                            vOuterIter[outer, rIndex - 1, :], vOuterIter[outer, 0, :], r[rIndex - 1], rIndex - 1,
                            deltaRrs, deltaPhi, solarRotFreq, alpha, rH, noOfLonPoints
                        )

                ################################################################
                # Change vOuterIter from velocity space to ln(v - 200) space
                ################################################################


                # Initialise chi
                chi = np.zeros((innerLoopIter + 1, noOfLonPoints))

                #vB = vOuterIter[outer, 0, :] + Bhalf.dot(chiB[outer, :])
                if useLogTrans:
                    vLogIter[outer, 0, :] = np.log(vOuterIter[outer, 0, :] - 200)
                    print(f"vLogPrior = {vLogPrior[w, 0, :].max()}")
                    vPriCalc = vLogIter[outer, 0, :] + BlogHalf.dot(chiB[outer, :])
                    print(f"vOuterIter[{outer}] + Bhalf*(chiB) = {np.abs(vPriCalc).max()}")

                    print(f"vPrior = {np.abs(vPrior[w, 0, :]).max()}")
                    vPriCalc = 200 + np.exp(vLogIter[outer, 0, :] + BlogHalf.dot(chiB[outer, :]))
                        # 200 + np.multiply(vOuterIter[outer, 0, :] - 200, np.exp(BlogHalf.dot(chiB[outer, :])))
                    print(
                        f"vOuterIter[{outer}] + Bhalf*(chiB) = {np.abs(vPriCalc).max()}")
                else:
                    print(f"vPrior = {np.abs(vPrior[w, 0, :]).max()}")
                    vPriCalc = vOuterIter[outer, 0, :] + Bhalf.dot(chiB[outer, :])
                    print(f"vOuterIter[{outer}] + Bhalf*(chiB) = {np.abs(vPriCalc).max()}")
                # plt.imshow(BlogHalf)
                # plt.show()

                # Calculate non-linear cost function
                #print("Initial non-linear cost function")
                if outer == 0:
                    costFuncOuterVal[outer] = bme.calcCostFuncNoLat(
                        np.transpose(Bhalf).dot(Bhalf), R, H, vOuterIter[outer, 0, :],
                        vPrior[w, 0, :], vOuterIter[outer, :, :], y, radObs, nRadObs
                    )

                    # bme.calcNonlinCostFuncPrecond(
                    #     chi[0, :], chiB[outer, :], R, H, vOuterIter[outer, :, :],
                    #     y, radObs, nRadObs, printCF=printCost
                    # )

                # Initialise innovation vector
                # if useLogTrans:
                #     innov = bme.calcInnovLog(
                #         y, H, vLogIter[outer, :, :], radObs, nRadObs,
                #         r, deltaRrs, deltaPhi, solarRotFreq, alpha, rH, noOfRadPoints, noOfLonPoints
                #     )
                # else:
                innov[outer, :] = bme.calcInnov(y, H, vOuterIter[outer, :, :], radObs, nRadObs)
                print(obsPosDf.loc["ACE"].keys())
                # fig = plt.figure()
                # plt.plot(range(len(y)), innov[outer - 1, :], color='r')
                # plt.plot(range(len(y)), innov[outer, :])
                # plt.plot(range(len(y)), y-vOuterIter[outer, obsPosDf.loc["ACE"]["radCoord"][0], :], color='m')
                # plt.show()
                # Initialise variables for the CG algorithm
                residuals = np.zeros((innerLoopIter + 1, noOfLonPoints))
                searchDir = np.zeros((innerLoopIter + 1, noOfLonPoints))
                a_j = np.zeros(innerLoopIter)
                b_j = np.zeros(innerLoopIter)

                # Initialise variables for the Lanczos algorithm
                q_j = np.zeros((innerLoopIter, noOfLonPoints))
                alpha_j = np.zeros(innerLoopIter)
                beta_j = np.zeros(innerLoopIter)

                for inner in range(innerLoopIter):
                    print(f"Loop {outer}.{inner}")
                    # Use CG/Lanczos algorithm to update and minimise cost function and estimate
                    # the eigenvalues/eigenvectors of the Hessian
                    # if inner == 0:
                    #     # Calculate preconditioned cost function
                    #     if useLogTrans:
                    #         costFuncInnerVal[0] = bme.calcCostFuncPrecond(
                    #             chi[0, :], chiB[outer, :], Bhalf, R, H, vOuterIter[outer, :, :], innov[outer, :], radObs,
                    #             nRadObs, r, rH, deltaR, deltaPhi, alpha, solarRotFreq, noOfRadPoints, noOfLonPoints,
                    #             printCF=printCost, useLogTrans=useLogTrans
                    #         )
                    #
                    #         # Calculate gradient at chi=0
                    #         gradient = bme.makeGradCGPrecond(
                    #             chi[0, :], chiB[outer, :], Bhalf, R, H, vOuterIter[outer, :, :], innov[outer, :], radObs,
                    #             nRadObs, r, rH, deltaR, deltaPhi, alpha, solarRotFreq, noOfRadPoints, noOfLonPoints,
                    #             useLogTrans=useLogTrans
                    #         )
                    #     else:
                    #         costFuncInnerVal[0] = bme.calcCostFuncPrecond(
                    #             chi[0, :], chiB[outer, :], Bhalf, R, H, vOuterIter[outer, :, :], innov[outer, :], radObs, nRadObs,
                    #             r, rH, deltaR, deltaPhi, alpha, solarRotFreq, noOfRadPoints, noOfLonPoints,
                    #             printCF=printCost
                    #         )
                    #
                    #         # Calculate gradient at chi=0
                    #         gradient = bme.makeGradCGPrecond(
                    #             chi[0, :], chiB[outer, :], Bhalf, R, H, vOuterIter[outer, :, :], innov[outer, :], radObs, nRadObs,
                    #             r, rH, deltaR, deltaPhi, alpha, solarRotFreq, noOfRadPoints, noOfLonPoints
                    #         )

                        # Set r_0 = p_0 = -grad(J(0))
                        # Check chi is zero for first iteration
                        #assert chi[0, :].all() == 0

                    # Minimise cost function using fmin_cg
                    if useLogTrans:
                        chiIter[outer + 1, :] = scipy.optimize.minimize(
                            fun=bme.calcCostFuncPrecond, x0=chi[0, :],
                            args=(
                                chiB[outer, :], BlogHalf, R, H, vOuterIter[outer, :, :], innov[outer, :], radObs,
                                nRadObs,
                                r, rH, deltaR, deltaPhi, alpha, solarRotFreq,
                                noOfRadPoints, noOfLonPoints, False, useLogTrans
                            ), method='CG', jac=bme.makeGradCGPrecond, options={'gtol': gTol, 'disp': True}
                        ).x
                        # chiIter[outer + 1, :] = lsm.fmin_cg(
                        #     f=bme.calcCostFuncPrecond, x0=chi[0, :], fprime=bme.makeGradCGPrecond, maxiter=1,
                        #     args=(
                        #         chiB[outer, :], BlogHalf, R, H, vOuterIter[outer, :, :], innov[outer, :], radObs, nRadObs,
                        #         r, rH, deltaR, deltaPhi, alpha, solarRotFreq,
                        #         noOfRadPoints, noOfLonPoints, False, useLogTrans
                        #     )
                        # )
                    else:
                        chiIter[outer+1, :] = lsm.fmin_cg(
                            f=bme.calcCostFuncPrecond, x0=chi[0, :], fprime=bme.makeGradCGPrecond,
                            args=(
                                chiB[outer, :], Bhalf, R, H, vOuterIter[outer, :, :], innov[outer, :], radObs, nRadObs,
                                r, rH, deltaR, deltaPhi, alpha, solarRotFreq,
                                noOfRadPoints, noOfLonPoints, False
                            )
                        )

                    #fmin_cg(f, x0, fprime=gradf, args=args)
                    #     # Set first search direction and residual as -gradient(chi=0)
                    #     residuals[0, :] = -gradient
                    #     searchDir[0, :] = -gradient
                    #
                    #     print(f"pk^T * gradJ = {np.transpose(searchDir[0, :]).dot(gradient)}")
                    #     # Calculate a0 (step length for first step)
                    #     resLS = scipy.optimize.line_search(
                    #         bme.calcCostFuncPrecond, bme.makeGradCGPrecond, chi[0, :], searchDir[0, :],
                    #         gfk=gradient, old_fval=costFuncInnerVal[0], c1=1e-4, c2=0.1, maxiter=5,
                    #         args=(
                    #             chiB[outer, :], Bhalf, R, H, vOuterIter[outer, :, :], innov, radObs, nRadObs,
                    #             r, rH, deltaR, deltaPhi, alpha, solarRotFreq,
                    #             noOfRadPoints, noOfLonPoints, False
                    #         )
                    #     )
                    #     a_j[0] = resLS[0]
                    #
                    #     if np.isnan(a_j[0]):
                    #         print("Ending inner loop early due to lack of convergence")
                    #         break
                    #     # a_j[0] = bme.lineSearchFunc(
                    #     #     1, chi[0, :], chiB[outer, :], searchDir[0, :], residuals[0, :],
                    #     #     vOuterIter[outer, :, :], Bhalf, innov, H, R,
                    #     #     r, deltaRrs, deltaPhi, solarRotFreq,
                    #     #     alpha, rH, noOfRadPoints, noOfLonPoints, radObs, nRadObs,
                    #     #     maxIter=5, tau=0.5, tolerance=1e-4
                    #     # )
                    #     # resOpt = scipy.optimize.minimize(
                    #     #     fun=bme.lineSearchFunc, x0=0,
                    #     #     args=(
                    #     #         chi[0, :], chiB[outer, :], searchDir[0, :], vOuterIter[outer, :, :], Bhalf, innov, H, R,
                    #     #         r, deltaRrs, deltaPhi, solarRotFreq,
                    #     #         alpha, rH, noOfRadPoints, noOfLonPoints, radObs, nRadObs
                    #     #     ), options={'disp': True}
                    #     # )
                    #     # a_j[0] = np.copy(resOpt.x)
                    #     # Calculate new state vector
                    #     chi[1, :] = chi[0, :] + (a_j[0] * searchDir[0, :])
                    #
                    #     costFuncInnerVal[1] = bme.calcCostFuncPrecond(
                    #         chi[0, :], chiB[outer, :], Bhalf, R, H, vOuterIter[outer, :, :], innov, radObs, nRadObs,
                    #         r, rH, deltaR, deltaPhi, alpha, solarRotFreq, noOfRadPoints, noOfLonPoints,
                    #         printCF=printCost
                    #     )
                    #     # Calculate new residual (negative of gradient)
                    #     gradJLin = bme.makeGradCGPrecond(
                    #         chi[1, :], chiB[outer, :], Bhalf, R, H, vOuterIter[outer, :, :], innov, radObs, nRadObs,
                    #         r, rH, deltaR, deltaPhi, alpha, solarRotFreq, noOfRadPoints, noOfLonPoints
                    #     )
                    #     residuals[1, :] = -gradJLin
                    #
                    #     # Calculate b_0
                    #     num = np.transpose(residuals[1, :]).dot(residuals[1, :])
                    #     denom = np.transpose(residuals[0, :]).dot(residuals[0, :])
                    #     b_j[0] = num / denom
                    #     #print(f"b_j[0] = {b_j[0]}")
                    #     # Calculate new search direction
                    #     searchDir[1, :] = residuals[1, :] + (b_j[0] * searchDir[0, :])
                    #
                    #     # Calculate Lanczos algorithm quantities
                    #     magResid = np.sqrt(residuals[0, :].dot(residuals[0, :]))
                    #     q_j[0, :] = residuals[0, :] / magResid
                    #
                    #     # Calculate alpha and beta values from CG scalars
                    #     alpha_j[0] = 1 / a_j[0]
                    #     beta_j[0] = np.sqrt(b_j[0]) / a_j[0]
                    # else:
                    #     # Calculate a_j (step length for j^th step)
                    #     print(f"pk^T * gradJ = {np.transpose(searchDir[inner, :]).dot(gradJLin)}")
                    #
                    #     resLS = scipy.optimize.line_search(
                    #         bme.calcCostFuncPrecond, bme.makeGradCGPrecond, chi[inner, :], searchDir[inner, :],
                    #         gfk=gradJLin, old_fval=costFuncInnerVal[inner], c1=1e-4, c2=0.1, maxiter=5,
                    #         args=(
                    #             chiB[outer, :], Bhalf, R, H, vOuterIter[outer, :, :], innov, radObs, nRadObs,
                    #             r, rH, deltaR, deltaPhi, alpha, solarRotFreq,
                    #             noOfRadPoints, noOfLonPoints, False
                    #         )
                    #     )
                    #     a_j[inner] = resLS[0]
                    #     if np.isnan(a_j[inner]):
                    #         print("Ending inner loop early due to lack of convergence")
                    #         break
                    #     # a_j[inner] = bme.lineSearchFunc(
                    #     #     1, chi[inner, :], chiB[outer, :], searchDir[inner, :], residuals[inner, :],
                    #     #     vOuterIter[outer, :, :], Bhalf, innov, H, R,
                    #     #     r, deltaRrs, deltaPhi, solarRotFreq,
                    #     #     alpha, rH, noOfRadPoints, noOfLonPoints, radObs, nRadObs,
                    #     #     maxIter=5, tau=0.5, tolerance=1e-4
                    #     # )
                    #     # resOpt = scipy.optimize.minimize(
                    #     #     fun=bme.lineSearchFunc, x0=0,
                    #     #     args=(
                    #     #         chi[inner, :], chiB[outer, :], searchDir[inner, :], vOuterIter[outer, :, :], Bhalf, innov, H, R,
                    #     #         r, deltaRrs, deltaPhi, solarRotFreq,
                    #     #         alpha, rH, noOfRadPoints, noOfLonPoints, radObs, nRadObs
                    #     #     ), options={'disp': True}
                    #     # )
                    #     # a_j[inner] = np.copy(resOpt.x)
                    #
                    #     # Calculate new state vector
                    #     chi[inner + 1, :] = chi[inner, :] + (a_j[inner] * searchDir[inner, :])
                    #     #print(f"a_j={a_j[inner]}")
                    #     #print(f"chi={chi[inner + 1]}")
                    #     costFuncInnerVal[inner + 1] = bme.calcCostFuncPrecond(
                    #         chi[inner + 1, :], chiB[outer, :], Bhalf, R, H, vOuterIter[outer, :, :], innov, radObs, nRadObs,
                    #         r, rH, deltaR, deltaPhi, alpha, solarRotFreq, noOfRadPoints, noOfLonPoints,
                    #         printCF=printCost
                    #     )
                    #     # Calculate new residual (negative of gradient)
                    #     gradJLin = bme.makeGradCGPrecond(
                    #         chi[inner + 1, :], chiB[outer, :], Bhalf, R, H, vOuterIter[outer, :, :], innov,
                    #         radObs, nRadObs, r, rH, deltaR, deltaPhi, alpha, solarRotFreq,
                    #         noOfRadPoints, noOfLonPoints
                    #     )
                    #     residuals[inner + 1, :] = -gradJLin
                    #
                    #     # Calculate b_0
                    #     num = np.transpose(residuals[inner + 1, :]).dot(residuals[inner + 1, :])
                    #     denom = np.transpose(residuals[inner, :]).dot(residuals[inner, :])
                    #     b_j[inner] = num / denom
                    #     #print(f"b_j[{inner}] = {b_j[inner]}")
                    #
                    #     # Calculate new search direction
                    #     searchDir[inner + 1, :] = residuals[inner + 1, :] + (b_j[inner] * searchDir[inner, :])

                        # Calculate Lanczos algorithm quantities
                #         magResid = np.sqrt(residuals[inner, :].dot(residuals[inner, :]))
                #         q_j[inner, :] = (-1**inner) * residuals[inner, :] / magResid
                #
                #         # Calculate alpha and beta values from CG scalars
                #         alpha_j[inner] = (1 / a_j[inner]) + (b_j[inner - 1] / a_j[inner - 1])
                #         beta_j[inner] = np.sqrt(b_j[inner]) / a_j[inner]
                #
                #
                # # Estimate the eigen-values of Hessian matrix from the inner-loop alpha, beta and q's
                # mainDiag = np.diag(a_j, k=0)
                # lowerDiag = np.diag(b_j[:(innerLoopIter - 1)], k=-1)
                # upperDiag = np.diag(b_j[:(innerLoopIter - 1)], k=1)
                # triDiag = mainDiag + lowerDiag + upperDiag
                #
                # # Calculate eigenvalues of triDiag
                # tEval, tEvect = scipy.linalg.eigh_tridiagonal(a_j, b_j[:(innerLoopIter - 1)])
                #
                # # Transform to get e-vectors of preconditioned Hessian
                # hessEvect = np.matrix(np.transpose(q_j).dot(tEvect))
                #
                # # Estimate the preconditioned inverse Hessian to get
                # # an estimate of the preconditioned posterior error cov with
                # # A = I + \sum_{i=0}^innerLoop [(tEval^{-1} - 1)(hessEvect*hessEvect^T)
                # invTEval = [1 / t for t in tEval]
                # sqHess = [
                #     np.transpose(hessEvect[i, :]) @ hessEvect[i, :] for i in range(innerLoopIter)
                # ]
                # sumValVect = np.sum(
                #     [(invTEval[i] - 1) * sqHess[i] for i in range(innerLoopIter)]
                # )
                #
                # postErrCov = np.identity(noOfLonPoints) + sumValVect
                # totalVar = np.trace(postErrCov)
                # print(f"Loop {outer}, Posterior total S.D. = {np.sqrt(totalVar)}")
                #
                # totalVarOuter[outer] = totalVar.copy()

                # Update outer-loop solar wind speed array
                if useLogTrans:
                    # plt.plot(range(noOfLonPoints), BlogHalf.dot(chiIter[outer + 1, :]))
                    # vLogIter[outer + 1, 0, :] = vLogIter[outer, 0, :] + BlogHalf.dot(chiIter[outer + 1, :])
                    # vOuterIter[outer + 1, 0, :] = 200 + np.exp(vLogIter[outer + 1, 0, :])
                    vOuterIter[outer + 1, 0, :] = (
                            vOuterIter[outer, 0, :] - np.exp(vLogIter[outer, 0, :])
                            + np.exp(vLogIter[outer, 0, :] + BlogHalf.dot(chiIter[outer + 1, :]))
                    )
                    # print(vOuterIter[outer, 0, :])
                    # print(200 + np.exp(BlogHalf.dot(chiIter[outer + 1, :])+np.log(vOuterIter[outer, 0, :]-200)))
                    # sys.exit()
                    # vOuterIter[outer + 1, 0, :] = (
                    #         200 + np.exp(BlogHalf.dot(chiIter[outer + 1, :])+np.log(vOuterIter[outer, 0, :]-200))
                    # )
                else:
                    vOuterIter[outer + 1, 0, :] = vOuterIter[outer, 0, :] + Bhalf.dot(chiIter[outer + 1, :])

                print(np.abs(vOuterIter[outer + 1, 0, :] - vOuterIter[outer, 0, :]).max())
                for rIndex in range(1, noOfRadPoints):
                    vOuterIter[outer + 1, rIndex, :] = bme.forwardRadModelNoLat(
                        vOuterIter[outer + 1, rIndex - 1, :], vOuterIter[outer + 1, 0, :], r[rIndex - 1], rIndex - 1,
                        deltaRrs, deltaPhi, solarRotFreq, alpha, rH, noOfLonPoints
                    )
                fig = plt.figure()
                plt.plot(range(noOfLonPoints), vOuterIter[0, obsPosDf.loc["ACE"]["radCoord"][0], :], color='r')
                # plt.plot(range(noOfLonPoints), obsCompDf.loc["STERA"]["y"][:], color='grey')
                # plt.plot(range(noOfLonPoints), obsCompDf.loc["STERB"]["y"][:], color='lightgrey')
                plt.plot(range(noOfLonPoints), obsCompDf.loc["ACE"]["y"][:], color='k')
                plt.plot(range(noOfLonPoints), vOuterIter[outer + 1, obsPosDf.loc["ACE"]["radCoord"][0], :])
                plt.show()
                # Calculate non-linear cost function
                print("Final non-linear cost function")
                costFuncOuterVal[outer + 1] = bme.calcCostFuncNoLat(
                    np.transpose(Bhalf).dot(Bhalf), R, H, vOuterIter[outer + 1, 0, :],
                    vPrior[w, 0, :],  vOuterIter[outer + 1, :, :], y,  radObs, nRadObs
                )
                print("")
                # bme.calcNonlinCostFuncPrecond(
                #     chi[-1, :], chiB[outer, :], R, H, vOuterIter[outer + 1, :, :], y, radObs, nRadObs,
                #     printCF=printCost
                # )
                vIter[1, :, :] = vOuterIter[outer + 1, :, :]

                #chiIter[outer + 1, :] = chi[-1, :].copy()
                chiB[outer + 1, :] = chiB[outer, :] - chiIter[outer + 1, :]

                # Check whether tolerance is exceeded
                if np.abs(costFuncOuterVal[outer + 1] - costFuncOuterVal[outer]) < gTol:
                    break

            # Place final outer-loop velocity into vIter array
            #vIter[1, 0, :] = vOuterIter[-1, 0, :]

        ####################END PRECONDITIONED DA####################################################
        else:
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
            resOpt = scipy.optimize.minimize(
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
        # Initialise dictionaries to hold plotting data at each observation location
        vPlotDf = pd.DataFrame(columns=["observations", "ensMean", "prior", "posterior"])
        rmseDf = pd.DataFrame(columns=["ensMean", "prior", "posterior"])

        vPriorObs = np.zeros((noOfLonPoints))
        vPostObs = np.zeros((noOfLonPoints))
        vEnsObs = np.zeros((noOfLonPoints))

        for obsName in obsCompDf.index:
            radCoordReq = radCoordDict[obsName]
            lonCoordReq = lonCoordDict[obsName]

            vPriorTemp = np.copy(vIter[0, radCoordReq, ::-1])
            vPostTemp = np.copy(vIter[-1, radCoordReq, ::-1])
            vEnsTemp = np.copy(vMASMean[w, radCoordReq, ::-1])
            vObsTemp = np.copy(obsCompDf.loc[obsName]["yPlot"][::-1])

            for i in range(noOfLonPoints):
                # Extract relevant solar wind velocities at STEREO A location
                obsTrans = int(np.mod(i + lonCoordReq, noOfLonPoints))
                vPriorObs[obsTrans] = np.copy(vPriorTemp[i])
                vPostObs[obsTrans] = np.copy(vPostTemp[i])
                vEnsObs[obsTrans] = np.copy(vEnsTemp[i])

            vPlotDf.loc[obsName] = [vObsTemp, vEnsObs, vPriorObs, vPostObs]

            # Calculate RMSEs
            rmsePriorTemp = bme.calcStateObsRMSE(vPriorObs, vObsTemp)
            rmsePostTemp = bme.calcStateObsRMSE(vPostObs, vObsTemp)
            rmseEnsTemp = bme.calcStateObsRMSE(vEnsObs, vObsTemp)

            # Append RMSEs into rmseAllWinDf
            rmseDf.loc[obsName] = [rmseEnsTemp, rmsePriorTemp, rmsePostTemp]
            for col in rmseDf.columns:
                try:
                    rmseAllWinDf.loc[obsName][col].append(rmseDf.loc[obsName][col])
                except:
                    rmseAllWinDf.loc[obsName] = [[x] for x in rmseDf.loc[obsName][:]]
                    break

            if makePlots:
                swFileDir = os.path.join(outputDir, 'plots', 'swOver27d')
                fig, ax = bme.plotSWspeedDict(
                    deltaPhiDeg, vPlotDf.loc[obsName], swFileDir, obsName, currMJD[w], fontSize=18, lWid=2.0
                )

        ###########################################################################
        # Output RMSEs for user
        ###########################################################################
        print('------------------------------------------------------')
        print('Window no. ' + str(w))
        print(' ')

        outFile.write('------------------------------------------------------\n')
        outFile.write(f'Window no. {w} \n\n')
        for obsName in rmseAllWinDf.index:
            print(f'RMSE {obsName} MASMean = {rmseAllWinDf.loc[obsName]["ensMean"][w]}')
            print(f'RMSE {obsName} Prior = {rmseAllWinDf.loc[obsName]["prior"][w]}')
            print(f'RMSE {obsName} Posterior = {rmseAllWinDf.loc[obsName]["posterior"][w]}')
            print(' ')

            outFile.write(f'RMSE {obsName} MASMean = {rmseAllWinDf.loc[obsName]["ensMean"][w]}')
            outFile.write(f'RMSE {obsName} Prior = {rmseAllWinDf.loc[obsName]["prior"][w]}')
            outFile.write(f'RMSE {obsName} Posterior = {rmseAllWinDf.loc[obsName]["posterior"][w]}')
            outFile.write(' ')

            ######################################################################################
            # Write ACE, STEREOA and STEREOB obs and prior, post, MAS ens in obs loc.
            # during solar rotation to file for later use
            # Data output is ASCENDING IN TIME
            ######################################################################################
            outObsFile = os.path.join(
                outputDir, obsName, f'obs_MJDstart{int(currMJD[w])}.txt'
            )
            with open(outObsFile, 'w') as fOut:
                np.savetxt(fOut, vPlotDf.loc[obsName]["observations"])

            outEnsFile = os.path.join(
                outputDir, obsName, f'ensMean_MJDstart{int(currMJD[w])}.txt'
            )
            with open(outEnsFile, 'w') as fOut:
                np.savetxt(fOut, vPlotDf.loc[obsName]["ensMean"])

            outPriorFile = os.path.join(
                outputDir, obsName, f'prior_MJDstart{int(currMJD[w])}.txt'
            )
            with open(outPriorFile, 'w') as fOut:
                np.savetxt(fOut, vPlotDf.loc[obsName]["prior"])

            outPostFile = os.path.join(
                outputDir, obsName, f'post_MJDstart{int(currMJD[w])}.txt'
            )
            with open(outPostFile, 'w') as fOut:
                np.savetxt(fOut, vPlotDf.loc[obsName]["posterior"])

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
