# -*- coding: utf-8 -*-
# Import necessary packages
import numpy as np
import os
import matplotlib.pyplot as plt
import sys


def gregToMJD(year, dayOfYear, hour):
    # Convert between Gregorian year, day of year and hour to MJD
    # Find month and day of month from day of year
    # Define length of months
    if np.mod(year, 4) == 0:
        # dayMonths=[31,29,31,30,31,30,31,31,30,31,30,31]
        daysOfStartMonth = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366]
    else:
        # dayMonths=[31,28,31,30,31,30,31,31,30,31,30,31]
        daysOfStartMonth = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]

    month = 0
    a = 0
    m = 1
    while (a == 0) or (m == 14):
        if m == 13:
            print('m=13: While loop overdone: Day no. not found')
            sys.exit()
        elif daysOfStartMonth[m - 1] < dayOfYear <= daysOfStartMonth[m]:
            month = m
            a = 1

        m = m + 1

    day = dayOfYear - daysOfStartMonth[month - 1]

    # Calculate Julian Date
    jdTerm1 = (367 * year)
    jdTerm2 = np.floor(7 * (year + np.floor((month + 9) / 12.0)) * 0.25)
    jdTerm3 = np.floor(0.75 * (np.floor(0.01 * (year + ((month - 9) / 7.0))) + 1))
    jdTerm4 = np.floor((275 * month) / 9.0)
    jdTerm5 = day + 1721028.5

    julianDayObs = jdTerm1 - jdTerm2 - jdTerm3 + jdTerm4 + jdTerm5

    # julianDayObs = (
    #        (367 * year) - np.floor(7 * (year + np.floor((month + 9) / 12.0)) * 0.25) - (
    #    np.floor(0.75 * (np.floor(0.01 * (year + ((month - 9) / 7.0))) + 1))
    # ) + np.floor((275 * month) / 9.0) + day + 1721028.5
    # )

    julianTimeObs = julianDayObs + (hour / 24.0)
    MJD = julianTimeObs - 2400000.5

    return MJD


def dateToMJD(date):
    # Convert between date (string in form DDMMYYYY) to MJD
    day = int(date[:2])
    month = int(date[2:4])
    year = int(date[4:])

    if not 0 < day < 32:
        print('Invalid day input in date, please check.')
        print(f'Day input = {day}')
        sys.exit()

    if not 0 < month < 13:
        print('Invalid month input in date')
        print(f'Month input = {month}')
        sys.exit()

    jdTerm1 = (367 * year)
    jdTerm2 = np.floor(7 * (year + np.floor((month + 9) / 12.0)) * 0.25)
    jdTerm3 = np.floor(0.75 * (np.floor(0.01 * (year + ((month - 9) / 7.0))) + 1))
    jdTerm4 = np.floor((275 * month) / 9.0)
    jdTerm5 = day + 1721028.5

    julianDayObs = jdTerm1 - jdTerm2 - jdTerm3 + jdTerm4 + jdTerm5
    MJD = julianDayObs - 2400000.5

    return MJD


def readObsFile(obsFile, mjdFile):
    # Extract Nearest Neighbour observations from obsFile
    # at the values specified in mjdFile

    # Read obs file
    with open(obsFile, 'r') as fObs:
        vObsLines = fObs.readlines()

    with open(mjdFile, 'r') as MJDfiles:
        MJDvalues = MJDfiles.readlines()

    # Reverse the MJD values
    MJDvalues.reverse()

    # Initialise variables
    obsAll = []
    obsMJDtime = []

    # Extract all observations from file
    for ln in vObsLines:
        s = ln.split()

        # Calculate time of observation and convert to MJD
        year = int(s[0].strip())
        doY = int(s[1].strip())
        hr = int(s[2].strip())
        MJDtime = gregToMJD(year, doY, hr)

        # Append MJDtime to obsMJDtime
        obsMJDtime.append(MJDtime)

        # Extract all observations and add to obsAll
        obsValue = float(s[3].strip())
        obsAll.append(obsValue)

    # Extract length of observations list
    lenObs = len(obsAll)

    # Filter to MJD times required (nearest neighbour)
    observations = []
    i = 0

    for lnMJD in MJDvalues:
        # i=j
        # If obsMJDtime is less than required MJD, move to next one until obsMJDtime is greater
        reqMJD = float(lnMJD.strip())
        while (obsMJDtime[i] < reqMJD) and (i < (lenObs + 2)):
            if i == (lenObs + 1):
                print('Date out of bounds')
                sys.exit()

            i = i + 1

        # Extract the observation nearest to required MJD
        observations.append(obsAll[i])

    return observations


def readObsFileAvg(obsFile, mjdFile, initTime, finalTime):
    # Extract average observations from obsFile
    # around the values specified in mjdFile
    # Read file containing the observations
    with open(obsFile, 'r') as fObs:
        vObsLines = fObs.readlines()

    # Read in the MJD values we wish to extract the obs. at
    with open(mjdFile, 'r') as MJDfiles:
        MJDvalues = MJDfiles.readlines()

    # Reverse the MJD values to get them in ascending order
    MJDvalues.reverse()
    noOfObsReq = len(MJDvalues)

    # Initialise variables
    obsAll = []
    obsMJDtime = []

    # Extract all observations from file
    for ln in vObsLines:
        s = ln.split()

        # Calculate time of observation
        year = int(s[0].strip())
        doY = int(s[1].strip())
        hr = int(s[2].strip())
        MJDtime = gregToMJD(year, doY, hr)

        # Append MJDtime to obsMJDtime
        obsMJDtime.append(MJDtime)

        # Extract all observations and add to obsAll
        obsValue = float(s[3].strip())
        obsAll.append(obsValue)

    # Extract length of observations list
    lenObs = len(obsAll)

    # Filter to MJD times required (averaging)
    observations = np.ones(noOfObsReq) * 9999
    cntMin = 0

    # Calculate the start and end times for the averaging period of each observation
    for k in range(noOfObsReq):
        if k == 0:
            # First observation MJD time
            # Extract first and second obs. MJD values
            MJDvalCurr = float(MJDvalues[0].strip())
            MJDvalNext = float(MJDvalues[1].strip())

            # Calculate half the time between second and first obs. MJD
            MJDNextDiff = (MJDvalNext - MJDvalCurr) / 2

            # Make start of average period either the initial model time or the the first obs. time
            # minus the mid-point of first and second MJD obs. time (i.e. make it symmetrical)
            startAvgPeriod = max(initTime, MJDvalCurr - MJDNextDiff)
            endAvgPeriod = (MJDvalNext + MJDvalCurr) / 2

        elif k == (noOfObsReq - 1):
            # Final observation MJD time

            # Extract penultimate and last obs. MJD values
            MJDvalPrev = float(MJDvalues[-2].strip())
            MJDvalCurr = float(MJDvalues[-1].strip())

            # Calculate mid-point of penultimate and final obs. time
            MJDPrevDiff = (MJDvalCurr - MJDvalPrev) / 2

            # Make start of average period either the final model time or the the last obs. time
            # minus the halved difference between penultimate and last MJD obs. time (i.e. make it symmetrical)
            startAvgPeriod = (MJDvalPrev + MJDvalCurr) / 2
            endAvgPeriod = min(MJDvalCurr + MJDPrevDiff, finalTime)
        else:
            # Extract previous, current and next obs. MJD values
            MJDvalPrev = float(MJDvalues[k - 1].strip())
            MJDvalCurr = float(MJDvalues[k].strip())
            MJDvalNext = float(MJDvalues[k + 1].strip())

            # Calculate mid-point between current obs. MJD and previous/next obs. MJD
            startAvgPeriod = (MJDvalPrev + MJDvalCurr) / 2
            endAvgPeriod = (MJDvalCurr + MJDvalNext) / 2

        # This loop will exit either when the obsMJD time is first greater than the start of the averaging period
        # or we've reached the end of the obs. file, in which case it breaks out of the loop
        while (obsMJDtime[cntMin] < startAvgPeriod) and (cntMin < (lenObs + 2)):
            if cntMin == (lenObs + 1):
                print('Date out of bounds')
                print(f'No observations for date requested in {obsFile}')
                print(f'Final date in obs. file is MJD: {MJDvalCurr}')
                print(f'Date requested: {startAvgPeriod}')
                print('Filling with 9999.9')
                observations[k] = 9999.9
                break

            cntMin = cntMin + 1

        # Once the start of the averaging period has been found,
        # we need to find the end of the averaging period in same way
        # As end of averaging period must be greater than start,
        # start upper counter from start of averaging period
        cntMax = cntMin
        while (obsMJDtime[cntMax] < endAvgPeriod) and (cntMax < (lenObs + 2)):
            if cntMax == (lenObs + 1):
                print('Date out of bounds')
                print(f'No observations for date requested in {obsFile}')
                print(f'Final date in obs. file is MJD: {MJDvalCurr}')
                print(f'Date requested: {endAvgPeriod}')
                print('Filling with 9999.9')
                observations[k] = 9999.9
                break

            cntMax = cntMax + 1

        # Now find the average observation value between startAvgPeriod and endAvgPeriod
        obsToAvg = obsAll[cntMin:cntMax]

        # Remove any erroneous values from obsToAvg (eg. 9999.9's)
        obsToAvg = [np.nan if (x > 5000) or (x < 0) else x for x in obsToAvg]

        if obsToAvg:
            observations[k] = np.nanmean(obsToAvg)
        else:
            observations[k] = 9999.9

    # Filter out any NaNs and replace with 9999.9 for BRaVDA to filter later
    observations = [9999.9 if np.isnan(x) else x for x in observations]

    return observations


def readObsLocFileWindows(obsLocFile, mjdCRFile, noOfWindows):
    # Extract the observations locations at the mid-point of the
    # window

    # Read obs. location file
    with open(obsLocFile, 'r') as fObs:
        vObsLines = fObs.readlines()

    # Read mjd values required
    with open(mjdCRFile, 'r') as MJDfiles:
        MJDvalues = MJDfiles.readlines()

    obsLon = []
    obsRad = []
    obsLat = []
    obsMJDtime = []

    # Extract all observation locations from obsLocFile
    for ln in vObsLines:
        s = ln.split()
        if s[0].strip() != 'YEAR':
            year = int(s[0].strip())
            doY = int(s[1].strip())
            MJDtime = gregToMJD(year, doY, 0)

            obsMJDtime.append(float(MJDtime))
            obsLon.append(float(s[4].strip()))
            obsRad.append(float(s[2].strip()))
            obsLat.append(float(s[3].strip()))

    lenObsLon = len(obsLon)
    obsLocLon = []
    obsLocRad = []
    obsLocLat = []
    i = 0

    # Extract required locations
    for lnMJD in MJDvalues:

        MJDsplit = lnMJD.split()
        winNo = int(MJDsplit[0].strip())
        obsTimeReq = float(MJDsplit[1].strip())

        if 0 <= winNo <= noOfWindows:
            while (obsMJDtime[i] < obsTimeReq) and (i < lenObsLon + 2):
                if i == lenObsLon + 1:
                    print('Date out of bounds')
                    sys.exit()

                i = i + 1

            # Remake obs lon, so it goes from -180 to 180deg
            if obsLon[i] > 180:
                obsLon[i] = obsLon[i] - 360

            obsLocLon.append(float(obsLon[i]))
            obsLocRad.append(float(obsRad[i]))
            obsLocLat.append(float(obsLat[i]))

    # Package up for returning
    obsOut = np.zeros((3, len(obsLocLon)))
    obsOut[0, :] = obsLocLon
    obsOut[1, :] = obsLocRad
    obsOut[2, :] = obsLocLat

    return obsOut


def getObsPos(
        obsLocFile, earthLocFile, mjdCRFile, noOfWindows,
        rS, innerRadRs, deltaRrs, deltaPhiDeg
):
    # Read in observation positions
    readObsPos = bme.readObsLocFileWindows(obsLocFile, mjdCRFile, noOfWindows)
    readEarthPos = bme.readObsLocFileWindows(earthLocFile, mjdCRFile, noOfWindows)

    # Extract difference between Earth and observation locations
    obsEarthPosDiff = readEarthPos - readObsPos

    # Extract radial positions
    obsRadAU = readObsPos[1, :]

    # Convert radial positions from AU to solar radii
    obsRadRs = obsRadAU * 215
    obsRadKm = obsRadRs * rS

    # Get radial coordinates of obs.
    obsRadCoord = ((obsRadKm - innerRadRs) / deltaRrs - 1).round().astype(int)

    # Extract longitude positions
    obsLon = np.mod(360 * np.ones(noOfWindows) - obsEarthPosDiff[0, :], 360)

    # Set up model coordinates for the longitudes
    obsLonCoord = (obsLon / deltaPhiDeg).round().astype(int)

    return obsRadAU, obsRadRs, obsRadKm, obsRadCoord, obsLon, obsLonCoord


def makeMJDfile(MJDstart, noOfLon, fileMJDOut, daySolarRot=27.2753):
    # Function to make MJD file for one solar rotation from MJD start of length daySolarRot

    # Define timestep as length of solar rotation divided by number of lon points during rotation
    MJDstep = daySolarRot / float(noOfLon)

    # Initialise MJD output variable with MJD startdate provided by user
    MJD = np.zeros(noOfLon)
    MJD[0] = np.copy(MJDstart)

    # Calculate evenly spread MJD values over solar rotation
    for lo in range(noOfLon - 1):
        MJD[lo + 1] = MJD[lo] + MJDstep

    # Write MJD values to required output file
    with open(fileMJDOut, "w") as MJDfile:
        for lo in range(noOfLon - 1, -1, -1):
            MJDfile.write(str(MJD[lo]) + '\n')

    return None


def createUnpertAndBMatrix(
        deltaPhiDeg, locRad, noOfLonPoints, ensFile,
        initMJDCR, currMJD, noOfMASens
):
    # Read ensemble file, extract the unperturbed first member to use as the prior
    # and generate the B matrix with the rest

    # Read input ensemble file
    vIn = np.loadtxt(ensFile)
    '''plt.plot(range(len(vIn[0, :])), vIn[0, :])
    plt.show()
    print(np.shape(vIn))'''

    #sys.exit()
    ################################
    # Rotate vIn to required Carrington longitude
    ###############################
    rotLonCoord = int(np.round(128 * (1 - ((currMJD - initMJDCR) / 27.2753))))
    vIn[:, :] = np.append(
        vIn[:, rotLonCoord:], vIn[:, :rotLonCoord], axis=1
    )
    ####################################################

    '''plt.plot(range(len(vIn[0, :])), vIn[0, :])
    plt.show()
    print(np.shape(vIn))
    sys.exit()'''
    # Extract the unperturbed ensemble member (the first one)
    unperturbEnsMem = vIn[0, :]

    # Generate perturbation ensemble (X-mean(X))
    meanEns = vIn[1:, :].mean(axis=0)
    pertEns = vIn[1:, :] - meanEns

    B = (1 / ((noOfMASens - 1) - 1)) * np.transpose(pertEns).dot(pertEns)

    '''plt.figure(3)
    im = plt.contourf(
        range(128),range(128), B,
        cmap=plt.cm.seismic,levels = range(-5700,5700,200)
    )
    cobar = plt.colorbar(im)
    cobar.set_label('km/s')
    plt.xlabel('Carrington longitude coordinate')
    plt.ylabel('Carrington longitude coordinate')
    plt.title('Localised B matrix')
    plt.show()'''

    if locRad != 0:
        locMat = np.zeros((np.shape(B)))

        for i in range(noOfLonPoints):
            for j in range(noOfLonPoints):
                dist1 = abs((i - j) * deltaPhiDeg)
                dist2 = 360 - abs((i - j) * deltaPhiDeg)
                distBetweenPoints = min(dist1, dist2)
                locMat[i, j] = np.exp(-(distBetweenPoints ** 2) / (2 * (locRad ** 2)))

        B = locMat * np.copy(B)

        '''plt.figure(4)
        im = plt.contourf(
            range(128),range(128), B,
            cmap=plt.cm.seismic,levels = range(-5700,5700,200)
        )
        cobar = plt.colorbar(im)
        cobar.set_label('km/s')
        plt.xlabel('Carrington longitude coordinate')
        plt.ylabel('Carrington longitude coordinate')
        plt.title('Localised B matrix')
        plt.show()'''
    '''plt.plot(range(len(unperturbEnsMem)), unperturbEnsMem)
    plt.show()
    print(np.shape(vIn))'''
    return unperturbEnsMem, meanEns[:noOfLonPoints], B[:noOfLonPoints, :noOfLonPoints]


def forwardRadModelNoLat(
        v, v0, r, rIndex, deltaR, deltaPhi, solarRotFreq, alpha, rH, noOfLonPoints
):
    # Run HUX model forward one timestep

    # Initialise output array
    velOut = np.zeros(noOfLonPoints)

    # Calculate CFL criteria for model run
    CFL = deltaR * solarRotFreq / deltaPhi

    # Execute the model equation for each longitude
    for phi in range(noOfLonPoints):
        # Extract required velocities from array and store in dummy variables
        vPhi = v[phi]
        vPhi_1 = v[phi - 1]
        v0Phi_1 = v0[phi - 1]

        if rIndex != 0:
            # Away from the inner boundary

            # Calculate terms in velOut
            coeffVO1 = CFL * (vPhi - vPhi_1) / vPhi_1

            # Calculate acceleration terms
            accelTerm1 = alpha * v0Phi_1 * (1 - np.exp(-(r - 1) / rH))
            accelTerm2 = alpha * v0Phi_1 * (1 - np.exp(-r / rH))

            # Calculate velOut
            velOut[phi - 1] = vPhi_1 + coeffVO1 - accelTerm1 + accelTerm2

        else:
            # If at the inner boundary
            # Calculate terms of velOut
            coeffVO1 = CFL * (vPhi - vPhi_1) / vPhi_1
            accelTerm = alpha * v0Phi_1 * (1 - np.exp(-r / rH))

            velOut[phi - 1] = vPhi_1 + coeffVO1 + accelTerm

    return velOut


def TLMRadNoLat(
        dv, v, r, rIndex, deltaR, deltaPhi, solarRotFreq, alpha, rH, noOfLonPoints
):
    # Calculate Tangent linear model for forwardRadModel

    # Initialise output TLM variable/perturbation
    dvOut = np.zeros(noOfLonPoints)

    # Calculate CFL criterion for use in model equations
    CFL = deltaR * solarRotFreq / deltaPhi

    # Run Tangent Linear Model equations
    if rIndex != 0:
        for phi in range(noOfLonPoints):
            # Extract required velocities and store in dummy variables
            vPhi_1 = v[phi - 1]
            vPhi = v[phi]

            # Extract required perturbations and store in dummy variables
            dvPhi_1 = dv[phi - 1]
            dvPhi = dv[phi]

            # Calculate coefficients for producing dvOut
            coeffVO1 = 1 - CFL * (vPhi / (vPhi_1 * vPhi_1))
            coeffVO2 = CFL * (1.0 / vPhi_1)

            # Calculate output perturbation
            dvOut[phi - 1] = (coeffVO1 * dvPhi_1) + (coeffVO2 * dvPhi)

    else:
        # If r=0, different differentiation result from vAcc component
        for phi in range(noOfLonPoints):
            # Extract required velocities and store in dummy variables
            vPhi_1 = v[phi - 1]
            vPhi = v[phi]

            # Extract required perturbations and store in dummy variables
            dvPhi_1 = dv[phi - 1]
            dvPhi = dv[phi]

            # Calculate coefficient of first term
            coeffVO11 = 1 - CFL * (vPhi / (vPhi_1 * vPhi_1))
            coeffVO12 = alpha * (1 - np.exp(-r / rH))
            coeffVO1 = coeffVO11 + coeffVO12

            # Calculate coefficient of second term
            coeffVO2 = CFL * (1.0 / vPhi_1)

            # Calculate dvOut
            dvOut[phi - 1] = (coeffVO1 * dvPhi_1) + (coeffVO2 * dvPhi)

    return dvOut


def adjRadNoLat(
        dv9, v, r, rIndex, deltaR, deltaPhi, solarRotFreq, alpha, rH, noOfLonPoints
):
    # Calculate Adjoint model of forwardRadModel

    # Initialise output adjoint state
    dv9Out = np.zeros(noOfLonPoints)  # +1))

    # Calculate CFL criterion for use in model equations
    CFL = deltaR * solarRotFreq / deltaPhi

    if rIndex != 0:
        for phi in range(noOfLonPoints):
            # Extract required velocities for adjoint equation
            vPhi = v[phi]
            vPhi_1 = v[phi - 1]
            vPhi_2 = v[phi - 2]

            # Extract required adjoint variables
            dv9Phi_1 = dv9[phi - 1]
            dv9Phi_2 = dv9[phi - 2]

            # Generate coefficients for dv9Out equation
            coeffAdj1 = 1 - CFL * (vPhi / (vPhi_1 * vPhi_1))
            coeffAdj2 = CFL * (1.0 / vPhi_2)

            # Calculate dv9Out
            dv9Out[phi - 1] = (coeffAdj1 * dv9Phi_1) + (coeffAdj2 * dv9Phi_2)

    else:
        # When r=0, there is an additional component in the v_acc that must be accounted for
        for phi in range(noOfLonPoints):
            # Extract required velocities for adjoint equation
            vPhi = v[phi]
            vPhi_1 = v[phi - 1]
            vPhi_2 = v[phi - 2]

            # Extract required adjoint variables
            dv9Phi_1 = dv9[phi - 1]
            dv9Phi_2 = dv9[phi - 2]

            # Generate coefficients for adjoint equation
            coeffAdj1 = CFL * (1.0 / vPhi_2)
            coeffAdj2 = 1 - CFL * (vPhi / (vPhi_1 * vPhi_1)) + alpha * (1 - np.exp(-r / rH))

            # Calculate dv9Out
            dv9Out[phi - 1] = (coeffAdj1 * dv9Phi_2) + (coeffAdj2 * dv9Phi_1)

    return dv9Out


def adjointTestNoLat(
        x, y, v, r, rIndex, deltaR, deltaPhi, solarRotFreq, alpha, rH, noOfLonPoints
):
    #######################################################################
    # Code to test if <Fx,y>=<x,F^T y> holds for adjoint and TLM scripts
    #######################################################################

    # Run TLM model on x and Adjoint model on y
    TLMx = TLMRadNoLat(
        x, v, r, rIndex, deltaR, deltaPhi, solarRotFreq, alpha, rH, noOfLonPoints
    )
    Adjy = adjRadNoLat(
        y, v, r, rIndex, deltaR, deltaPhi, solarRotFreq, alpha, rH, noOfLonPoints
    )

    # Reshape variables
    xReshape = x.reshape(noOfLonPoints, 1)
    yReshape = y.reshape(noOfLonPoints, 1)

    # Do dot products
    yDotTLMx = (TLMx.transpose()).dot(yReshape)
    adjyDotx = (xReshape.transpose()).dot(Adjy)

    # Calculate difference between them
    tlmAdjDiff = yDotTLMx - adjyDotx

    # Print out results of adjoint test
    print(f'<y,Fx> = {yDotTLMx}')
    print(f'<F*y,x> = {adjyDotx}')
    print(f'<y,Fx>-<F*y,x> = {tlmAdjDiff}')

    return tlmAdjDiff


def obsOp(lengthState, obsToBeTakenSpace, lonCoord):
    # Generates observation operator H_{obsTime}(lengthState,noOfObs)
    lengthObsSpace = len(obsToBeTakenSpace)
    lengthLon = len(lonCoord)
    H = np.zeros((lengthObsSpace * lengthLon, lengthState))

    # Use identity to start with
    for lo in range(lengthLon):
        for i in range(lengthObsSpace):
            j = obsToBeTakenSpace[i]
            ind1 = int((lo * lengthObsSpace) + i)
            ind2 = int(np.mod(j + lonCoord[lo], lengthState))
            H[ind1, ind2] = 1

    return H


def calcCostFuncNoLat(
        B, R, H, x, xb, v, y, radObs, nRadObs
):
    #################################################################
    # Function to calculate the cost function and print value
    #################################################################

    # Initialise variables
    costFuncObs = 0

    # Calculate background component of cost function
    backDiff = x - xb
    Binv = np.linalg.pinv(B)
    costFuncBack = np.transpose(backDiff).dot(Binv).dot(backDiff)

    # Calculate observational component of cost function
    for r in range(len(radObs)):
        if r == 0:
            # Extract required variables for current observation
            obsRadReq = int(radObs[0])
            nrUp = int(nRadObs[0])

            # Calculate the components of the observational component of cost function
            innov = y[0:nrUp] - H[0:nrUp, :].dot(v[obsRadReq, :])
            Rinv = np.linalg.pinv(R[0:nrUp, 0:nrUp])

            # Calculate observational component of cost function for FIRST observation
            costFuncObs = np.transpose(innov).dot(Rinv).dot(innov)
        else:
            # Extract required variables for current observation
            obsRadReq = int(radObs[r])
            nrDown = int(nRadObs[r - 1])
            nrUp = int(nRadObs[r])

            # Calculate the components of the observational component of cost function
            innov = y[nrDown:nrUp] - H[nrDown:nrUp, :].dot(v[obsRadReq, :])
            Rinv = np.linalg.pinv(R[nrDown:nrUp, nrDown:nrUp])

            # Add current observational component to costFuncobs to build the sum
            costFuncObs = costFuncObs + (
                np.transpose(innov).dot(Rinv).dot(innov)
            )

    # Add together the two components of the cost function
    costFunc = costFuncBack + costFuncObs

    # Print the cost function's components values and cost function's value
    print(f'Jb = {costFuncBack}')
    print(f'Jo = {costFuncObs}')
    print(f'costFunc = {costFunc}')

    return costFunc


def calcCostFuncForCGNoLat(
        xb, B, R, H, xBIter0, y, radObs, nRadObs, r, rH,
        deltaR, deltaPhi, alpha, solarRotFreq, noOfRadPoints, noOfLonPoints
):
    ############################################################################
    # Function to calculate the cost function for use in minimisation algorithm
    # (includes forwardRadModel run within it)
    ############################################################################
    # Initialise variables
    costFuncObs = 0

    # Calculate background component of cost function
    backDiff = xb - xBIter0
    Binv = np.linalg.pinv(B)
    costFuncBack = np.transpose(backDiff).dot(Binv).dot(backDiff)

    # Initialise solar wind speed variable
    velNonLin = np.zeros((noOfRadPoints, noOfLonPoints))
    velNonLin[0, :] = xb

    # Define velocity at initial radius
    vInitRad = velNonLin[0, :]

    # Run forward model
    for rIndex in range(1, noOfRadPoints):
        prevRadCoord = rIndex - 1
        vPrevRad = velNonLin[prevRadCoord, :]
        prevRad = r[prevRadCoord]

        velNonLin[rIndex, :] = forwardRadModelNoLat(
            vPrevRad, vInitRad, prevRad, prevRadCoord, deltaR, deltaPhi,
            solarRotFreq, alpha, rH, noOfLonPoints
        )

    lenRo = len(radObs)
    # Calculate observational component of cost function
    for i in range(lenRo):
        if i == 0:
            # Extract required variables for current observation
            obsRadReq = int(radObs[0])
            nrUp = int(nRadObs[0])

            # Calculate the components of the observational component of cost function
            innov = y[:nrUp] - H[:nrUp, :].dot(velNonLin[obsRadReq, :])
            Rinv = np.linalg.pinv(R[:nrUp, :nrUp])

            # Calculate observational component of cost function for FIRST observation
            costFuncObs = np.transpose(innov).dot(Rinv).dot(innov)

        else:
            # Extract required variables for current observation
            obsRadReq = int(radObs[i])
            nrDown = int(nRadObs[i - 1])
            nrUp = int(nRadObs[i])

            # Calculate the components of the observational component of cost function
            innov = y[nrDown:nrUp] - H[nrDown:nrUp, :].dot(velNonLin[obsRadReq, :])
            Rinv = np.linalg.pinv(R[nrDown:nrUp, nrDown:nrUp])

            # Add current observational component to costFuncobs to build the sum
            costFuncObs = costFuncObs + np.transpose(innov).dot(Rinv).dot(innov)

    # Add together the two components of the cost function
    costFunc = costFuncBack + costFuncObs

    return costFunc


def makeGradCGNoLat(
        xb, B, R, H, xBIter0, y, radObs, nRadObs, r, rH, deltaR, deltaPhi, alpha, solarRotFreq,
        noOfRadPoints, noOfLonPoints
):
    ##############################################################################
    # Function that make the gradient of the cost function, for use within the
    # minimisation algorithm
    ###############################################################################

    # Define variable to hold solar wind speed
    velNonLin = np.zeros((noOfRadPoints, noOfLonPoints))

    # Initialise the solar wind speed variable with the background solar wind
    velNonLin[0, :] = xb

    # Initialise adjoint variable
    dv9 = np.zeros((noOfRadPoints + 1, noOfLonPoints))

    # Extract initial radius' velocity for adjoint calculations
    vInitRad = velNonLin[0, :]

    # Run forward model
    for rIndex in range(1, noOfRadPoints):
        prevRadCoord = rIndex - 1
        vPrevRad = velNonLin[prevRadCoord, :]
        prevRad = r[prevRadCoord]

        velNonLin[rIndex, :] = forwardRadModelNoLat(
            vPrevRad, vInitRad, prevRad, prevRadCoord, deltaR, deltaPhi,
            solarRotFreq, alpha, rH, noOfLonPoints
        )

    # Run adjoint model
    for rIndex in range(noOfRadPoints - 1, 0, -1):
        vCurrRad = velNonLin[rIndex]

        if rIndex in radObs:  # noOfRadPoints-1: #radObsLocation.any()==rIndex:
            i = list(radObs).index(rIndex)

            if i == 0:
                # Extract required variables and store in dummy variables
                nrUp = int(nRadObs[0])
                dv9Plus1 = dv9[rIndex + 1, :]

                # Calculate first term of adjoint equation by running adjoint model
                dv9Coeff1 = adjRadNoLat(
                    dv9Plus1, vCurrRad, r[rIndex], rIndex,
                    deltaR, deltaPhi, solarRotFreq, alpha, rH,
                    noOfLonPoints
                )

                # Calculate the observational addition to the adjoint equation
                innov = y[:nrUp] - H[:nrUp, :].dot(vCurrRad)
                Rinv = np.linalg.pinv(R[:nrUp, :nrUp])
                dv9Coeff2 = np.transpose(H[:nrUp, :]).dot(Rinv).dot(innov)

                # Sum the two terms to complete update the adjoint variable
                dv9[rIndex, :] = dv9Coeff1 + dv9Coeff2

            else:
                # Extract required variables and store in dummy variables
                nrDown = int(nRadObs[i - 1])
                nrUp = int(nRadObs[i])
                dv9Plus1 = dv9[rIndex + 1, :]

                # Calculate first term of adjoint equation by running adjoint model
                dv9Coeff1 = adjRadNoLat(
                    dv9Plus1, vCurrRad, r[rIndex], rIndex,
                    deltaR, deltaPhi, solarRotFreq, alpha, rH,
                    noOfLonPoints
                )

                # Calculate the observational addition to the adjoint equation
                innov = y[nrDown:nrUp] - H[nrDown:nrUp, :].dot(vCurrRad)
                Rinv = np.linalg.pinv(R[nrDown:nrUp, nrDown:nrUp])
                dv9Coeff2 = np.transpose(H[nrDown:nrUp, :]).dot(Rinv).dot(innov)

                # Sum the two terms to complete update the adjoint variable
                dv9[rIndex, :] = dv9Coeff1 + dv9Coeff2

        else:
            # If no observations are at the current radius, run the adjoint model
            # to get previous adjoint variable solution
            # Extract required adjoint variable
            dv9Plus1 = dv9[rIndex + 1, :]

            # Run adjoint model
            dv9[rIndex, :] = adjRadNoLat(
                dv9Plus1, vCurrRad, r[rIndex], rIndex,
                deltaR, deltaPhi, solarRotFreq, alpha, rH,
                noOfLonPoints
            )

    # Calculate the gradient of the cost function
    # Calculate component of gradient associated with the background states
    backDiff = xb - xBIter0
    Binv = np.linalg.pinv(B)
    gradJcoeff1 = Binv.dot(backDiff)

    # Calculate final adjoint model run
    gradJcoeff2 = adjRadNoLat(
        dv9[1, :], vInitRad, r[0], 0,
        deltaR, deltaPhi, solarRotFreq, alpha, rH,
        noOfLonPoints
    )

    # Calculate gradient of cost function
    gradJ = np.transpose(np.matrix(gradJcoeff1 - gradJcoeff2))

    # Extract required column and convert back to an array for output
    gradJout = np.array(gradJ)[:, 0]

    return gradJout


def calcStateObsRMSE(state, obs):
    # Calculate squared differences between state and obs.
    stateObsDiff = state - obs

    # Calculate Mean Square Errors
    mse = np.nanmean(stateObsDiff * stateObsDiff)

    # Calculate Root Mean Squared Error
    rmseOut = np.sqrt(mse)

    return rmseOut


def plotSWspeed(
        deltaPhiDeg, vPlot, outputDir, spacecraftName, currMJD, fontSize=18, lWid=2.0
):
    fig, ax = plt.subplots(
        nrows=1, ncols=1, figsize=(12, 6),
        # constrained_layout=True
    )

    # plt.figure(figsize=(12, 6))
    ax.plot(
        np.arange(0, 359, deltaPhiDeg), vPlot[0, :],
        color='k', linewidth=lWid, label=f'{spacecraftName} Data'
    )
    ax.plot(
        np.arange(0, 359, deltaPhiDeg), vPlot[1, :],
        color='m', linewidth=lWid, label='MAS ens. mean'
    )
    ax.plot(
        np.arange(0, 359, deltaPhiDeg), vPlot[2, :],
        color='b', linewidth=lWid, label='Prior'
    )
    ax.plot(
        np.arange(0, 359, deltaPhiDeg), vPlot[3, :],
        color='g', linewidth=lWid, label='Posterior'
    )
    ax.set_xlabel('Time (days)', fontsize=fontSize)
    ax.set_ylabel('Speed (km/s)', fontsize=fontSize)
    ax.set_xlim(0, 360)
    ax.set_ylim(240, 800)

    # Set fontsize for axis ticks
    ax.tick_params(axis='both', which='major', labelsize=fontSize)
    ax.tick_params(axis='both', which='minor', labelsize=int(0.8 * fontSize))

    # Set x and y ticks
    ax.set_yticks(np.arange(300, 801, 100))
    ax.set_yticklabels(['300', '400', '500', '600', '700', '800'])

    ax.set_xticks(np.append(np.arange(0, 361, 53.3), 360))
    ax.set_xticklabels(['0', '4', '8', '12', '16', '20', '24', '27'])

    # ax.set_xticks(
    #    np.append(np.arange(0, 361, 53.3), 360), ('0', '4', '8', '12', '16', '20', '24', '27'),
    #    fontsize=fontSize)
    # if w == 0:
    plt.title(fr'Solar wind speed at {spacecraftName}', fontsize=fontSize)
    ax.legend()
    # plt.tight_layout()
    saveFileLoc = os.path.join(outputDir, f'SWspeed{spacecraftName}_MJDstart{int(currMJD)}.pdf')
    plt.savefig(saveFileLoc)
    # plt.show()

    return fig, ax
