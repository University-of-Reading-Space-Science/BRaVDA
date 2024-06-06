import os
import random as rdm
import datetime
import runBravdaExport as rbe
import pandas as pd




def bravdafunction(forecastdate, configFile, obsToAssim = 'C', usecustomens = True,
                   runoutputdir = '', plottimeseries = True, corona = 'MAS',
                   precondState = False, useLogTrans = False):
    rdm.seed(20000)
    #################################
    # Data to be provided by user
    #################################
    # Find current directory
    currentDir = os.path.dirname(os.path.abspath(rbe.__file__))
    #make this the owrking directory?


    
    ###########################################################
    # Define relevant file/directories to runBravda function
    ###########################################################
    # Specify location of config file containing all necessary files and directories
    # (relative to the location of the currentDir)
    # configFile = os.path.join(currentDir, 'configBravda.dat')
    
    # Specify location of file containing necessary HUX variables
    if corona == 'MAS':
        huxVarFile = os.path.join(currentDir, 'huxVar_MAS.dat')
    elif corona == 'WSA':
        huxVarFile = os.path.join(currentDir, 'huxVar_WSA.dat')
    
    # Specify output directory
    if runoutputdir:
        outputDir = runoutputdir
    else:
        outputDir = os.path.join(currentDir, 'BRaVDA_test', 'Solar_max_test')
    
    ####################################################################
    # Choose which date to run BRaVDA from and for how many consecutive
    # 27 day periods
    ####################################################################
    # Initial date (inclusive) (string in the form DDMMYYYY required)
    initDateTime = forecastdate - datetime.timedelta(days=27.27)
    initDate = initDateTime.strftime('%d%m%Y')
    
    # Number of consecutive 27-day periods(/solar Rotations) to run
    noOfConsecWindows = 1
    #############################################################################
    
    ########################################################################################
    # Define variables for the generation of the prior and prior error covariance matrix, B
    ########################################################################################
    # Number of MAS ensemble members to generate prior
    # (default is 501 members)
    noOfMASens = 501
    
    # Give localisation radius for generating B as a function of longitude (in degrees); 0= no localisation
    # (default is 15 degrees)
    locRad = 15.0
    ##########################################################################################
    
    #################################################################
    # Set up which obs are assimilated and how R is set-up
    #################################################################
    # Specify which observations are to be assimilated
    # Write a string of the letters A,B and/or C: A=STERA, B=STERB, C=ACE)
    # eg. ABC means assimilated all obs. streams, AC means assimilated STER A and ACE only
    #obsToAssim = ["ACE"]
    
    # Specify whether R be generated as a percentage of the prior solar wind or as a constant supplied by the user
    # (Write B for background solar wind or C for constant standard deviation, then provide 3 constants seperated
    # by a space indicating the percentage of background (0,1]
    # or constant>0 for STERA, STERB and ACE [even if not assim.])
    setupOfR = 'B 0.1 0.1 0.1'

    ############################################################################
    # Set up tolerance in minimisation function in DA
    # gtol for minimisation algorithm.
    # Default is 1e-5, set to values>1 to skip over the minimisation algorithm
    ############################################################################
    gTol = 1e-5
    
    ############################################################################
    # Specify whether you wish to make plots of solar wind speed over 27 days
    ############################################################################
    makePlots = plottimeseries

    rbe.runBravDA(configFile, huxVarFile, outputDir, obsToAssim, setupOfR, 
                  initDate, noOfConsecWindows, noOfMASens,
                  locRad, gTol, makePlots, usecustomens=usecustomens, 
                  precondState=precondState, useLogTrans=useLogTrans)

if __name__=="__main__":

    start_date = datetime.datetime(2023, 4, 14)
    end_date = datetime.datetime(2024, 4, 14)
    dates = pd.date_range(start_date, end_date, freq='D')

    # Specify forecastDate
    fcDate = datetime.datetime(2013, 8, 4)
    currentDir = os.path.dirname(os.path.abspath(rbe.__file__))
    configFile = os.path.join(currentDir, 'configBravda.dat')
    precondState = True
    useLogTrans = False
    """
    for date in dates:
        # Run BRaVDA function
        bravdafunction(
            date, obsToAssim=["OMNI"], usecustomens=False, runoutputdir='', plottimeseries=True,
            corona='MAS', precondState=precondState, useLogTrans=useLogTrans
        )
    """

    bravdafunction(
        fcDate, configFile, obsToAssim=["STERA", "STERB", "OMNI"], usecustomens=False, runoutputdir='',
        plottimeseries=True, corona='MAS', precondState=precondState, useLogTrans=useLogTrans
    )