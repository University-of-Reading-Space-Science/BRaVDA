import os
import random as rdm
import matplotlib.pyplot as plt
from datetime import datetime as dt
import runBravdaExport as rbe

rdm.seed(20000)
plt.close('all')

#################################
# Data to be provided by user
#################################
# Find current directory
currentDir = os.path.dirname(os.path.abspath(__file__))

###########################################################
# Define relevant file/directories to runBravda function
###########################################################
# Specify location of config file containing all necessary files and directories
# (relative to the location of the currentDir)
configFile = os.path.join(currentDir, 'configBravda.dat')

# Specify location of file containing necessary HUX variables
huxVarFile = os.path.join(currentDir, 'huxVar.dat')

# Specify output directory
outputDir = os.path.join(currentDir, 'testBC')

####################################################################
# Choose which date to run BRaVDA from and for how many consecutive
# 27 day periods
####################################################################
# Initial date (inclusive) (string in the form DDMMYYYY required)
initDateTime = dt(2010, 8, 10)
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
obsToUse = 'BC'

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
makePlots = True

rbe.runBravDA(configFile, huxVarFile, outputDir, obsToUse, setupOfR,
              initDate, noOfConsecWindows, noOfMASens, locRad, gTol, makePlots)
