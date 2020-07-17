#Import necessary packages
from pylab import *
from sympy import Symbol, pprint
from collections import *
from random import *
from scipy import optimize, linalg
import os
import shutil
from matplotlib import *
from matplotlib.pyplot import *
import subprocess
from math import *
import sys
import pickle
import json

from STEREOmethodsExport import *

seed(20000)
close('all')

#Find current directory
currentDir=os.path.dirname(os.path.abspath(__file__))
print('Parameter file: '+str(currentDir+'\\varSWDA.def'))

#########################################################
#Read in file containing all constants/directory names ((should be kept in same directory as this script))
#########################################################
fvSWDA=open(currentDir+'\\varSWDA.def','r')
vSWDALines=fvSWDA.readlines()
fvSWDA.close()

outputDir=str(vSWDALines[1].strip())
print('Output directory: '+str(vSWDALines[1]).strip())
if not os.path.isdir(outputDir):
    os.makedirs(outputDir)

#Copy varSWDA.dat file used for this run to the output directory so the user 
#can see what parameters were used to generate the most *recent* output in that output dir
shutil.copy2(currentDir+'\\varSWDA.def',outputDir)

#Make all necessary directories in output directory to save files
if not os.path.isdir(outputDir+'meanMAS/'):
    os.makedirs(outputDir+'meanMAS/')
if not os.path.isdir(outputDir+'prior/'):
    os.makedirs(outputDir+'prior/')
if not os.path.isdir(outputDir+'posterior/'):
    os.makedirs(outputDir+'posterior/')
if not os.path.isdir(outputDir+'ACE/'):
    os.makedirs(outputDir+'ACE/')
if not os.path.isdir(outputDir+'STERA/'):
    os.makedirs(outputDir+'STERA/')
if not os.path.isdir(outputDir+'STERB/'):
    os.makedirs(outputDir+'STERB/')     
if not os.path.isdir(outputDir+'RMSE/'):
    os.makedirs(outputDir+'RMSE/')     

#Make all necessary directories in output directory to save plots
if not os.path.isdir(outputDir+'Plots/'):
    os.makedirs(outputDir+'Plots/')
if not os.path.isdir(outputDir+'Plots/RMSE/'):
    os.makedirs(outputDir+'Plots/RMSE/')
if not os.path.isdir(outputDir+'Plots/speedsOverCR/'):
    os.makedirs(outputDir+'Plots/speedsOverCR/')
if not os.path.isdir(outputDir+'Plots/speedsOverCR/ACE/'):
    os.makedirs(outputDir+'Plots/speedsOverCR/ACE/')
if not os.path.isdir(outputDir+'Plots/speedsOverCR/STERA/'):
    os.makedirs(outputDir+'Plots/speedsOverCR/STERA/')
if not os.path.isdir(outputDir+'Plots/speedsOverCR/STERB/'):
    os.makedirs(outputDir+'Plots/speedsOverCR/STERB/')        

#########################################
#Define constants
#########################################
#Set Initial/Final Carrington Rotations to be run
initdate = vSWDALines[3].strip()
initMJD= dateToMJD(initdate)
noOfWindows= int(vSWDALines[5].strip())#2085
finalMJD=initMJD+(27*noOfWindows)

print('Initial MJD: '+str(initMJD))
print('Final MJD  : '+str(finalMJD))
print('No. Windows: '+str(noOfWindows))

noOfLatPoints=1 #In case we wish to expand prop. model to different lat's

#Longitudinal constants
noOfLonPoints=int(vSWDALines[7].strip())
deltaPhiDeg=float(360.0/noOfLonPoints)
deltaPhi=deltaPhiDeg*(pi/180.0)

#Radial constants (rS=solar radius)
rS=700000 #695700
deltaR=float(vSWDALines[13].strip())*rS
r=arange(float(vSWDALines[9].strip())*rS,float(vSWDALines[11].strip())*rS,deltaR)
noOfRadPoints=len(r)

#Constants for solar wind propagation model
solarRotFreq=(2*pi)/(25.38*24*3600.0)
alpha=float(vSWDALines[37])
rH=float(vSWDALines[39])*rS

#######################################
#Read in which obs. to assimilate
#######################################
#Observations assimilated
#Input a string of A, B and C where A=STERA, B=STERB and C='ACE' eg. for all assim, input obsToUse='ABC'
obsToUse=str(vSWDALines[31].strip()).upper()

obsToUse=''.join(sorted(obsToUse)) #Order the characters alphabetically

print('Observations to be assimilated: '+str(obsToUse))

########################################################
#Initialise arrays for use later in script
########################################################
#Observation vectors
yAllA=9999*ones((noOfWindows*noOfLonPoints))
yAllB=9999*ones((noOfWindows*noOfLonPoints))
yAllC=9999*ones((noOfWindows*noOfLonPoints))

#Vectors to hold solar wind speeds
#vPrior=zeros((noOfWindows,noOfRadPoints,noOfLatPoints,noOfLonPoints))
vBackground=zeros((noOfWindows,noOfRadPoints,noOfLatPoints,noOfLonPoints))
vAnalysis=zeros((noOfWindows,noOfRadPoints,noOfLatPoints,noOfLonPoints))
vMASMean=zeros((noOfWindows,noOfRadPoints,noOfLatPoints,noOfLonPoints))

#If user wishes to plot all windows consecutively, this array holds the solar wind speeds sequentially
vBackAll=zeros((noOfRadPoints,noOfLatPoints,noOfLonPoints*noOfWindows))
#vPriorAll=zeros((noOfRadPoints,noOfLatPoints,noOfLonPoints*noOfWindows))
vAnaAll=zeros((noOfRadPoints,noOfLatPoints,noOfLonPoints*noOfWindows))
vForecastAll=zeros((noOfRadPoints,noOfLatPoints,noOfLonPoints*noOfWindows))

#Initialise arrays to store error information about the solar wind speeds
squareDiffForecast=zeros((noOfWindows))
#squareDiffPrior2=zeros((noOfWindows))
squareDiffAna2=zeros((noOfWindows))
squareDiffMASMean=zeros((noOfWindows))

squareDiffSTERAMASMean=zeros((noOfWindows))
squareDiffSTERAForecast=zeros((noOfWindows))
squareDiffSTERAAna2=zeros((noOfWindows))

squareDiffSTERBMASMean=zeros((noOfWindows))
squareDiffSTERBForecast=zeros((noOfWindows))
squareDiffSTERBAna2=zeros((noOfWindows))

squareDiffACEMASMean=zeros((noOfWindows))
squareDiffACEForecast=zeros((noOfWindows))
squareDiffACEAna2=zeros((noOfWindows))

RMSEForecast=zeros((noOfWindows))
#RMSEPrior2=zeros((noOfWindows))
RMSEAna2=zeros((noOfWindows))
RMSEMASMean=zeros((noOfWindows))

RMSESTERAForecast=zeros((noOfWindows))
#RMSEPrior2=zeros((noOfWindows))
RMSESTERAAna2=zeros((noOfWindows))
RMSESTERAMASMean=zeros((noOfWindows))

RMSESTERBForecast=zeros((noOfWindows))
#RMSEPrior2=zeros((noOfWindows))
RMSESTERBAna2=zeros((noOfWindows))
RMSESTERBMASMean=zeros((noOfWindows))

RMSEACEForecast=zeros((noOfWindows))
#RMSEPrior2=zeros((noOfWindows))
RMSEACEAna2=zeros((noOfWindows))
RMSEACEMASMean=zeros((noOfWindows))

##################################################################################
#Make MJD and initial ensemble files and store them in an array to be read in later
##################################################################################
fileMJD=[] #Variable to hold MJD file names
initMJDCR=[] #Variable to hold CR of initial MJD so we can extract from the correct MAS ens file
MJDLines=[]
CRLines=[]

currMJD=[] #Variables to hold each window's MJD and CR starting points
currCR=[]
midMJD=[] #Variable to hold mid-point MJD of each 27-day period

#Read in table of MJD values and Carrington rotations
MJDstartFile=currentDir+'\\'+str(vSWDALines[17].strip())#+'\\2039CRMJDstart.csv'
fMJDstart=open(MJDstartFile,'r')
fMJDLines=fMJDstart.readlines()
fMJDstart.close()

#Read in MJD values and CR values and store them in an array
for l in range(len(fMJDLines)):
    splitLine=(fMJDLines[l].strip()).split(',')
    
    CRLines.append(int(splitLine[0]))
    MJDLines.append(float(splitLine[1]))

###########################################################################################
#Find MJD starting/mid points for all windows, find which Carrington Rotation they're in 
#and make MJD files for each longitude point in each 27 day period
###########################################################################################
for w in range(noOfWindows):
    #Calculate current window's initial MJD
    currMJD.append(initMJD+(27*w))
    
    #Calculate current window's mid-point MJD
    midMJD.append(currMJD[-1]+13.5)
    
    #Calculate which CR currMJD is in
    #Find index where currMJD first exceeds MJD lines, currMJD will be in the previous line's CR
    reqIndex = next(x for x, val in enumerate(MJDLines) if val > currMJD[-1])
    currCR.append(CRLines[reqIndex-1])
    
    #Store name of required files for creating/reading later
    fileMJD.append(outputDir+'\\MJDfiles\\MJD_'+str(int(currMJD[-1]))+'.dat')
    
    #Check if MJD file exists
    if not os.path.isfile(fileMJD[-1]):
        #If file does not exist, make it
        if not os.path.isdir(outputDir+'\\MJDfiles\\'):
            #Check if directory exists, if not make it
            os.makedirs(outputDir+'\\MJDfiles\\')
        
        makeMJDfile(currMJD[-1],noOfLonPoints,fileMJD[-1],daySolarRot=27)

#####################################################################################
#Check existence of and make initial ensemble files if required for each window
#####################################################################################
fileInitSpecDir=outputDir+'MJDfiles\\MJDstart_'+str(int(currMJD[0]))+'_nWindows_'+str(noOfWindows)+'\\' #Directory to hold solar wind ens files specific to this run
fileInit=[] #Variable to hold required initial solar wind speed ensemble file names for this run
for w in range(noOfWindows):
    fileInit.append(fileInitSpecDir+'vin_ensemble_MJDstart'+str(int(currMJD[w]))+'.dat')
    #Check if  initial ensemble file exists, if not make it
    if not os.path.isfile(fileInit[w]):
        #If file does not exist, make it
        if not os.path.isdir(fileInitSpecDir+'\\'):
            #Check if directory exists, if not make it
            os.makedirs(fileInitSpecDir+'\\')
        
        #Make initial speed ensemble
        makeVinEns27d(currMJD[w],currCR[w],noOfLonPoints,outputDir+'\\MJDfiles\\',currentDir+'\\'+str(vSWDALines[19].strip()),\
                      fileInitSpecDir,daySolarRot=27,noOfMASens=int(vSWDALines[15].strip()))

############################################################################
#Check if mid-points file exists for these windows, if not, make it
############################################################################
mjdCRFile=outputDir+'\\MJDfiles\\MJDMidPoints_'+str(int(currMJD[0]))+'_'+str(int(currMJD[-1]))+'.dat'  
if not os.path.isfile(mjdCRFile):
    #filename=currentDir+'\\'+str(vSWDALines[17].strip())+'\\MJD_'+str(int(initMJD))+'_'+str(int(currMJD[-1]))+'.dat'
    MJDMidFile = open(mjdCRFile, "w")
    for i in range(noOfWindows):
        MJDMidFile.write(str(i)+'\t'+str(midMJD[i])+'\n')
    MJDMidFile.close()

##################################################################################
#Read in mid-points and locations, then extract the longitudes of STEREO A and B that are relevant for this run
##################################################################################    
steraLocFile=currentDir+'\\'+str(vSWDALines[21].strip())
sterbLocFile=currentDir+'\\'+str(vSWDALines[23].strip())

sterALon=readObsLocFileWindows(steraLocFile, mjdCRFile,noOfWindows)
sterBLon=readObsLocFileWindows(sterbLocFile, mjdCRFile,noOfWindows)
print('sterALon='+str(sterALon))
print('sterBLon='+str(sterBLon))    

###############################################################################
#Read in file names of observation files
###############################################################################
fileObsA=currentDir+'\\'+str(vSWDALines[25].strip())
fileObsB=currentDir+'\\'+str(vSWDALines[27].strip())
fileObsACE=currentDir+'\\'+str(vSWDALines[29].strip())

outFile=open(outputDir+'output.txt','w')
for w in range(noOfWindows):
    print('Window No.: '+str(w)+'/'+str(noOfWindows))
    #print('Carr. Rot.: '+str(initCarrRot+w))
    
    #########################################
    #Make mean and B matrix
    #########################################
    #Read in initial ensemble file name and loc rad
    locRad=float(vSWDALines[35].strip())
    
    #Generate prior mean and prior covariance matrix
    meanPrior,B=createMeanAndBMatrix(deltaPhiDeg,locRad,noOfLonPoints,fileInit[w],int(vSWDALines[15].strip()),outputDir)
        
    ################################################
    #Initialise ensemble for this CR
    ################################################
    #Initialise Background and MAS ensemble vectors
    forwardStateBackground=zeros((noOfRadPoints,noOfLatPoints*noOfLonPoints))
    forwardStateMASMean=zeros((noOfRadPoints,noOfLatPoints*noOfLonPoints))
    
    B0InitState=copy(linalg.sqrtm(B))
    
    #Generate a random pertubation vector
    normalVector=zeros((len(forwardStateBackground[0,:])))
    for i in range(len(forwardStateBackground[0,:])):
        normalVector[i]=gauss(0,1)
    
    #Transform perturbation matrix into the distribution we require
    initEns=(transpose(real(B0InitState))).dot(normalVector)
    
    ##############################################
    #Generate prior and MAS Mean state
    ##############################################
    #MAS Mean is just mean of distribution
    forwardStateMASMean[0,:]=copy(meanPrior[:])
    
    #Generate Prior state as a pertubation away from the mean of the prior distribution
    for i in range(len(forwardStateBackground[0,:])):                
        #Prior
        #If you want to shift the prior, do so here in number of lon points shifted
        priorShift=0
        #If perturbation pushes ensemble outside of acceptable values for solar wind speed, adjust accordingly
        if (meanPrior[i-priorShift]+initEns[i]<200):
            forwardStateBackground[0,i]=200+0.1*(abs(200-abs(meanPrior[i-priorShift]+initEns[i])))
        elif (meanPrior[i]+initEns[i]>800):
            forwardStateBackground[0,i]=800-0.1*(abs(800-(meanPrior[i-priorShift]+initEns[i])))
        else:
            forwardStateBackground[0,i]=meanPrior[i-priorShift]+initEns[i]
    
    #Reshape for use in forward model
    vBackground[w,0,:,:]=forwardStateBackground[0,:].reshape(noOfLatPoints,noOfLonPoints)
    vMASMean[w,0,:,:]=forwardStateMASMean[0,:].reshape(noOfLatPoints,noOfLonPoints)
    
    ##########################################################################
    #Run solar wind propagation model to get estimates of the solar wind throughout domain
    ##########################################################################
    for rIndex in range(1,noOfRadPoints):
        #Run prior state forward
        forwardStateBackground[rIndex,:]=forwardRadModel(vBackground[w,rIndex-1,:,:],vBackground[w,0,:,:],r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
        vBackground[w,rIndex,:,:]=copy(forwardStateBackground[rIndex,:].reshape(noOfLatPoints,noOfLonPoints))
        
        #Run MAS Mean forward
        forwardStateMASMean[rIndex,:]=forwardRadModel(vMASMean[w,rIndex-1,:,:],vMASMean[w,0,:,:],r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
        vMASMean[w,rIndex,:,:]=copy(forwardStateMASMean[rIndex,:].reshape(noOfLatPoints,noOfLonPoints))
    
    #Store in prior in sequential matrix   
    vBackAll[:,:,(w*noOfLonPoints):((w+1)*noOfLonPoints)]=copy(vBackground[w,:,:,:])

    
    #############################################################################
    #Generate all observations
    #############################################################################
    #Initialise variables
    #Observation uncertainty at Earth
    noOfObsA=noOfLonPoints
    noOfObsB=noOfLonPoints
    noOfObsC=noOfLonPoints
    
    obsToBeTakenA=range(noOfLonPoints)
    obsToBeTakenB=range(noOfLonPoints)
    obsToBeTakenC=range(noOfLonPoints)
    
    obsNotTakenA=[]
    obsNotTakenB=[]
    obsNotTakenC=[]
    
    yA=zeros((noOfObsA))
    yB=zeros((noOfObsB))
    yC=zeros((noOfObsC))
    
    #STEREO satellite locations
    sterALonCoord=round(sterALon[w]/deltaPhiDeg)
    sterBLonCoord=round(sterBLon[w]/deltaPhiDeg)
    ACELonCoord=0
    #print(sterALonCoord)
    #print(sterBLonCoord)
    
    ###################################
    #Read observation files
    ###################################
    yA=readObsFile(fileObsA,fileMJD[w]) 
    yB=readObsFile(fileObsB,fileMJD[w])
    yC=readObsFile(fileObsACE,fileMJD[w])
    
    yA=copy(yA[::-1])
    yB=copy(yB[::-1])
    yC=copy(yC[::-1])
    
    yAPlot=copy(yA)#[::-1])
    yBPlot=copy(yB)#[::-1])
    yCPlot=copy(yC)#[::-1])
    
    ####################################################################################################################
    #Filter out observations that are unrealistic/were not recorded (negative SW speeds and extremely large SW speeds)
    ####################################################################################################################
    #############
    #STEREO A
    #############
    #Initialise temporary variables
    noOfObsATemp=copy(noOfObsA)
    yATemp=list(copy(yA[:]))
    obsToBeTakenATemp=list(copy(obsToBeTakenA))
    
    #Check if obs. speed negative or greater than 2000km/s remove observation
    for i in range(noOfObsA):
        if ((abs(yA[i])>2000) or (yA[i]<0)):
            #Remove observation from list
            yATemp.remove(yA[i])
            
            #Update plotting variable as NaN
            yAPlot[i]=NaN
            
            #Reduce number of obs. to be taken and record which longitude obs is removed from (and which obs. are to be taken)
            noOfObsATemp=noOfObsATemp-1
            obsToBeTakenATemp.remove(obsToBeTakenA[i])
            obsNotTakenA.append(obsToBeTakenA[i])
    
    ################
    #STEREO B
    ################
    #Initialise temporary variables
    noOfObsBTemp=copy(noOfObsB)
    yBTemp=list(copy(yB[:]))
    obsToBeTakenBTemp=list(copy(obsToBeTakenB))

    #Check if obs. speed negative or greater than 2000km/s remove observation
    for i in range(noOfObsB):
        if ((abs(yB[i])>2000) or (yB[i]<0)):
            #Remove observation from list
            yBTemp.remove(yB[i])
            
            #Update plotting variable as NaN
            yBPlot[i]=NaN
            
            #Reduce number of obs. to be taken and record which longitude obs is removed from (and which obs. are to be taken)
            noOfObsBTemp=noOfObsBTemp-1
            obsToBeTakenBTemp.remove(obsToBeTakenB[i])
            obsNotTakenB.append(obsToBeTakenB[i])
    
    ###############
    #ACE
    ###############
    #Initialise temporary variables
    noOfObsCTemp=copy(noOfObsC)
    yCTemp=list(copy(yC[:]))
    obsToBeTakenCTemp=list(copy(obsToBeTakenC))

    #Check if obs. speed negative or greater than 2000km/s remove observation
    for i in range(noOfObsC):
        if ((abs(yC[i])>2000) or (yC[i]<0)):
            #Remove observation from list
            yCTemp.remove(yC[i])
            
            #Update plotting variable as NaN
            yCPlot[i]=NaN
            
            #Reduce number of obs. to be taken and record which longitude obs is removed from (and which obs. are to be taken)
            noOfObsCTemp=noOfObsCTemp-1
            obsToBeTakenCTemp.remove(obsToBeTakenC[i])     
            obsNotTakenC.append(obsToBeTakenC[i])
            
    ######################################
    #Update obs. vectors
    ######################################
    #Update all observations vector
    yAllA[w*noOfLonPoints:(w+1)*noOfLonPoints]=copy(yAPlot)#[::-1])
    yAllB[w*noOfLonPoints:(w+1)*noOfLonPoints]=copy(yBPlot)#[::-1])
    yAllC[w*noOfLonPoints:(w+1)*noOfLonPoints]=copy(yCPlot)#[::-1])
    
    #Update no. of obs. at each satellite location
    noOfObsA=copy(noOfObsATemp)
    noOfObsB=copy(noOfObsBTemp)
    noOfObsC=copy(noOfObsCTemp)
    
    #Update number of obs. to be taken
    obsToBeTakenA=copy(obsToBeTakenATemp)
    obsToBeTakenB=copy(obsToBeTakenBTemp)
    obsToBeTakenC=copy(obsToBeTakenCTemp)
    
    #Update the observation vector
    yA=copy(yATemp)
    yB=copy(yBTemp)
    yC=copy(yCTemp)
    
    #########################################################
    #Depending on what observations the user has specified to be assimilated
    #Update variables according to the radius at which the STEREO A/B satellites are
    #and the number of observations that occur at each radius
    #########################################################
    if (obsToUse=='A'):
        radObs=zeros((1))
        nRadObs=zeros((len(radObs)))
        
        #Input radius of observations
        radObs[0]=int(215-30-1)
        
        #Input number of observations at each radius
        nRadObs[0]=noOfObsA
        
    elif (obsToUse=='B'):
        radObs=zeros((1))
        nRadObs=zeros((len(radObs)))
        
        #Input radius of observations
        radObs[0]=int(215-30-1)
    
        #Input number of observations at each radius
        nRadObs[0]=noOfObsB
    elif (obsToUse=='C'):
        radObs=zeros((1))
        nRadObs=zeros((len(radObs)))
        
        #Input radius of observations
        radObs[0]=int(213-30-1)
        
        #Input number of observations at each radius
        nRadObs[0]=noOfObsC
    elif (obsToUse=='AB'):
        radObs=zeros((1))
        nRadObs=zeros((len(radObs)))
        
        #Input radius of observations
        radObs[0]=int(215-30-1)
        
        #Input number of observations at each radius
        nRadObs[0]=noOfObsA+noOfObsB
    elif (obsToUse=='AC'):
        radObs=zeros((2))
        nRadObs=zeros((len(radObs)))
        
        #Input radius of observations
        radObs[0]=int(215-30-1)
        radObs[1]=int(213-30-1)
        
        #Input number of observations at each radius
        nRadObs[0]=noOfObsA
        nRadObs[1]=noOfObsA+noOfObsC
    elif (obsToUse=='BC'):
        radObs=zeros((2))
        nRadObs=zeros((len(radObs)))
        
        #Input radius of observations
        radObs[0]=int(215-30-1)
        radObs[1]=int(213-30-1)
        
        #Input number of observations at each radius
        nRadObs[0]=noOfObsB
        nRadObs[1]=noOfObsB+noOfObsC
    elif (obsToUse=='ABC'):
        radObs=zeros((2))
        nRadObs=zeros((len(radObs)))
        
        #Input radius of observations
        radObs[0]=int(215-30-1)
        radObs[1]=int(213-30-1)
        
        #Input number of observations at each radius
        nRadObs[0]=noOfObsA+noOfObsB
        nRadObs[1]=noOfObsA+noOfObsB+noOfObsC
    else:
        print('[31] must equal "A","B","C","AB","AC","BC" or "ABC" (or some permutation of these) without spaces')
        print('where A corresponds to STEREO A, B corresponds to STEREO B and C corresponds to ACE data being assimilated')
        print('Please update [31] accordingly. System will now exit...')
        sys.exit()
        
    noOfObsTotal=int(nRadObs[-1])
    
    #Print number of observations
    print ('len(yA)='+str(len(yA)))
    print ('len(yB)='+str(len(yB)))
    print ('len(yC)='+str(len(yC)))
    print('Total number of observations: '+str(noOfObsTotal))
    
    #####################################
    #Initialise observation variables
    #####################################
    #Initialise obs.
    y=zeros((noOfObsTotal))
    
    #Initialise observation operator variables
    HA=zeros((noOfObsA,noOfLonPoints))
    HB=zeros((noOfObsB,noOfLonPoints))
    HC=zeros((noOfObsC,noOfLonPoints))
    
    #Initialise observation error covariance matrix variables
    obsUncA=zeros((noOfObsA))
    obsUncB=zeros((noOfObsB))
    obsUncC=zeros((noOfObsC))
    RA=zeros((noOfObsA,noOfObsA))   #Obs. error covar. for STERA
    RB=zeros((noOfObsB,noOfObsB))   #Obs. error covar. for STERB
    RC=zeros((noOfObsC,noOfObsC))   #Obs. error covar. for ACE
    
    R=zeros((noOfObsTotal,noOfObsTotal))    #Obs. error covar. for all obs. to be assimilated
    
    ###################################################        
    #Generate observation operators
    ###################################################    
    H=zeros((noOfObsTotal,noOfLonPoints))
    HA=obsOp(noOfLonPoints,obsToBeTakenA,[sterALonCoord])
    HB=obsOp(noOfLonPoints,obsToBeTakenB,[sterBLonCoord])
    HC=obsOp(noOfLonPoints,obsToBeTakenC,[0])
    
    ######################################################################################################
    #Make observation error covariance matrices (assumed diagonal, i.e. all observations are uncorrelated)
    ######################################################################################################
    #Read in whether a constant obs. covariance matrix is required or if 
    #the uncertainty should be a percentage of the prior solar wind speed at the observation radius
    splitLine=(vSWDALines[33].strip()).split(' ')
    print(splitLine)
    
    if str(splitLine[0])=='B':
        #Generate observation errors proportional to mean of solar wind speed at obs. radius
        obsUncA=(float(splitLine[1])*vBackground[w,int(radObs[0]),0,:].mean())*ones(noOfObsA)
        obsUncB=(float(splitLine[2])*vBackground[w,int(radObs[0]),0,:].mean())*ones(noOfObsB)
        if len(radObs)==2:
            obsUncC=(float(splitLine[3])*vBackground[w,int(radObs[1]),0,:].mean())*ones(noOfObsC)
        else:
            obsUncC=(float(splitLine[3])*vBackground[w,int(radObs[0]),0,:].mean())*ones(noOfObsC)
    elif str(splitLine[0])=='C':
        #Generate observation errors as a constant supplied by the user
        obsUncA=float(splitLine[1])*ones(noOfObsA)
        obsUncB=float(splitLine[2])*ones(noOfObsB)
        obsUncC=float(splitLine[3])*ones(noOfObsC)
    else:
        print('First character should be a "B" or "C"')
        print('where B corresponds to a observation error standard deviation proportional to mean prior solar wind at obs. radius')
        print('and C corresponds to a constant observation error standard deviation being used.')
        print('Please update [33] accordingly. System will now exit...')
        sys.exit()
    
    #Assume observations at different satellites are not correlated
    for i in range(noOfObsA):
        RA[i,i]=obsUncA[i]*obsUncA[i]
    for i in range(noOfObsB):
        RB[i,i]=obsUncB[i]*obsUncB[i]
    for i in range(noOfObsC):
        RC[i,i]=obsUncC[i]*obsUncC[i]
    
    '''
    ################################################
    #Script to read observation error covariance from a file (code that may be used in a future update)
    ################################################
    #Observation error covariances read from file
    RfromFile=zeros((3*noOfLonPoints,3*noOfLonPoints))
    
    rFile=open('<directory with R matrix in>\\<file with R matrix in>.txt','r')
    RfromFile=loadtxt(rFile)
    rFile.close()
    
    RfromFileTemp=RfromFile.copy()
    obsNotTakenTotal=concatenate((obsNotTakenA, obsNotTakenB+noOfLonPoints*ones(len(obsNotTakenB)),obsNotTakenC+(2*noOfLonPoints*ones(len(obsNotTakenC)))), axis=None)
    #obsNotTakenTotal.append()
    #obsNotTakenTotal.append()
    
    print(obsNotTakenTotal)
    for i in sorted(obsNotTakenTotal,reverse=True):
        RfromFileTemp=delete(RfromFile,i,0)
        RfromFileTemp=delete(RfromFile,i,1)
    R=RfromFileTemp.copy()
    
    #sys.exit()
    '''
    #########################################################################################
    #Update the generic DA variables depending upon what obs. the user wishes to assimilate
    #########################################################################################
    #Input appropriate parts into y, R and H to be used in DA
    if (obsToUse=='A'):
        #Make Observations
        y[:noOfObsA]=copy(yA)
                
        #Input appropriate components into full observation covariance matrix
        R[:noOfObsA,:noOfObsA]=copy(RA)

        #Input appropriate components into full observation covariance matrix
        H[:noOfObsA,:]=copy(HA)
        
    elif (obsToUse=='B'):
        #Make Observations
        y[:noOfObsB]=copy(yB)
                
        #Input appropriate components into full observation covariance matrix
        R[:noOfObsB,:noOfObsB]=copy(RB)
        
        #Input appropriate components into full observation covariance matrix
        H[:noOfObsB,:]=copy(HB)
    elif (obsToUse=='C'):
        #Make Observations
        y[:noOfObsC]=copy(yC)
        
        #Input appropriate components into full observation covariance matrix
        R[:noOfObsC,:noOfObsC]=copy(RC)
                
        #Input appropriate components into full observation covariance matrix
        H[:noOfObsC,:]=copy(HC)
    elif (obsToUse=='AB'):
        #Make Observations
        y[:noOfObsA]=copy(yA)
        y[noOfObsA:(noOfObsA+noOfObsB)]=copy(yB)
        
        #Input appropriate components into full observation covariance matrix
        R[:noOfObsA,:noOfObsA]=copy(RA)
        R[noOfObsA:(noOfObsA+noOfObsB),noOfObsA:(noOfObsA+noOfObsB)]=copy(RB)
        
        
        #Input appropriate components into full observation covariance matrix
        H[:noOfObsA,:]=copy(HA)
        H[noOfObsA:(noOfObsA+noOfObsB),:]=copy(HB)
    elif (obsToUse=='AC'):
        #Make Observations
        y[:noOfObsA]=copy(yA)
        y[noOfObsA:]=copy(yC)
                
        #Input appropriate components into full observation covariance matrix
        R[:noOfObsA,:noOfObsA]=copy(RA)
        R[noOfObsA:,noOfObsA:]=copy(RC)
                
        #Input appropriate components into full observation covariance matrix
        H[:noOfObsA,:]=copy(HA)
        H[noOfObsA:,:]=copy(HC)
    elif (obsToUse=='BC'):
        #Make Observations
        y[:noOfObsB]=copy(yB)
        y[noOfObsB:]=copy(yC)
        
        
        #Input appropriate components into full observation covariance matrix
        R[:noOfObsB,:noOfObsB]=copy(RB)
        R[noOfObsB:,noOfObsB:]=copy(RC)
        
        #Input appropriate components into full observation covariance matrix
        H[:noOfObsB,:]=copy(HB)
        H[noOfObsB:,:]=copy(HC)
    elif (obsToUse=='ABC'):
        #Make Observations
        y[:noOfObsA]=copy(yA)
        y[noOfObsA:(noOfObsA+noOfObsB)]=copy(yB)
        y[(noOfObsA+noOfObsB):]=copy(yC)
        
        #Input appropriate components into full observation covariance matrix
        R[:noOfObsA,:noOfObsA]=copy(RA)
        R[noOfObsA:(noOfObsA+noOfObsB),noOfObsA:(noOfObsA+noOfObsB)]=copy(RB)
        R[(noOfObsA+noOfObsB):,(noOfObsA+noOfObsB):]=copy(RC)
       
        #Input appropriate components into full observation covariance matrix
        H[:noOfObsA,:]=copy(HA)
        H[noOfObsA:(noOfObsA+noOfObsB),:]=copy(HB)
        H[(noOfObsA+noOfObsB):,:]=copy(HC)
    else:
        print('obsToUse must equal "A","B","C","AB","AC","BC" or "ABC" (or some permutation of these) without spaces')
        print('where A corresponds to STEREO A, B corresponds to STEREO B and C corresponds to ACE data being assimilated')
        print('Please update obsToUse accordingly. System will now exit...')
        sys.exit()
    
    #############################################################################
    #Data Assimilation
    #############################################################################
    #Specify number of iterations to be used per window (only useful if running inner and outer loops in minimisation algorithm)
    #If not performing inner/outer loops for minimisation of cost function, noOfIter=1    
    noOfIter=1
    
    for i in range(noOfIter):
        if i==0:
            #Initialise variables in initial iteration
            vIter=zeros((noOfIter+1,noOfRadPoints,noOfLatPoints,noOfLonPoints))
            forwardStateIter=zeros((noOfIter+1,noOfRadPoints,noOfLatPoints*noOfLonPoints))
            costFuncVal=zeros((noOfIter+1))
            
            #Set the initial solar wind speed as equal to the prior solar wind speed
            vIter[0,0,:,:]=copy(vBackground[w,0,:,:])
            forwardStateIter[0,0,:]=vIter[0,0,:,:].reshape(noOfLatPoints*noOfLonPoints)
            
            #Run initial solar wind speed out into the heliosphere (from 30rS -> 215rS)
            for rIndex in range(1,noOfRadPoints):
                #forwardStateIter[0,rIndex,:]=forwardRadModel(vIter[0,rIndex-1,:,:],v0,r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
                forwardStateIter[0,rIndex,:]=forwardRadModel(vIter[0,rIndex-1,:,:],vIter[0,0,:,:],r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
                vIter[0,rIndex,:,:]=copy(forwardStateIter[0,rIndex,:].reshape(noOfLatPoints,noOfLonPoints))
            
            #Initialise state vector variables at inner radius
            x=forwardStateIter[0,0,:]
            xb=forwardStateIter[0,0,:]
            
            #Calculate cost function at initial iteration
            costFuncVal[i]=calcCostFunc(B,R,H,forwardStateIter[0,0,:],xb,vIter[0,:,:,:],y,radObs,nRadObs,noOfObsTotal,noOfLatPoints,noOfLonPoints)
        else:
            #Update background state vector to previous iteration's analysis state
            xb=forwardStateIter[i,0,:]
        
        ###########################################################################################
        #Minimise the cost function (change minimisation methodology here)
        ###########################################################################################
        #Perform minimisation of cost function to obtain analysis state
        resOpt=optimize.minimize(fun=calcCostFuncForCG, x0=xb, args=(B,R,H,xb,vIter[i,:,:,:],y,radObs,nRadObs,noOfObsTotal,\
                                                                     r,rH,deltaR,deltaPhi,alpha,solarRotFreq,noOfRadPoints,noOfLatPoints,noOfLonPoints,i),\
                                method='BFGS',jac=makeGradCG,options={'gtol': 1e-5, 'disp': True})#,'maxiter':20})
        
        ###########################################################################
        #Extract the  analysis solar wind speed in both speed matrix and state vector form
        ###########################################################################
        vIter[i+1,0,:,:]=copy((resOpt.x).reshape(noOfLatPoints,noOfLonPoints))
        forwardStateIter[i+1,0,:]=copy(resOpt.x)
        
        #####################################################################
        #Run analysed DA state out into the heliosphere using the numerical model
        #####################################################################
        for rIndex in range(1,noOfRadPoints):
            #print('rIndex='+str(rIndex))
            forwardStateIter[i+1,rIndex,:]=forwardRadModel(vIter[i+1,rIndex-1,:,:],vIter[i+1,0,:,:],r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
            vIter[i+1,rIndex,:,:]=copy(forwardStateIter[i+1,rIndex,:].reshape(noOfLatPoints,noOfLonPoints))
        
        print ('After Analysis')
        ###################################################################################
        #Calculate cost function after DA analysis and store in cost function variable
        ###################################################################################
        costFuncVal[i+1]=calcCostFunc(B,R,H,forwardStateIter[i+1,0,:],xb,vIter[i+1,:,:,:],y,radObs,nRadObs,noOfObsTotal,noOfLatPoints,noOfLonPoints)
    
    ################################################
    #Extract variables into larger arrays 
    ################################################
    #Store analysis for current window in file        
    vAnalysis[w,:,:,:]=vIter[-1,:,:,:]
    vForecastAll[:,:,(w*noOfLonPoints):((w+1)*noOfLonPoints)]=vIter[0,:,:,:]
    vAnaAll[:,:,(w*noOfLonPoints):((w+1)*noOfLonPoints)]=vIter[-1,:,:,:]
    
    #############################################################
    #Extract variables for plotting at each observation location
    #############################################################
    #Initialise variables
    lengthVPlot=len(arange(-90,449,deltaPhiDeg))
    vPlot=zeros((3,noOfLonPoints))
    vPlotTemp=zeros((3,noOfLonPoints))
    
    #Initialise plotting variables for each obs. location
    #Column 0=Obs. values, column 1=MAS Mean speed, column 2=Prior speed, column 3= Posterior speed
    vACEPlot=ones((4,noOfLonPoints))*NaN
    vSterAPlot=ones((4,noOfLonPoints))*NaN
    vSterBPlot=ones((4,noOfLonPoints))*NaN
    
    ######################################################
    #Extract variables for ACE satellite and order variables as ascending in time (-Carr. Rot.)
    ######################################################
    #Extract velocities at ACE radius (and reorder)
    vPlot[0,:]=copy(vMASMean[w,213-30-1,0,::-1])
    vPlot[1,:]=copy(vIter[0,213-31,0,::-1])
    vPlot[2,:]=copy(vIter[-1,213-31,0,::-1])
    
    #Extract velocities at STERA/STERB radii (and reorder)
    vPlotTemp[0,:]=copy(vMASMean[w,215-30-1,0,::-1])
    vPlotTemp[1,:]=copy(vIter[0,215-31,0,::-1])
    vPlotTemp[2,:]=copy(vIter[-1,215-31,0,::-1])
    
    ######################################################
    #Copy observation data (and reorder)
    ######################################################
    vSterAPlot[0,:]=copy(yAPlot[::-1])
    vSterBPlot[0,:]=copy(yBPlot[::-1])
    vACEPlot[0,:]=copy(yCPlot[::-1])
    
    
    for i in range(noOfLonPoints):        
        #Extract relevant solar wind velocities at STEREO A location
        vSterAPlot[1,int(mod(i+sterALonCoord,noOfLonPoints))]=copy(vPlotTemp[0,i])
        vSterAPlot[2,int(mod(i+sterALonCoord,noOfLonPoints))]=copy(vPlotTemp[1,i])
        vSterAPlot[3,int(mod(i+sterALonCoord,noOfLonPoints))]=copy(vPlotTemp[2,i])
        
        #Extract relevant solar wind velocities at STEREO B location
        vSterBPlot[1,int(mod(i+sterBLonCoord,noOfLonPoints))]=copy(vPlotTemp[0,i])
        vSterBPlot[2,int(mod(i+sterBLonCoord,noOfLonPoints))]=copy(vPlotTemp[1,i])
        vSterBPlot[3,int(mod(i+sterBLonCoord,noOfLonPoints))]=copy(vPlotTemp[2,i])
        
        #Extract relevant solar wind velocities at ACE location
        vACEPlot[1,int(mod(i,noOfLonPoints))]=copy(vPlot[0,i])
        vACEPlot[2,int(mod(i,noOfLonPoints))]=copy(vPlot[1,i])
        vACEPlot[3,int(mod(i,noOfLonPoints))]=copy(vPlot[2,i])
    
    ########################################################################
    #Calculate RMSEs of Prior, forecast and Posterior compared to obs. data
    ########################################################################
    #Variables for calculating number of obs. at each obs. location (should be same as noOfObsA, noOfObsB and noOfObsC)
    noOfSTERAobs2=0
    noOfSTERBobs2=0
    noOfACEobs2=0
    
    #Calculate the sum of the squared errors at each obs. location
    for lo in range(noOfLonPoints):
        #Check if STERA observation is non-nan at current longitude
        if isnan(vSterAPlot[0,lo])==0:
            #Update number of non-nan obs at STERA
            noOfSTERAobs2=noOfSTERAobs2+1
            
            #Calculate at square diff. of STERA obs and calculated velocities
            squareDiffSTERAMASMean[w]=squareDiffSTERAMASMean[w]+(vSterAPlot[1,lo]-vSterAPlot[0,lo])**2
            squareDiffSTERAForecast[w]=squareDiffSTERAForecast[w]+(vSterAPlot[2,lo]-vSterAPlot[0,lo])**2
            squareDiffSTERAAna2[w]=squareDiffSTERAAna2[w]+(vSterAPlot[3,lo]-vSterAPlot[0,lo])**2
        
        #Check if STERB observation is non-nan at current longitude
        if isnan(vSterBPlot[0,lo])==0:
            #Update number of non-nan obs at STERB
            noOfSTERBobs2=noOfSTERBobs2+1
            
            #Calculate at square diff. of STERB obs and calculated velocities
            squareDiffSTERBMASMean[w]=squareDiffSTERBMASMean[w]+(vSterBPlot[1,lo]-vSterBPlot[0,lo])**2
            squareDiffSTERBForecast[w]=squareDiffSTERBForecast[w]+(vSterBPlot[2,lo]-vSterBPlot[0,lo])**2
            squareDiffSTERBAna2[w]=squareDiffSTERBAna2[w]+(vSterBPlot[3,lo]-vSterBPlot[0,lo])**2
        
        #Check if ACE observation is non-nan at current longitude
        if isnan(vACEPlot[0,lo])==0:
            #Update number of non-nan obs at ACE
            noOfACEobs2=noOfACEobs2+1
            
            #Calculate at square diff. of ACE obs and calculated velocities
            squareDiffACEMASMean[w]=squareDiffACEMASMean[w]+(vACEPlot[1,lo]-vACEPlot[0,lo])**2
            squareDiffACEForecast[w]=squareDiffACEForecast[w]+(vACEPlot[2,lo]-vACEPlot[0,lo])**2
            squareDiffACEAna2[w]=squareDiffACEAna2[w]+(vACEPlot[3,lo]-vACEPlot[0,lo])**2
    
    #Take sqrt and divide by number of obs at each location to get RMSEs
    #STEREO A
    RMSESTERAMASMean[w]=sqrt(squareDiffSTERAMASMean[w]/noOfSTERAobs2)
    RMSESTERAForecast[w]=sqrt(squareDiffSTERAForecast[w]/noOfSTERAobs2)
    RMSESTERAAna2[w]=sqrt(squareDiffSTERAAna2[w]/noOfSTERAobs2)
    
    #STEREO B
    RMSESTERBMASMean[w]=sqrt(squareDiffSTERBMASMean[w]/noOfSTERBobs2)
    RMSESTERBForecast[w]=sqrt(squareDiffSTERBForecast[w]/noOfSTERBobs2)
    RMSESTERBAna2[w]=sqrt(squareDiffSTERBAna2[w]/noOfSTERBobs2)
    
    #ACE
    RMSEACEMASMean[w]=sqrt(squareDiffACEMASMean[w]/noOfACEobs2)
    RMSEACEForecast[w]=sqrt(squareDiffACEForecast[w]/noOfACEobs2)
    RMSEACEAna2[w]=sqrt(squareDiffACEAna2[w]/noOfACEobs2)
    
    
    ##############################################################################
    #Plot solar wind speeds at STEREO A, STEREO B and ACE
    ##############################################################################
    figure(8000+w)
    plot(arange(0,359,deltaPhiDeg),vSterAPlot[0,:],color='k',linewidth=3.0,label='STEREO A Data')
    plot(arange(0,359,deltaPhiDeg),vSterAPlot[1,:],color='m',linewidth=3.0,label='MAS ens. mean')
    plot(arange(0,359,deltaPhiDeg),vSterAPlot[2,:],color='b',linewidth=3.0,label='Prior')
    plot(arange(0,359,deltaPhiDeg),vSterAPlot[3,:],color='g',linewidth=3.0,label='Posterior')
    xlabel('Time (days)',fontsize=18)
    ylabel('Speed (km/s)',fontsize=18)
    xlim(0,360)
    ylim(240,800)
    yticks(arange(300,801,100), ('300','400','500','600','700','800'),fontsize=18)
    #xticks(arange(0,361,120), ('0', '9', '18', '27'),fontsize=18)
    xticks(append(arange(0,361,53.3),360), ('0','4','8','12','16','20','24','27'),fontsize=18)
    #if w==0:
    title(r'Solar wind speed at STEREO A',fontsize=18)
    legend()
    tight_layout()
    savefig(outputDir+'Plots/speedsOverCR/STERA/SWspeedSTEREOA_MJDstart'+str(int(currMJD[w]))+'.png')
    #show()
    
    figure(9000+w)
    plot(arange(0,359,deltaPhiDeg),vSterBPlot[0,:],color='k',linewidth=3.0,label='STEREO B Data')
    plot(arange(0,359,deltaPhiDeg),vSterBPlot[1,:],color='m',linewidth=3.0,label='MAS ens. mean')
    plot(arange(0,359,deltaPhiDeg),vSterBPlot[2,:],color='b',linewidth=3.0,label='Prior')
    plot(arange(0,359,deltaPhiDeg),vSterBPlot[3,:],color='g',linewidth=3.0,label='Posterior')
    xlabel('Time (days)',fontsize=18)
    ylabel('Speed (km/s)',fontsize=18)
    xlim(0,360)
    ylim(240,800)
    yticks(arange(300,801,100),('300','400','500','600','700','800'),fontsize=18)
    #xticks(arange(0,361,120), ('0', '9', '18', '27'),fontsize=18)
    xticks(append(arange(0,361,53.3),360), ('0','4','8','12','16','20','24','27'),fontsize=18)
    #if w==0:
    title(r'Solar wind speed at STEREO B',fontsize=18)
    legend()
    tight_layout()
    savefig(outputDir+'Plots/speedsOverCR/STERB/SWspeedSTEREOB_MJDstart'+str(int(currMJD[w]))+'.png')
    #show()
    
    figure(9500+w)
    plot(arange(0,359,deltaPhiDeg),vACEPlot[0,:],color='k',linewidth=3.0,label='ACE Data')
    plot(arange(0,359,deltaPhiDeg),vACEPlot[1,:],color='m',linewidth=3.0,label='MAS ens. mean')
    plot(arange(0,359,deltaPhiDeg),vACEPlot[2,:],color='b',linewidth=3.0,label='Prior')
    plot(arange(0,359,deltaPhiDeg),vACEPlot[3,:],color='g',linewidth=3.0,label='Posterior')
    xlabel('Time (days)',fontsize=18)
    ylabel('Speed (km/s)',fontsize=18)
    xlim(0,360)
    ylim(240,800)
    yticks(arange(300,801,100),('300','400','500','600','700','800'),fontsize=18)
    #xticks(arange(0,361,120), ('0', '9', '18', '27'),fontsize=18)
    xticks(append(arange(0,361,53.3),360), ('0','4','8','12','16','20','24','27'),fontsize=18)
    #if w==0:
    title(r'Solar wind speed at ACE',fontsize=18)
    legend()
    tight_layout()
    savefig(outputDir+'Plots/speedsOverCR/ACE/SWspeedACE_MJDstart'+str(int(currMJD[w]))+'.png')
    #show()
    
    ###########################################################################
    #Output RMSEs for user
    ###########################################################################
    print('------------------------------------------------------')
    print('Window no. '+str(w))
    print(' ')
    #print('RMSE Prior    = '+str(RMSEPrior2[w]))
    print('RMSE STEREO A MASMean = '+str(RMSESTERAMASMean[w]))
    print('RMSE STEREO A Forecast = '+str(RMSESTERAForecast[w]))
    print('RMSE STEREO A Posterior= '+str(RMSESTERAAna2[w]))
    print(' ')
    #print ('RMSE Prior    = '+str(RMSEPrior2[w]))
    print('RMSE STEREO B MASMean = '+str(RMSESTERBMASMean[w]))
    print('RMSE STEREO B Forecast = '+str(RMSESTERBForecast[w]))
    print('RMSE STEREO B Posterior= '+str(RMSESTERBAna2[w]))
    print(' ')
    #print('RMSE Prior    = '+str(RMSEPrior2[w]))
    print('RMSE ACE MASMean = '+str(RMSEACEMASMean[w]))
    print('RMSE ACE Forecast = '+str(RMSEACEForecast[w]))
    print('RMSE ACE Posterior= '+str(RMSEACEAna2[w]))
    
    ###########################################################################
    #Write RMSEs into output file for further use
    ###########################################################################
    outFile.write('------------------------------------------------------\n')
    outFile.write('Window no. '+str(w)+'\n\n')
    #print('RMSE Prior    = '+str(RMSEPrior2[w]))
    outFile.write('RMSE STEREO A MASMean = '+str(RMSESTERAMASMean[w])+'\n')
    outFile.write('RMSE STEREO A Forecast = '+str(RMSESTERAForecast[w])+'\n')
    outFile.write('RMSE STEREO A Posterior= '+str(RMSESTERAAna2[w])+'\n\n')
    #print(' ')
    #print('RMSE Prior    = '+str(RMSEPrior2[w]))
    outFile.write('RMSE STEREO B MASMean = '+str(RMSESTERBMASMean[w])+'\n')
    outFile.write('RMSE STEREO B Forecast = '+str(RMSESTERBForecast[w])+'\n')
    outFile.write('RMSE STEREO B Posterior= '+str(RMSESTERBAna2[w])+'\n\n')
    #print(' ')
    #print('RMSE Prior    = '+str(RMSEPrior2[w]))
    outFile.write('RMSE ACE MASMean = '+str(RMSEACEMASMean[w])+'\n')
    outFile.write('RMSE ACE Forecast = '+str(RMSEACEForecast[w])+'\n')
    outFile.write('RMSE ACE Posterior= '+str(RMSEACEAna2[w])+'\n')
    outFile.write('------------------------------------------------------\n\n')   
    
    #############################################################################
    #Write MASMean, prior and posterior arrays to file for later use
    #############################################################################
    outWindowMASMeanFile=open(outputDir+'meanMAS/meanMAS_MJDstart'+str(int(currMJD[w]))+'.txt','w')
    outWindowPriorFile=open(outputDir+'prior/prior_MJDstart'+str(int(currMJD[w]))+'.txt','w')
    outWindowPostFile=open(outputDir+'posterior/post_MJDstart'+str(int(currMJD[w]))+'.txt','w')
    
    #Write ACE, STEREOA and STEREOB values during solar rotation to file for later use
    outWindowACEFile=open(outputDir+'ACE/ACE_MJDstart'+str(int(currMJD[w]))+'.txt','w')
    outWindowSTERAFile=open(outputDir+'STERA/STERA_MJDstart'+str(int(currMJD[w]))+'.txt','w')
    outWindowSTERBFile=open(outputDir+'STERB/STERB_MJDstart'+str(int(currMJD[w]))+'.txt','w')
    #outWindowSDPriorFile=open(outputDir+'sdPrior_CR'+str(initCarrRot+w)+'.txt','w')
    #outWindowSDPostFile=open(outputDir+'sdPost_CR'+str(initCarrRot+w)+'.txt','w')
    
    savetxt(outWindowMASMeanFile,vMASMean[w,:,0,::-1])
    savetxt(outWindowPriorFile,vBackground[w,:,0,::-1])
    savetxt(outWindowPostFile,vAnalysis[w,:,0,::-1])
    
    savetxt(outWindowACEFile,yCPlot[::-1])
    savetxt(outWindowSTERAFile,yAPlot[::-1])
    savetxt(outWindowSTERBFile,yBPlot[::-1])
    
    outWindowMASMeanFile.close()
    outWindowPriorFile.close()
    outWindowPostFile.close()
    outWindowACEFile.close()
    outWindowSTERAFile.close()
    outWindowSTERBFile.close()
    #sys.exit()

outFile.close()

####################################
#Plot RMSEs
####################################
figure(91000)
plot(midMJD,RMSEACEMASMean,color='m',marker='x',markersize=12,label='MASMean RMSE')
plot(midMJD,RMSEACEForecast,color='b',marker='x',markersize=12,label='Prior RMSE')
plot(midMJD,RMSEACEAna2,color='g',marker='x',markersize=12,label='Posterior RMSE')
xlabel('MJD',fontsize=18)
ylabel('RMSE (km/s)',fontsize=18)
title(r'RMSE at $213r_S$',fontsize=18)
xlim(currMJD[0],currMJD[-1]+27)
legend()
tight_layout()
savefig(outputDir+'Plots/RMSE/RMSE_ACE.png')    

figure(92000)
plot(midMJD,RMSESTERAMASMean,color='m',marker='x',markersize=12,label='MASMean RMSE')
plot(midMJD,RMSESTERAForecast,color='b',marker='x',markersize=12,label='Prior RMSE')
plot(midMJD,RMSESTERAAna2,color='g',marker='x',markersize=12,label='Posterior RMSE')
xlabel('MJD',fontsize=18)
ylabel('RMSE (km/s)',fontsize=18)
title(r'RMSE at STEREO A',fontsize=18)
xlim(currMJD[0],currMJD[-1]+27)
legend()
tight_layout()
savefig(outputDir+'Plots/RMSE/RMSE_STERA.png') 

figure(93000)
plot(midMJD,RMSESTERBMASMean,color='m',marker='x',markersize=12,label='MASMean RMSE')
plot(midMJD,RMSESTERBForecast,color='b',marker='x',markersize=12,label='Prior RMSE')
plot(midMJD,RMSESTERBAna2,color='g',marker='x',markersize=12,label='Posterior RMSE')
xlabel('MJD',fontsize=18)
ylabel('RMSE (km/s)',fontsize=18)
title(r'RMSE at STEREO B',fontsize=18)
xlim(currMJD[0],currMJD[-1]+27)
legend()
tight_layout()
savefig(outputDir+'Plots/RMSE/RMSE_STERB.png') 

#############################################
#Write RMSEs to file for later use
#############################################
outWindowRMSEMASMeanACEFile=open(outputDir+'RMSE/RMSE_MASMeanACE_MJD'+str(int(currMJD[0]))+'_'+str(int(currMJD[-1]+27))+'.txt','w')
outWindowRMSEPriorACEFile=open(outputDir+'RMSE/RMSE_priorACE_MJD'+str(int(currMJD[0]))+'_'+str(int(currMJD[-1]+27))+'.txt','w')
outWindowRMSEPostACEFile=open(outputDir+'RMSE/RMSE_postACE_MJD'+str(int(currMJD[0]))+'_'+str(int(currMJD[-1]+27))+'.txt','w')

outWindowRMSEMASMeanSTERAFile=open(outputDir+'RMSE/RMSE_MASMeanSTERA_MJD'+str(int(currMJD[0]))+'_'+str(int(currMJD[-1]+27))+'.txt','w')
outWindowRMSEPriorSTERAFile=open(outputDir+'RMSE/RMSE_priorSTERA_MJD'+str(int(currMJD[0]))+'_'+str(int(currMJD[-1]+27))+'.txt','w')
outWindowRMSEPostSTERAFile=open(outputDir+'RMSE/RMSE_postSTERA_MJD'+str(int(currMJD[0]))+'_'+str(int(currMJD[-1]+27))+'.txt','w')

outWindowRMSEMASMeanSTERBFile=open(outputDir+'RMSE/RMSE_MASMeanSTERB_MJD'+str(int(currMJD[0]))+'_'+str(int(currMJD[-1]+27))+'.txt','w')
outWindowRMSEPriorSTERBFile=open(outputDir+'RMSE/RMSE_priorSTERB_MJD'+str(int(currMJD[0]))+'_'+str(int(currMJD[-1]+27))+'.txt','w')
outWindowRMSEPostSTERBFile=open(outputDir+'RMSE/RMSE_postSTERB_MJD'+str(int(currMJD[0]))+'_'+str(int(currMJD[-1]+27))+'.txt','w')

savetxt(outWindowRMSEMASMeanACEFile,RMSEACEMASMean)
savetxt(outWindowRMSEPriorACEFile,RMSEACEForecast)
savetxt(outWindowRMSEPostACEFile,RMSEACEAna2)

savetxt(outWindowRMSEMASMeanSTERAFile,RMSESTERAMASMean)
savetxt(outWindowRMSEPriorSTERAFile,RMSESTERAForecast)
savetxt(outWindowRMSEPostSTERAFile,RMSESTERAAna2)

savetxt(outWindowRMSEMASMeanSTERBFile,RMSESTERBMASMean)
savetxt(outWindowRMSEPriorSTERBFile,RMSESTERBForecast)
savetxt(outWindowRMSEPostSTERBFile,RMSESTERBAna2)    

outWindowRMSEMASMeanACEFile.close()
outWindowRMSEPriorACEFile.close()
outWindowRMSEPostACEFile.close()

outWindowRMSEMASMeanSTERAFile.close()
outWindowRMSEPriorSTERAFile.close()
outWindowRMSEPostSTERAFile.close()

outWindowRMSEMASMeanSTERBFile.close()
outWindowRMSEPriorSTERBFile.close()
outWindowRMSEPostSTERBFile.close()