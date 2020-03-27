# -*- coding: utf-8 -*-
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
import datetime as dt
#from jdcal import gcal2jd#, jd2gca
#import julian

def initBoundCond(initFile):
    #Read input ensemble file
    fVin=open(initFile,'r')
    vInLines=fVin.readlines()
    fVin.close()

    #print len(vInLines)
    vIn=zeros((576,len(vInLines)))
    for i in range(len(vInLines)):
        #print i
        #print vInLines[i]
        vIn[:,i]=(vInLines[i]).split('   ')
    
    return vIn[0,:]
    
def readObsFile(obsFile, mjdFile):
    #Read STEREOobs file
    fObs=open(obsFile,'r')
    vObsLines=fObs.readlines()
    fObs.close()
    
    MJDfiles=open(mjdFile)
    MJDvalues=MJDfiles.readlines()
    MJDfiles.close()
    
    MJDvalues.reverse()
    
    obsAll=[]
    obsMJDtime=[]
    #print len(vObsLines)
    for ln in vObsLines:
        s=ln.split()
        #print(s)
        MJDtime=gregToMJD(int(s[0].strip()),int(s[1].strip()),int(s[2].strip()))
        #print(MJDtime)
        
        obsMJDtime.append(MJDtime)
        obsAll.append(float(s[3].strip()))
    
    #print obsMJDtime[-25:-1]
    observations=[]
    i=0
    #j=0
    obsNo=0    
    for lnMJD in MJDvalues:
        #i=j
        while ((obsMJDtime[i]<float(lnMJD.strip())) and (i<(len(obsAll)+2))):
            if i==(len(obsAll)+1):
                print('Date out of bounds')
                sys.exit()
                
            i=i+1
        #print i
        observations.append(obsAll[i])
        #obsNo=obsNo+1
        
    return observations

def readObsFileCom(obsFile, mjdFile):
    #Read STEREOobs file
    fObs=open(obsFile,'r')
    vObsLines=fObs.readlines()
    fObs.close()
    
    MJDfiles=open(mjdFile)
    MJDvalues=MJDfiles.readlines()
    MJDfiles.close()
    
    MJDvalues.reverse()
    
    obsAll=[]
    obsMJDtime=[]
    #print len(vObsLines)
    for ln in vObsLines:
        s=ln.split(',')
        #print(s)
        MJDtime=gregToMJD(int(s[0].strip()),int(s[1].strip()),int(s[2].strip()))
        #print(MJDtime)
        
        obsMJDtime.append(MJDtime)
        obsAll.append(float(s[3].strip()))
    
    #print obsMJDtime[-25:-1]
    observations=[]
    i=0
    #j=0
    obsNo=0    
    for lnMJD in MJDvalues:
        #i=j
        while ((obsMJDtime[i]<float(lnMJD.strip())) and (i<(len(obsAll)+2))):
            if i==(len(obsAll)+1):
                print('Date out of bounds')
                sys.exit()
                
            i=i+1
        #print i
        observations.append(obsAll[i])
        #obsNo=obsNo+1
        
    return observations

def readObsLocFile(obsLocFile, mjdCRFile, CRMin, CRMax):
    #Read STEREOobs file
    fObs=open(obsLocFile,'r')
    vObsLines=fObs.readlines()
    fObs.close()
    
    MJDfiles=open(mjdCRFile)
    MJDvalues=MJDfiles.readlines()
    MJDfiles.close()
    
    #MJDvalues.reverse()
    
    obsAll=[]
    obsMJDtime=[]
    #print len(vObsLines)
    for ln in vObsLines:
        s=ln.split()
        if (s[0].strip()!='YEAR'):
            #print(s)  
            MJDtime=gregToMJD(int(s[0].strip()),int(s[1].strip()),0)
            #print(s[0].strip()+', '+s[1].strip()+', '+str(MJDtime))
        
            obsMJDtime.append(float(MJDtime))
            obsAll.append(float(s[4].strip()))
    
    #print len(obsAll)
    obsLoc=[]
    i=0
    #j=0
    obsNo=0    
    for lnMJD in MJDvalues:
        #i=j
        #print(lnMJD)
        MJDsplit=lnMJD.split()
        #print(str(CRMin)+'/'+MJDsplit[0].strip()+'/'+str(CRMax))
        if (CRMin<=int(MJDsplit[0].strip())<=CRMax):
            #print(str(obsMJDtime[i])+'/'+MJDsplit[1].strip())
            while ((obsMJDtime[i]<float(MJDsplit[1].strip())) and (i<(len(obsAll)+2))):
                if i==(len(obsAll)+1):
                    print('Date out of bounds')
                    sys.exit()       
                i=i+1
            print(obsAll[i]) 
            #print i
            if obsAll[i]>180:
                obsAll[i]=obsAll[i]-360

            obsLoc.append(float(obsAll[i]))
    #obsNo=obsNo+1
        
    return obsLoc

def readObsLocFileWindows(obsLocFile, mjdCRFile, noOfWindows):
    #Read STEREOobs file
    fObs=open(obsLocFile,'r')
    vObsLines=fObs.readlines()
    fObs.close()
    
    MJDfiles=open(mjdCRFile)
    MJDvalues=MJDfiles.readlines()
    MJDfiles.close()
    
    #MJDvalues.reverse()
    
    obsAll=[]
    obsMJDtime=[]
    #print len(vObsLines)
    for ln in vObsLines:
        s=ln.split()
        if (s[0].strip()!='YEAR'):
            #print(s)  
            MJDtime=gregToMJD(int(s[0].strip()),int(s[1].strip()),0)
            #print(s[0].strip()+', '+s[1].strip()+', '+str(MJDtime))
        
            obsMJDtime.append(float(MJDtime))
            obsAll.append(float(s[4].strip()))
    
    #print len(obsAll)
    obsLoc=[]
    i=0
    #j=0
    obsNo=0    
    for lnMJD in MJDvalues:
        #i=j
        #print(lnMJD)
        MJDsplit=lnMJD.split()
        #print(str(CRMin)+'/'+MJDsplit[0].strip()+'/'+str(CRMax))
        if (0<=int(MJDsplit[0].strip())<=noOfWindows):
            #print(str(obsMJDtime[i])+'/'+MJDsplit[1].strip())
            while ((obsMJDtime[i]<float(MJDsplit[1].strip())) and (i<(len(obsAll)+2))):
                if i==(len(obsAll)+1):
                    print('Date out of bounds')
                    sys.exit()       
                i=i+1
            #print(obsAll[i]) 
            #print i
            if obsAll[i]>180:
                obsAll[i]=obsAll[i]-360

            obsLoc.append(float(obsAll[i]))
    #obsNo=obsNo+1
        
    return obsLoc
    
def gregToMJD(year,dayOfYear,hour):
    #Convert between year, day of year and hour to MJD
    #Find month and day of month from day of year
    #Define length of months
    if mod(year,4)==0:
        #dayMonths=[31,29,31,30,31,30,31,31,30,31,30,31]
        daysOfStartMonth=[0,31,60,91,121,152,182,213,244,274,305,335,366]
    else:
        #dayMonths=[31,28,31,30,31,30,31,31,30,31,30,31]
        daysOfStartMonth=[0,31,59,90,120,151,181,212,243,273,304,334,365]
    
    a=0
    m=1
    while ((a==0) or (m==14)):
        if (m==13):
            print('m=13: While loop overdone: Day no. not found')
            sys.exit()
        elif (daysOfStartMonth[m-1]<dayOfYear<=daysOfStartMonth[m]):
            month=m
            a=1
        
        m=m+1
    
    day=dayOfYear-daysOfStartMonth[month-1]
    
    #timeObs=dt.datetime(year,month,day,hour,0,0)
    
    #julianTimeObs=gcal2jd(year,month,day,hour)
    #Julian day init
    #JYear=-4713
    #JMonth=1
    #JDay=1
    #JHour=12
    #jdate=dt.datetime(JYear,JMonth,JDay,JHour,0,0)
    
    #year=2017
    #month=12
    #day=18
    #hour=21
    #timeObs=dt.datetime(year,month,day,hour,0,0)
    
    julianDayObs=(367*year)-floor(7*(year+floor((month+9)/12.0))*0.25)-floor(0.75*(floor(0.01*(year+((month-9)/7.0)))+1))+floor((275*month)/9.0)+day+1721028.5
    
    julianTimeObs=julianDayObs+(hour)/24.0
    MJD=julianTimeObs-2400000.5
    #sys.exit()
    #julian.to_jd(timeObs,format='mjd')
    
    return MJD

def dateToMJD(date):
    #Convert between date (string in form DDMMYYYY) to MJD
    day=int(date[:2])
    month=int(date[2:4])
    year=int(date[4:])
    
    if not 0<day<32:
        print('Invalid day input in date, please check.')
        print('Day input='+str(day))
        sys.exit()
    
    if not 0<month<13:
        print('Invalid month input in date')
        print('Month input='+str(month))
        sys.exit()
        
    if not 2004<year<2025:
        print('This program only performs DA over the period of the STEREO missions and cannot assimilate observations from the future. Please check year input in date')
        print('Year input='+str(year))
        sys.exit()
        
    julianDayObs=(367*year)-floor(7*(year+floor((month+9)/12.0))*0.25)-floor(0.75*(floor(0.01*(year+((month-9)/7.0)))+1))+floor((275*month)/9.0)+day+1721028.5
    
    #julianTimeObs=julianDayObs+(hour)/24.0
    MJD=julianDayObs-2400000.5
    #sys.exit()
    #julian.to_jd(timeObs,format='mjd')
    
    return MJD


def makeMJDfile(MJDstart,noOfLon,fileMJDOut,daySolarRot=27):
    #Function to make MJD file for one solar rotation from MJD start of length daySolarRot
    
    #Define timestep as length of solar rotation divided by number of lon points during rotation
    MJDstep=(daySolarRot)/float(noOfLon)
    #print initMJD
    #print finalMJD
    
    #Initialise MJD output variable with MJD startdate provided by user
    MJD=zeros((noOfLon))
    MJD[0]=copy(MJDstart)
 
    #Calculate evenly spread MJD values over solar rotation
    for lo in range(noOfLon-1):
        MJD[lo+1]=MJD[lo]+MJDstep

    #MJD.reverse()
    #print MJD[-1]
    #sys.exit()
    
    #Write MJD values to required output file
    MJDfile = open(fileMJDOut, "w")
    for lo in range(noOfLon-1,-1,-1):
        MJDfile.write(str(MJD[lo])+'\n')
    
    MJDfile.close()
    #######################no return############
    
def makeVinEns27d(initMJD,initCR,noOfLon,MJDdir,MASensDir,outputDir,daySolarRot=27,noOfMASens=576):
    #Function to make initial ensemble for required time period from Matt Owens' MAS ensemble members by smashing the two consecutive CR's
    #ensemble over which the rotation takes place
    
    MJDstart=initMJD
    finalMJD=initMJD+daySolarRot
    CRneeded=[int(initCR),int(initCR)+1]
    
    #MJDstartArray=initMJD
    fileMJD=MJDdir+'\\MJD_'+str(int(MJDstart))+'.dat'
    fMJD=open(fileMJD,'r')
    MJDcoord=fMJD.readlines()
    fMJD.close()

    startCoord=0
    for i in range(len(MJDcoord)):
        #print(float(MJDcoord[i].strip()))
        if (float(MJDcoord[i].strip())<MJDstart):
            startCoord=startCoord+1
    
    noOfCoordsReq=128
    endCoord=startCoord+noOfCoordsReq
    
    #Define file names
    file1=MASensDir+'\\vin_ensemble_CR'+str(CRneeded[0])+'.dat'
    file2=MASensDir+'\\vin_ensemble_CR'+str(CRneeded[1])+'.dat'
    
    fileOut=outputDir+'\\vin_ensemble_MJDstart'+str(int(MJDstart))+'.dat'
    
    #Read in MAS ensemble files for each CR
    f1=open(file1,'r')
    vInLines1=f1.readlines()
    f1.close()
    
    f2=open(file2,'r')
    vInLines2=f2.readlines()
    f2.close()
    
    #print len(vInLines)
    #Convert strings to velocities to be used in the ensemble
    vIn=zeros((576,len(vInLines1)+len(vInLines2)))  #One long vector containing all carrington longitude values over the two CRs of interest
    vOut=zeros((576,noOfCoordsReq))  #A vector to extract only the relevant carrinton longitudes for our solar Rotation
    
    #Read in solar wind speed ensembles
    for i in range(len(vInLines1)+len(vInLines2)):
        #print i
        #print vInLines[i]
        if i<len(vInLines1):
            vIn[:,i]=(vInLines1[i]).split('   ')[1:]
        else:
            vIn[:,i]=(vInLines2[i-len(vInLines1)]).split('   ')[1:]
    
    #Extract relevant values
    vOut=copy(vIn[:,startCoord:endCoord])
    
    #Save the new MAS ensemble for our solar rotation in the output file
    fOut=open(fileOut,'w')
    savetxt(fOut,transpose(vOut))
    fOut.close()
    
    ################################no return##########################
    
def createMeanAndBMatrix(deltaPhiDeg,locRad,noOfLonPoints,initFile,noOfMASens,outputDir):
    #Read input ensemble file
   # fVin=open(initFile,'r')
   # vInLines=fVin.readlines()
    #fVin.close()

    #print len(vInLines)
    vIn=zeros((576,noOfLonPoints))
    vIn=transpose(loadtxt(initFile))
    
    #for i in range(len(vInLines)):
        #print i
        #print vInLines[i]
    #    vIn[:,i]=(vInLines[i]).split('   ')[1:]
    
    B=zeros((noOfLonPoints,noOfLonPoints))
    
    #Generate perturbation ensemble (X-mean(X))
    pertEns=zeros((shape(vIn)))
    meanEns=zeros((noOfLonPoints))
    meanEns=vIn.mean(axis=0)
    
    pertEns=vIn-meanEns
    
    B=(1/(576.0-1))*transpose(pertEns).dot(pertEns)
    
    '''
    figure(3)
    im=contourf(range(128),range(128),B,cmap=cm.seismic,levels=range(-21000,21000,200))
    cobar=colorbar(im)
    cobar.set_label('km/s')
    xlabel('Carrington longitude coordinate')
    ylabel('Carrington longitude coordinate')
    title('Unlocalised B matrix')
    savefig(outputDir+'B_unloc.png')
    show()'''
    
    if locRad!=0:
        locMat=zeros((shape(B)))
    
        for i in range(noOfLonPoints):
            for j in range(noOfLonPoints):
                distBetweenPoints=min(abs((i-j)*deltaPhiDeg),360-abs((i-j)*deltaPhiDeg))
                locMat[i,j]=exp(-(distBetweenPoints**2)/(2*(locRad**2)))
                
        B=locMat*copy(B) 
        '''
        figure(4)
        im=contourf(range(128),range(128),B,cmap=cm.seismic,levels=range(-21000,21000,200))
        cobar=colorbar(im)
        cobar.set_label('km/s')
        xlabel('Carrington longitude coordinate')
        ylabel('Carrington longitude coordinate')
        title('Localised B matrix')
        savefig(outputDir+'B_loc.png')
        show()               
        '''
    #sys.exit()
    return meanEns[:noOfLonPoints], B[:noOfLonPoints,:noOfLonPoints]
    
def forwardRadModel(v,v0,r,rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints):
    #print 'Hi'
    velIn=zeros((noOfLatPoints,noOfLonPoints))
    velOut=zeros((noOfLatPoints,noOfLonPoints))
    
    velIn=v.reshape(noOfLatPoints,noOfLonPoints)
    
    #Execute the model equation for each longitude
    for phi in range(noOfLonPoints):
        if rIndex!=0:
            velOut[:,phi-1]=copy(velIn[:,phi-1])+(deltaR*solarRotFreq/deltaPhi)*((velIn[:,phi]-velIn[:,phi-1])/velIn[:,phi-1])-alpha*v0[:,phi-1]*(1-exp(-(r-1)/rH))+alpha*v0[:,phi-1]*(1-exp(-r/rH))
            #print (deltaR*solarRotFreq/deltaPhi)*((velIn[:,phi]-velIn[:,phi-1])/velIn[:,phi-1])
        else:
            velOut[:,phi-1]=copy(velIn[:,phi-1])+(deltaR*solarRotFreq/deltaPhi)*((velIn[:,phi]-velIn[:,phi-1])/velIn[:,phi-1])+alpha*v0[:,phi-1]*(1-exp(-r/rH))

    return velOut.reshape(noOfLatPoints*noOfLonPoints)

def forwardTimeModel(v,noOfLatPoints,noOfLonPoints):
    velIn=zeros((noOfLatPoints,noOfLonPoints))
    velOut=zeros((noOfLatPoints,noOfLonPoints))
    
    velIn=v
    
    for phi in range(noOfLonPoints):
        velOut[:,phi-1]=velIn[:,phi]
    
    return velOut
    
def TLMRad(forwardState,v,v0,r,rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints):
    #Calculate Fv
    #Tangent linear model for forwardRadModel
    
    velIn=zeros((noOfLatPoints,noOfLonPoints))
    velOut=zeros((noOfLatPoints,noOfLonPoints))
    velIn=v
    stateVect=forwardState.reshape(noOfLatPoints,noOfLonPoints)

    CFL=deltaR*solarRotFreq/deltaPhi
    if rIndex!=0:
        for phi in range(noOfLonPoints):
            #Test code if stateVect=v_r
            #velOut[:,phi]=velIn[:,phi]+(alpha*r/rH)*v0[:,phi]*exp(-r/rH)
            velOut[:,phi-1]=(1-CFL*(velIn[:,phi]/(velIn[:,phi-1]*velIn[:,phi-1])))*stateVect[:,phi-1]+CFL*(1.0/velIn[:,phi-1])*stateVect[:,phi]
            '''
            #Check if phi=0 to ensure r component is skipped in stateVect
            if phi==0:
                velOut[:,-1]=(1-CFL*(velIn[:,0]/(velIn[:,-1]*velIn[:,-1])))*stateVect[:,noOfLonPoints-1]+CFL*(1.0/velIn[:,-1])*stateVect[:,0]+(alpha*v0[:,-1]/rH)*stateVect[:,noOfLonPoints]*exp(-r/rH)
            else:
                velOut[:,phi-1]=(1-CFL*(velIn[:,phi]/(velIn[:,phi-1]*velIn[:,phi-1])))*stateVect[:,phi-1]+CFL*(1.0/velIn[:,phi-1])*stateVect[:,phi]+(alpha*v0[:,phi-1]/rH)*stateVect[:,noOfLonPoints]*exp(-r/rH)
            '''
    else:
        #If r=0, different differentiation result from vAcc component
        for phi in range(noOfLonPoints):
            #Test code if stateVect=v_r
            #velOut[:,phi]=velIn[:,phi]+alpha*velIn[:,phi]*(1-exp(-r/rH))+(alpha*r/rH)*velIn[:,phi]*exp(-r/rH)
            if ((velIn[:,phi-1]*velIn[:,phi-1]).any())==0:
                print('velFail (phi-1)=: '+str(phi))
            if (rH==0):
                print('rH Fail')
            velOut[:,phi-1]=(1-CFL*(velIn[:,phi]/(velIn[:,phi-1]*velIn[:,phi-1]))+alpha*(1-exp(-r/rH)))*stateVect[:,phi-1]+CFL*(1.0/velIn[:,phi-1])*stateVect[:,phi]
            '''
            #Check if phi=0 to ensure r component is skipped in stateVect
            if phi==0:
                velOut[:,-1]=(1-CFL*(velIn[:,0]/(velIn[:,-1]*velIn[:,-1]))+alpha*(1-exp(-r/rH)))*stateVect[:,noOfLonPoints-1]+CFL*(1.0/velIn[:,-1])*stateVect[:,0]+(alpha*v0[:,-1]/rH)*stateVect[:,noOfLonPoints]*exp(-r/rH)
            else:
                velOut[:,phi-1]=(1-CFL*(velIn[:,phi]/(velIn[:,phi-1]*velIn[:,phi-1]))+alpha*(1-exp(-r/rH)))*stateVect[:,phi-1]+CFL*(1.0/velIn[:,phi-1])*stateVect[:,phi]+(alpha*v0[:,phi-1]/rH)*stateVect[:,noOfLonPoints]*exp(-r/rH)
            '''
    return velOut.reshape(noOfLatPoints*noOfLonPoints)

def adjRad(lambdaVar,v,v0,r,rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints):
    #Calculate F^T v
    #Adjoint code for forwardRad model
    velIn=zeros((noOfLatPoints,noOfLonPoints))
    adjOut=zeros((noOfLatPoints,noOfLonPoints))#+1))
    
    velIn=v
    adjVar=lambdaVar.reshape(noOfLatPoints,noOfLonPoints)
    
    CFL=deltaR*solarRotFreq/deltaPhi
    #solarRotFreq=1./(27.0*24.0*3600.0)
    
    if rIndex!=0:
        for phi in range(noOfLonPoints):
            '''
            #Testing code when adjVar=v_r
            if phi==0:
                adjOut[:,noOfLonPoints-1]=velIn[:,-1]+CFL*((velIn[:,-1]-velIn[:,0])/velIn[:,-1])
            else:
                adjOut[:,phi-1]=velIn[:,phi-1]+CFL*((velIn[:,phi-1]-velIn[:,phi])/velIn[:,phi-1])
            adjOut[:,noOfLonPoints]=adjOut[:,noOfLonPoints]+(alpha/rH)*(v0[:,phi]*velIn[:,phi])
            '''
            adjOut[:,phi-1]=(1-CFL*(velIn[:,phi]/(velIn[:,phi-1]*velIn[:,phi-1])))*adjVar[:,phi-1]+(CFL*(1/velIn[:,phi-2]))*adjVar[:,phi-2]
            '''
            #Check whether or not phi=0 in order to ensure final adjoint velocity is in correct location
            if phi==0:
                adjOut[:,noOfLonPoints-1]=(1-CFL*(velIn[:,0]/(velIn[:,-1]*velIn[:,-1])))*adjVar[:,-1]+(CFL*(1/velIn[:,-2]))*adjVar[:,-2]
            else:
                adjOut[:,phi-1]=(1-CFL*(velIn[:,phi]/(velIn[:,phi-1]*velIn[:,phi-1])))*adjVar[:,phi-1]+(CFL*(1/velIn[:,phi-2]))*adjVar[:,phi-2]
            '''
            #Add contribution for radius component of adjoint for current longitude
            #adjOut[:,noOfLonPoints]=adjOut[:,noOfLonPoints]+((alpha/rH)*exp(-r/rH)*v0[:,phi]*adjVar[:,phi])
            
    else:
        #When r=0, there is an additional component in the v_acc that must be accounted for
        for phi in range(noOfLonPoints):
            '''
            #Testing code when adjVar=v_r
            if phi==0:
                adjOut[:,noOfLonPoints-1]=velIn[:,-1]+CFL*((velIn[:,-1]-velIn[:,0])/velIn[:,-1])+(alpha*v0[:,-1]*(1-exp(-r/rH)))
            else:
                adjOut[:,phi-1]=velIn[:,phi-1]+CFL*((velIn[:,phi-1]-velIn[:,phi])/velIn[:,phi-1])+(alpha*v0[:,phi-1]*(1-exp(-r/rH)))
            adjOut[:,noOfLonPoints]=adjOut[:,noOfLonPoints]+(alpha/rH)*(v0[:,phi]*velIn[:,phi])
            '''
            adjOut[:,phi-1]=(CFL*(1.0/velIn[:,phi-2])*adjVar[:,phi-2])+(1-CFL*(velIn[:,phi]/(velIn[:,phi-1]*velIn[:,phi-1]))+alpha*(1-exp(-r/rH)))*adjVar[:,phi-1]
            '''
            #Check whether or not phi=0 in order to ensure final adjoint velocity is in correct location
            if phi==0:
                adjOut[:,noOfLonPoints-1]=(1-CFL*(velIn[:,0]/(velIn[:,-1]*velIn[:,-1]))+alpha*(1-exp(-r/rH)))*adjVar[:,-1]+(CFL*(1/velIn[:,-2]))*adjVar[:,-2]
            else:
                adjOut[:,phi-1]=(CFL*(1.0/velIn[:,phi-2])*adjVar[:,phi-2])+(1-CFL*(velIn[:,phi]/(velIn[:,phi-1]*velIn[:,phi-1]))+alpha*(1-exp(-r/rH)))*adjVar[:,phi-1]
                
            #Add contribution for radius component of adjoint for current longitude
            adjOut[:,noOfLonPoints]=adjOut[:,noOfLonPoints]+((alpha/rH)*exp(-r/rH)*v0[:,phi]*adjVar[:,phi])    
            '''
    #adjOut[:,noOfLonPoints]=((alpha/rH)*exp(r/rH))*(v0*adjVar).sum(axis=1)
    return adjOut.reshape(noOfLatPoints*noOfLonPoints)
    
def adjointTest(xIn,yIn,v,v0,r,rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints):
    #Code to test if <Fx,y>=<x,F^T y> holds for adjoint and TLM scripts
    #x=zeros((noOfLatPoints,noOfLonPoints+1))
    x=zeros((noOfLatPoints,noOfLonPoints))
    y=zeros((noOfLatPoints,noOfLonPoints))
    
    x=xIn.reshape((noOfLatPoints*noOfLonPoints))
    y=yIn.reshape((noOfLatPoints*noOfLonPoints))
    
    #print xIn[5,:]
    #print yIn[5,:]
    TLMx=zeros((noOfLatPoints,noOfLonPoints))
    Adjy=zeros((noOfLatPoints,noOfLonPoints))
    #Adjy=zeros((noOfLatPoints,noOfLonPoints+1))
    
    #Run TLM model on x and Adjoint model on y
    TLMx=TLMRad(x,v,v0,r,rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
    Adjy=adjRad(y,v,v0,r,rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
    
    #print 'TLMx='+str(TLMx[5,:])
    #print 'Adjy='+str(Adjy[5,:])
    
    #Reshape matrices into vectors in order to do dot products
    TLMxReshape=TLMx#.reshape(noOfLatPoints*noOfLonPoints,1)
    AdjyReshape=Adjy#.reshape(noOfLatPoints*(noOfLonPoints),1)
    #AdjyReshape=Adjy.reshape(noOfLatPoints*(noOfLonPoints+1),1)
    
    xReshape=x.reshape(noOfLatPoints*(noOfLonPoints),1)
    #xReshape=x.reshape(noOfLatPoints*(noOfLonPoints+1),1)
    yReshape=y.reshape(noOfLatPoints*noOfLonPoints,1)
    
    #Do dot products
    yDotTLMx=(TLMxReshape.transpose()).dot(yReshape)
    AdjyDotx=(xReshape.transpose()).dot(AdjyReshape)
    
    #Print out results of test
    print ('<y,Fx>= '+str(yDotTLMx))
    print ('<F*y,x>='+str(AdjyDotx))
    print ('<y,Fx>-<F*y,x>='+str(yDotTLMx-AdjyDotx))
    
    return yDotTLMx-AdjyDotx

def obsOp(lengthState, obsToBeTakenSpace,lonCoord):
    #Generates observation operator H_{obsTime}(lengthState,noOfObs)
    lengthObsSpace=len(obsToBeTakenSpace)
    lengthLon=len(lonCoord)
    H=zeros((lengthObsSpace*lengthLon,lengthState))
                
    #Use identity to start with
    for lo in range(lengthLon):
        for i in range(lengthObsSpace):
            j=obsToBeTakenSpace[i]
            H[int((lo*lengthObsSpace)+i),int(mod(j+lonCoord[lo],lengthState))]=1
            #print str(lo*(lengthObsSpace)+i)+', '+str(mod(j-lonCoord[lo],lengthState))
    '''        
    figure(50)        
    imshow(H)
    show()
    sys.exit()
    '''
    #print('obsOp='+str(gauss(1,2)))    
    return H

def makeGrad(B,R,H,x,xb,v0,y,radObs,noOfObs,r,rH,deltaR,deltaPhi,alpha,solarRotFreq,noOfRadPoints,noOfLatPoints,noOfLonPoints):
    #Define variables that are required
    gridDim=noOfLatPoints*noOfLonPoints
    
    #stateDim=noOfRadPoints*noOfLatPoints*noOfLonPoints
    #radObs=215
    radObsLocation=radObs
    #zeros((obsRadDim,noOfObs))
    #obsDim=zeros((len(radObs)))
    #obsDim[0]=noOfObs
    
    velNonLin=zeros((noOfRadPoints,noOfLatPoints,noOfLonPoints))
    forwardStateVar=zeros((noOfRadPoints,gridDim))
    velNonLin[0,:,:]=copy(x.reshape(noOfLatPoints,noOfLonPoints))
    forwardStateVar[0,:]=x
    #Create adjoint variables
    #adjOut=zeros((noOfRadPoints+1,noOfLatPoints*noOfLonPoints))
    adjVar=zeros((noOfRadPoints+1,gridDim))
    
    #Create variable to store gradient of cost function
    gradJ=zeros((gridDim))
    
    #Run forward model
    for rIndex in range(1,noOfRadPoints):
        #print rIndex
        forwardStateVar[rIndex,:]=forwardRadModel(velNonLin[rIndex-1,:,:],v0,r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
        #velTLM[rIndex,:,:]=TLMRad(forwardStateVect[rIndex-1,:,:],velNonLin[rIndex-1,:,:],v0,r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
        velNonLin[rIndex,:,:]=copy(forwardStateVar[rIndex,:].reshape(noOfLatPoints,noOfLonPoints))
    
    #Run adjoint model
    for rIndex in range(noOfRadPoints-1,-1,-1):
        #print rIndex
        if rIndex==noOfRadPoints-1: #radObsLocation.any()==rIndex:
            #i=radObsLocation.where(rIndex)
            #print 'H='+str(H.transpose())
            #print 'y='+str(y)
            #print 'Hv='+str(H.dot(velNonLin[rIndex,:,:].reshape(gridDim)))
            #print 'R='+str(R[0,0])
            #print 'R^{-1}='+str(pinv(R))
            adjVar[rIndex,:]=adjRad(adjVar[rIndex+1,:],velNonLin[rIndex,:,:],v0,r[rIndex],rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)+H.transpose().dot(pinv(R)).dot(y-H.dot(velNonLin[rIndex,:,:].reshape(gridDim)))
        else:
            adjVar[rIndex,:]=adjRad(adjVar[rIndex+1,:],velNonLin[rIndex,:,:],v0,r[rIndex],rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
    
    if radObsLocation.any()==0:
        i=radObsLocation.where(0)
        #gradJ=(B^(-1)).dot(x-xb)-adjRad(adjVar[rIndex,:],velNonLin,v0,r,rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)-H.transpose().dot(R^(-1)).dot(y[i]-H.dot(velNonLin[rIndex,:,:]))
    else:
        gradJ=(inv(B)).dot(x-xb)-adjRad(adjVar[0,:],velNonLin[0,:,:],v0,r[0],0,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
        
    #print gradJ
    
    return gradJ

def anaAdjSoln(B,R,H,xb,v0,y,radObs,noOfObs,r,rH,deltaR,deltaPhi,alpha,solarRotFreq,noOfRadPoints,noOfLatPoints,noOfLonPoints,iterNo):
    #Define variables that are required
    gridDim=noOfLatPoints*noOfLonPoints
    
    #stateDim=noOfRadPoints*noOfLatPoints*noOfLonPoints
    #radObs=215
    radObsLocation=radObs
    #zeros((obsRadDim,noOfObs))
    #obsDim=zeros((len(radObs)))
    #obsDim[0]=noOfObs
    
    velNonLin=zeros((noOfRadPoints,noOfLatPoints,noOfLonPoints))
    forwardStateVar=zeros((noOfRadPoints,gridDim))
    velNonLin[0,:,:]=copy(xb.reshape(noOfLatPoints,noOfLonPoints))
    forwardStateVar[0,:]=xb
    #Create adjoint variables
    #adjOut=zeros((noOfRadPoints+1,noOfLatPoints*noOfLonPoints))
    adjVar=zeros((noOfRadPoints+1,gridDim))
    
    #Create variable to store gradient of cost function
    gradJ=zeros((gridDim))
    
    #Run forward model
    for rIndex in range(1,noOfRadPoints):
        #print rIndex
        forwardStateVar[rIndex,:]=forwardRadModel(velNonLin[rIndex-1,:,:],v0,r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
        velNonLin[rIndex,:,:]=copy(forwardStateVar[rIndex,:].reshape(noOfLatPoints,noOfLonPoints))
    '''
    for rIndex in range(1,noOfRadPoints):
        #print rIndex
        forwardStateVar[rIndex,:]=forwardRadModel(velNonLin[rIndex-1,:,:],v0,r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
        #velTLM[rIndex,:,:]=TLMRad(forwardStateVect[rIndex-1,:,:],velNonLin[rIndex-1,:,:],v0,r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
        velNonLin[rIndex,:,:]=copy(forwardStateVar[rIndex,:].reshape(noOfLatPoints,noOfLonPoints))
    '''
    #Run adjoint model
    for rIndex in range(noOfRadPoints-1,0,-1):
        #print rIndex
        if rIndex==noOfRadPoints-1: #radObsLocation.any()==rIndex:
            #i=radObsLocation.where(rIndex)
            #print 'H='+str(H.transpose())
            #print 'y='+str(y)
            #print 'y-Hv='+str(y-H.dot(velNonLin[rIndex,:,:].reshape(gridDim)))
            #print 'R='+str(R[0,0])
            #print 'R^{-1}='+str(pinv(R))
            #print H.transpose().dot(pinv(R)).dot(y-H.dot(velNonLin[-1,:,:].reshape(gridDim)))
            #sys.exit()
            #print H.transpose().dot(pinv(R)).dot(y-H.dot(velNonLin[rIndex,:,:].reshape(gridDim)))
            #sys.exit()
            adjVar[rIndex,:]=adjRad(adjVar[rIndex+1,:],velNonLin[rIndex,:,:],v0,r[rIndex],rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)+H.transpose().dot(pinv(R)).dot(y-H.dot(velNonLin[rIndex,:,:].reshape(gridDim)))
        else:
            adjVar[rIndex,:]=adjRad(adjVar[rIndex+1,:],velNonLin[rIndex,:,:],v0,r[rIndex],rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
    
    levels=arange(-1,1.5005,0.01)
    fig=figure(54000+iterNo)
    im=contourf(arange(0,359,2.81),range(30,30+noOfRadPoints+1),adjVar,cmap=cm.seismic,levels=levels)
    colorbar(im)
    #adjVar=zeros((noOfRadPoints+1,gridDim))    
    #show()
    savefig('C:\\Users\\Matthew\\Desktop\\Work\\Python\\spaceWeatherModel\\output\\CR2100\\LargeErrors\\adjVar_'+str(iterNo)+'.png')
    close(fig)
    
    #Update state
    x=xb+B.dot(adjRad(adjVar[1,:],velNonLin[0,:,:],v0,r[0],0,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints))
    #print 'BF^T \lambda = '+str(xb+B.dot(adjRad(adjVar[1,:],velNonLin[0,:,:],v0,r[0],0,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)))
    #print 'B[0,:]='+str(B[0,:])
    #print '\lambda_0='+str(adjRad(adjVar[1,:],velNonLin[0,:,:],v0,r[0],0,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints))
    #sys.exit()
    #print x

    return x
    
def calcCostFunc(B,R,H,x,xb,v,y,radObs,nRadObs,noOfObs,noOfLatPoints,noOfLonPoints):
    costFuncObs=0
    costFuncBack=(transpose(x-xb).dot(pinv(B))).dot(x-xb)
    gridDim=noOfLatPoints*noOfLonPoints
    #costFuncObs=(transpose(y-H.dot(v[182,:,:].reshape(gridDim))).dot(pinv(R))).dot(y-H.dot(v[182,:,:].reshape(gridDim)))
    
    for r in range(len(radObs)):
        if r==0:
            costFuncObs=(transpose(y[0:int(nRadObs[0])]-H[0:int(nRadObs[0]),:].dot(v[int(radObs[0]),:,:].reshape(gridDim))).dot(pinv(R[0:int(nRadObs[0]),0:int(nRadObs[0])]))).dot(y[0:int(nRadObs[0])]-H[0:int(nRadObs[0]),:].dot(v[int(radObs[0]),:,:].reshape(gridDim)))
        else:
            costFuncObs=costFuncObs+(transpose(y[int(nRadObs[r-1]):int(nRadObs[r])]-H[int(nRadObs[r-1]):int(nRadObs[r]),:].dot(v[int(radObs[r]),:,:].reshape(gridDim))).dot(pinv(R[int(nRadObs[r-1]):int(nRadObs[r]),int(nRadObs[r-1]):int(nRadObs[r])]))).dot(y[int(nRadObs[r-1]):int(nRadObs[r])]-H[int(nRadObs[r-1]):int(nRadObs[r]),:].dot(v[int(radObs[r]),:,:].reshape(gridDim)))
    
    #(transpose(y-H.dot(velNonLin[-1,:,:].reshape(gridDim))).dot(pinv(R))).dot(y-H.dot(velNonLin[-1,:,:].reshape(gridDim)))
    #z=v[-1,:,:].reshape(gridDim)
    #print 'y-Hx='+str(z[:5])
    print ('Jb='+str(costFuncBack))
    print ('Jo='+str(costFuncObs))
    print ('costFunc='+str(costFuncBack+costFuncObs))
    
    return costFuncBack+costFuncObs
    
def makeGradCG(xb,B,R,H,xBIter0,v,y,radObs,nRadObs,noOfObs,r,rH,deltaR,deltaPhi,alpha,solarRotFreq,noOfRadPoints,noOfLatPoints,noOfLonPoints,iterNo):
    #Define variables that are required
    gridDim=noOfLatPoints*noOfLonPoints
    
    #stateDim=noOfRadPoints*noOfLatPoints*noOfLonPoints
    #radObs=215
    radObsLocation=radObs
    #zeros((obsRadDim,noOfObs))
    #obsDim=zeros((len(radObs)))
    #obsDim[0]=noOfObs
    
    velNonLin=zeros((noOfRadPoints,noOfLatPoints,noOfLonPoints))
    forwardStateVar=zeros((noOfRadPoints,gridDim))
    #velNonLinAlphaTest=zeros((noOfRadPoints,noOfLatPoints,noOfLonPoints))
    #forwardStateVarAlphaTest=zeros((noOfRadPoints,gridDim))
    velNonLin[0,:,:]=copy(xb.reshape(noOfLatPoints,noOfLonPoints))
    forwardStateVar[0,:]=xb
    #Create adjoint variables
    #adjOut=zeros((noOfRadPoints+1,noOfLatPoints*noOfLonPoints))
    adjVar=zeros((noOfRadPoints+1,gridDim))
    
    #Create variable to store gradient of cost function
    gradJ=zeros((gridDim))
    
    #Run forward model
    for rIndex in range(1,noOfRadPoints):
        #print rIndex
        forwardStateVar[rIndex,:]=forwardRadModel(velNonLin[rIndex-1,:,:],velNonLin[0,:,:],r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
        velNonLin[rIndex,:,:]=copy(forwardStateVar[rIndex,:].reshape(noOfLatPoints,noOfLonPoints))
        
    #velNonLin=copy(v)
    #forwardStateVector=copy(v.reshape(noOfRadPoints,noOfLatPoints*noOfLonPoints))
    '''
    for rIndex in range(1,noOfRadPoints):
        #print rIndex
        forwardStateVar[rIndex,:]=forwardRadModel(velNonLin[rIndex-1,:,:],v0,r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
        #velTLM[rIndex,:,:]=TLMRad(forwardStateVect[rIndex-1,:,:],velNonLin[rIndex-1,:,:],v0,r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
        velNonLin[rIndex,:,:]=copy(forwardStateVar[rIndex,:].reshape(noOfLatPoints,noOfLonPoints))
    '''
    #Run adjoint model
    for rIndex in range(noOfRadPoints-1,0,-1):
        #print rIndex
        if rIndex in radObs: #noOfRadPoints-1: #radObsLocation.any()==rIndex:
            i=list(radObs).index(rIndex)
            #print 'rIndex='+str(rIndex)
            #i=radObsLocation.where(rIndex)
            #print 'H='+str(H.transpose())
            #print 'y='+str(y)
            #print 'y-Hv='+str(y-H.dot(velNonLin[rIndex,:,:].reshape(gridDim)))
            #print 'R='+str(R[0,0])
            #print 'R^{-1}='+str(pinv(R))
            #print H.transpose().dot(pinv(R)).dot(y-H.dot(velNonLin[-1,:,:].reshape(gridDim)))
            #sys.exit()
            #print H.transpose().dot(pinv(R)).dot(y-H.dot(velNonLin[rIndex,:,:].reshape(gridDim)))
            #sys.exit()
            if i==0:
                adjVar[rIndex,:]=adjRad(adjVar[rIndex+1,:],velNonLin[rIndex,:,:],velNonLin[0,:,:],r[rIndex],rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)+H[:int(nRadObs[0]),:].transpose().dot(pinv(R[:int(nRadObs[0]),:int(nRadObs[0])])).dot(y[:int(nRadObs[0])]-H[:int(nRadObs[0]),:].dot(velNonLin[rIndex,:,:].reshape(gridDim)))
                #print 'y-Hx='+str(y[:int(nRadObs)]-H[:int(nRadObs[0]),:].dot(velNonLin[rIndex,:,:].reshape(gridDim)))
                #+H[:int(nRadObs[0]),:].transpose().dot((pinv(R[:int(nRadObs[0]),:int(nRadObs[0])]))).dot(y[:int(nRadObs[0])]-H[:int(nRadObs[0]),:].dot(velNonLin[rIndex,:,:].reshape(gridDim)))#.dot(pinv(R)).dot(y-H.dot(velNonLin[rIndex,:,:].reshape(gridDim)))
            else:
                adjVar[rIndex,:]=adjRad(adjVar[rIndex+1,:],velNonLin[rIndex,:,:],velNonLin[0,:,:],r[rIndex],rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)+H[int(nRadObs[i-1]):int(nRadObs[i]),:].transpose().dot(pinv(R[int(nRadObs[i-1]):int(nRadObs[i]),int(nRadObs[i-1]):int(nRadObs[i])])).dot(y[int(nRadObs[i-1]):int(nRadObs[i])]-H[int(nRadObs[i-1]):int(nRadObs[i]),:].dot(velNonLin[rIndex,:,:].reshape(gridDim)))#.dot(pinv(R)).dot(y-H.dot(velNonLin[rIndex,:,:].reshape(gridDim)))
        else:
            adjVar[rIndex,:]=adjRad(adjVar[rIndex+1,:],velNonLin[rIndex,:,:],velNonLin[0,:,:],r[rIndex],rIndex,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
    '''
    levels=arange(-1,1.5005,0.01)
    fig=figure(54000+iterNo)
    im=contourf(arange(0,359,2.81),range(30,30+noOfRadPoints+1),adjVar,cmap=cm.seismic,levels=levels)
    colorbar(im)
    #adjVar=zeros((noOfRadPoints+1,gridDim))    
    #show()
    savefig('C:\\Users\\Matthew\\Desktop\\Work\\Python\\spaceWeatherModel\\output\\CR2100\\LargeErrors\\adjVar_'+str(iterNo)+'.png')
    close(fig)
    '''
    gradJ=transpose(matrix(pinv(B).dot(xb-xBIter0)-adjRad(adjVar[1,:],velNonLin[0,:,:],velNonLin[0,:,:],r[0],0,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)))
    
    return array(gradJ)[:,0]
    
def calcCostFuncForCG(xb,B,R,H,xBIter0,v,y,radObs,nRadObs,noOfObs,r,rH,deltaR,deltaPhi,alpha,solarRotFreq,noOfRadPoints,noOfLatPoints,noOfLonPoints,iterNo):
    costFuncObs=0
    costFuncBack=(transpose(xb-xBIter0).dot(pinv(B))).dot(xb-xBIter0)
    gridDim=noOfLatPoints*noOfLonPoints
    
    velNonLin=zeros((noOfRadPoints,noOfLatPoints,noOfLonPoints))
    forwardStateVar=zeros((noOfRadPoints,gridDim))
    
    #velNonLin=copy(v)
    forwardStateVar[0,:]=xb
    velNonLin[0,:,:]=xb.reshape(noOfLatPoints,noOfLonPoints)
    
    #Run forward model
    for rIndex in range(1,noOfRadPoints):
        #print rIndex
        forwardStateVar[rIndex,:]=forwardRadModel(velNonLin[rIndex-1,:,:],velNonLin[0,:,:],r[rIndex-1],rIndex-1,deltaR,deltaPhi,solarRotFreq,alpha,rH,noOfRadPoints,noOfLatPoints,noOfLonPoints)
        velNonLin[rIndex,:,:]=copy(forwardStateVar[rIndex,:].reshape(noOfLatPoints,noOfLonPoints))
    
    #print (x-xb)[:15]
    #print pinv(B)[:15,0]
    for i in range(len(radObs)):
        if i==0:
            #print 'Hi'
            costFuncObs=(transpose(y[:int(nRadObs[0])]-H[:int(nRadObs[0]),:].dot(velNonLin[int(radObs[0]),:,:].reshape(gridDim))).dot(pinv(R[:int(nRadObs[0]),:int(nRadObs[0])]))).dot(y[:int(nRadObs[0])]-H[:int(nRadObs[0]),:].dot(velNonLin[int(radObs[0]),:,:].reshape(gridDim)))
            #print 'y-Hx(CG)='+str(y[:int(nRadObs[0])]-H[:int(nRadObs[0]),:].dot(velNonLin[int(radObs[0]),:,:].reshape(gridDim)))
            #print 'costFuncObsCG='+str(costFuncObs)
        else:
            #print 'Bob'
            costFuncObs=costFuncObs+(transpose(y[int(nRadObs[i-1]):int(nRadObs[i])]-H[int(nRadObs[i-1]):int(nRadObs[i]),:].dot(velNonLin[int(radObs[i]),:,:].reshape(gridDim))).dot(pinv(R[int(nRadObs[i-1]):int(nRadObs[i]),int(nRadObs[i-1]):int(nRadObs[i])]))).dot(y[int(nRadObs[i-1]):int(nRadObs[i])]-H[int(nRadObs[i-1]):int(nRadObs[i]),:].dot(velNonLin[int(radObs[i]),:,:].reshape(gridDim)))
    
    #print '1: '+str(costFuncObs)
    #costFuncObs2=(transpose(y-H.dot(velNonLin[-1,:,:].reshape(gridDim))).dot(pinv(R))).dot(y-H.dot(velNonLin[-1,:,:].reshape(gridDim)))
    #print '2: '+str(costFuncObs2)
    #costFuncObs=(transpose(y-H.dot(velNonLin[-1,:,:].reshape(gridDim))).dot(pinv(R))).dot(y-H.dot(velNonLin[-1,:,:].reshape(gridDim)))
    #(transpose(y-H.dot(velNonLin[-1,:,:].reshape(gridDim))).dot(pinv(R))).dot(y-H.dot(velNonLin[-1,:,:].reshape(gridDim)))
    #z=v[-1,:,:].reshape(gridDim)
    #print 'y-Hx='+str(z[:5])
    
    #print 'Non_lin_Jb='+str(costFuncBack)
    #print 'Non_lin_Jo='+str(costFuncObs)
    #print 'Non_lin_costFunc='+str(costFuncBack+costFuncObs)
    
    return costFuncBack+costFuncObs