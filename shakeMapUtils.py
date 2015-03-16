# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 08:59:06 2014

Nick Horspool
GNS Science
December 2014

Functions to generate N Ground Motion Fields that have
spatially correlated residuals. 

1) Read shakemap uncertainty.xml file
2) Extract event-specific intra-event uncertainty (if calculated), or use GMPE one
3) For each site, randomly sample from a normal distribution with mean=0,
and standard deviation=intra-event uncertainty. This will give a log(residual).
4) Apply a spatial correlation model to the residuals (e.g. Jayaram and Baker)
5) Add residuals to mean and output data
6) Repeat N times

"""

import xml.etree.ElementTree as ET
import numpy as np
import scipy as sp
import math

def parseUncertaintyXML(xmlFile):
    """
    Function that reads in the shakemap uncertainty.xml file and returns an
    array of data
    Input: uncertainty.xml file
    Output: array of stdevs: x,y,PGA,PGV,MMI,PSA03,PSA10,PSA30 and dictionary
    of event specific intraevent uncertainty
    """
    
    intraEventDict = {}    
    gmStdDevs = []
    
    tree = ET.parse(xmlFile)
    root = tree.getroot()
    
    for child in root:
        name = child.tag.split('}')[1]
        # Get event specific data and add to dictionary
        if name == 'event_specific_uncertainty':
            
            IMT = child.attrib['name']
            sigma = float(child.attrib['value'])
            intraEventDict[IMT] = sigma

        if name == 'grid_data':
            
            # Read in large block of text from grid_data node
            # and reformat
            data = child.text.split('\n')
            # First and last line is blank when parsed for some reason
            for line in data[1:-1]:
                siteStdDevs = line.split(' ')
                gmStdDevs.append(siteStdDevs)

    return intraEventDict, gmStdDevs
    
    
def parseGridXML(xmlFile):
    """
    Function that reads in the shakemap grid.xml file and returns an
    array of data
    Input: grid.xml file
    Output: array of GM medians x,y,PGA,PGV,MMI,PSA03,PSA10,PSA30
    """
       
    gmMedians = []
    lons = []
    lats = []
    
    tree = ET.parse(xmlFile)
    root = tree.getroot()
    
    for child in root:
        name = child.tag.split('}')[1]

        if name == 'grid_data':
            # Read in large block of text from grid_data node
            # and reformat
            data = child.text.split('\n')
            # First line is blank
            for line in data[1:]:
                siteMedians = line.split(' ')
                #print siteMedians[2]
                if len(siteMedians) > 2:
                    gmMedians.append(float(siteMedians[2])/100.0)
                    lons.append(float(siteMedians[0]))
                    lats.append(float(siteMedians[1]))
                

    return np.array(lons), np.array(lats), np.array(gmMedians)
    

def parseStationData(xmlFile):
    """
    Function that reads in the shakemap stationlist.xml file and returns an
    array of lons and lats
    Input: stationlist.xml file
    Output: lons, lats (as numpy arrays)
    """
    
    lons = []
    lats = []
    
    tree = ET.parse(xmlFile)
    root = tree.getroot()
    
    for child in root:
        
        # Get event specific data and add to dictionary
        if child.tag == 'stationlist':
            
            for step_child in child:

                lons.append(float(step_child.attrib['lon']))
                lats.append(float(step_child.attrib['lat']))
            
    
    return np.array(lons), np.array(lats)

def prepare_coords(lons1, lats1, lons2, lats2):
    """
    Convert two pairs of spherical coordinates in decimal degrees
    to numpy arrays of radians. Makes sure that respective coordinates
    in pairs have the same shape.
    
    Code origin: Openquake
    """
    
    lons1 = np.array(np.radians(lons1))
    lats1 = np.array(np.radians(lats1))
    assert lons1.shape == lats1.shape
    lons2 = np.array(np.radians(lons2))
    lats2 = np.array(np.radians(lats2))
    assert lons2.shape == lats2.shape
    
    return lons1, lats1, lons2, lats2

def getDistances(lons1, lats1, lons2, lats2):
    """
    Get distances between two points on earths surface
    Code origin: Openquake
    """
    
    EARTH_RADIUS = 6371.0
    
    lons1, lats1, lons2, lats2 = prepare_coords(lons1, lats1, lons2, lats2)
    distance = np.arcsin(np.sqrt(
        np.sin((lats1 - lats2) / 2.0) ** 2.0
        + np.cos(lats1) * np.cos(lats2)
        * np.sin((lons1 - lons2) / 2.0) ** 2.0
    ).clip(-1., 1.))
    
    return (2.0 * EARTH_RADIUS) * distance


def getDistanceMatrix(lons,lats):
    """
    Calculate a distance matrix of sites.
    Input: 1D numpy array of length N of lons and lats where N is number of sites
    Code origin: Openquake
    """
    
    distances = getDistances(
            lons.reshape(lons.shape + (1, )),
            lats.reshape(lats.shape + (1, )),
            lons,
            lats
        )
        
    return np.matrix(distances, copy=False)

def getCorrelationModel(distances, imt, vs30Clustering):
    """
    Returns the correlation model for a set of distances, given the
    appropriate period
    :param numpy.ndarray distances:
        Distance matrix
    :param float period:
        Period of spectral acceleration
    :param vs30Clustering:
        True or False, determines correlation model
    Code Origin: Openquake
    """
    if imt == 'PGA':
        period = 0
        
    else:
        period = imt

    # formulae are from page 1700
    if period < 1:
        if not vs30Clustering:
            # case 1, eq. (17)
            b = 8.5 + 17.2 * period
        else:
            # case 2, eq. (18)
            b = 40.7 - 15.0 * period
    else:
        # both cases, eq. (19)
        b = 22.0 + 3.7 * period

    # eq. (20)
    return np.exp((- 3.0 / b) * distances)
    
    
def getSiteIntraSigma(eventIntras, gmpeStdDev):
    """
    Returns intra sigma term for each site that will be sampled in spatial correlation
    """
    
    eventIntraSigma = eventIntras['pga']
    lons = []
    lats = []
    siteIntraSigmas = []
    
    # Loop over each site to get data out
    for i in xrange(0, len(gmpeStdDev)):
            
        totalSigma = float(gmpeStdDev[i][2])
        gmpeInterSigma = float(gmpeStdDev[i][8])
        gmpeIntraSigma = float(gmpeStdDev[i][9])
        
        # Scaling factor between event specific and gmpe derived
        sigmaScalingFactor = eventIntraSigma / gmpeIntraSigma
        
        lons.append(float(gmpeStdDev[i][0]))
        lats.append(float(gmpeStdDev[i][1]))
    
        # Where no bias has been applied, we seperate out the other term, then
        # recombine just with gmpe Intra sigma to get total Intra to sample   
        if eventIntraSigma == -1:
        
            otherSigma = sqrt(totalSigma**2 - (gmpeIntraSigma**2 + gmpeInterSigma**2))
            
            # Check as if small -ve will nan   
            if math.isnan(otherSigma):
                otherSigma = 0
                    
            newTotalSigma = sqrt(otherSigma**2 + gmpeIntraSigma**2)          
              
        # where we have a bias and event specific term, we seperate out the shakemap
        # other sigma, then recombine with the event specific intra sigma
        else:
            
            # Scale the total sigma by the scaling factor (direct method)
            newTotalSigma = totalSigma * sigmaScalingFactor            
            
            #otherSigma = sqrt(totalSigma**2 - gmpeIntraSigma**2)
    
            # Check as if small -ve will nan        
            #if math.isnan(otherSigma):
                
                #otherSigma = 0
    
            #newTotalSigma = sqrt(otherSigma**2 + eventIntraSigma**2)

        siteIntraSigmas.append(newTotalSigma)        
        
    return lons, lats, siteIntraSigmas       
        
    
if __name__ == "__main__":
    
    # Number of realisations
    N = 1 
    eventIntraSigma = 0.45
    
    IMT = ['pga']    
    
    stationLons, stationLats = parseStationData('stationlist.xml')
    
    
    eventIntras, gmpeStdDev = parseUncertaintyXML('uncertainty.xml')
    lons1, lats1, gmMedians = parseGridXML('grid.xml')
    
    lons, lats, siteIntraSigmas = getSiteIntraSigma(eventIntras, gmpeStdDev)
 
    lons=np.array(lons)
    lats=np.array(lats)
    
""" 
   means = np.array([[0.16, 0.15,0.19]])
    means= means.T
    
    residuals = eventIntraSigma * sp.stats.norm.rvs(size=(len(lats),N))
    #residuals = np.array([[0.1],[0.0],[0.0]])
    
    ###
    # Do some work
    
 
        
    
    distMatrix = getDistanceMatrix(lons,lats)
    corModel = getCorrelationModel(distMatrix, IMT[0].upper(), True)    
    
    corma = np.linalg.cholesky(corModel)
    
    intraResiduals = np.dot(corma, residuals)
    
    gms = means + intraResiduals
    
    print exp(gms)
    """
    
""" 
we want to have 1D array of lons, lats, mean GMs

    
    # CASE 1:
    # If bias corrected (eventIntra != -1) 
    # then totalSigma = otherSigma + gmpeIntra
    # and we simply need to work out what otherSigma is:
    otherSigma = sqrt(totalSigma - gmpeIntraSigma)  
    
    # if there is a small negative the sqrt will result in a nan so change to 
    # zero for this
    if math.isnan(otherSigma):
        otherSigma = 0:
            
    newTotalSigma = sqrt(otherSigma**2 + eventIntraSigma**2)
    # Then sample this newTotalSigma
    
    # CASE 2:
    # If no bias correction (eventIntra = -1)
    # then total sigma = other + gmpeInter + gmpeIntra
    
    # Calculate the "other component" of the total sigma
    otherSigma = sqrt(totalSigma**2 - (gmpeIntraSigma**2 + gmpeInterSigma**2))
    # if there is a small negative the sqrt will result in a nan so change to 
    # zero for this
    if math.isnan(otherSigma):
        otherSigma = 0:
            
    # Combine with event specific sigma to get new totalSigma that we want to 
    # sample for correlation. This include "shakemap other uncertainty" and
    # GMPE uncertainty
    newTotalSigma = sqrt(otherSigma**2 + eventIntra**2)
    
    # So now we have total sigma to sample to get residuals
    
"""