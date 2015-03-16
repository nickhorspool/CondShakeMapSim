# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 07:51:39 2015

@author: horspool
"""

import shakeMapUtils
import smtk.hazard.conditional_simulation as csim
from openquake.hazardlib.site import Site, SiteCollection
from openquake.hazardlib.geo.point import Point
import numpy as np
import matplotlib.pylab as plt
import os
import configCSS

print "Initialsing CSS" 

nSims = configCSS.nSims

print "Number of simulations:", nSims

# Input files

# ShakeMap grid file
shakeMapGrid = configCSS.scenarioDir + '/output/grid.xml'

# ShakeMap uncertainty file
shakeMapUncertaintyGrid = configCSS.scenarioDir + '/output/uncertainty.xml'

# ShakeMap station list file
shakeMapStations = configCSS.scenarioDir + '/output/stationlist.xml'

# Get event specific intraevent uncertainty
shakeMapUncertaintyDict, gmStdDevs = shakeMapUtils.parseUncertaintyXML(shakeMapUncertaintyGrid)
shakeMapUncertainty = shakeMapUncertaintyDict['pga']

# Get median ground motions from shakemap and grid locations and turn to OQ site collection
modLons, modLats, gmMedians = shakeMapUtils.parseGridXML(shakeMapGrid)

nSitesMod=len(modLons)

modelledSites = SiteCollection([
        Site(Point(modLons[i], modLats[i]), configCSS.vs30, True, configCSS.z1pt0, configCSS.z2pt5, i)
        for i in range(0, nSitesMod)])
           
           
# Get observed sites from shakemap and convert to OQ site collection
obsLons, obsLats = shakeMapUtils.parseStationData(shakeMapStations)  

nSitesObs=len(obsLons)

observedSites = SiteCollection([
        Site(Point(obsLons[i], obsLats[i]), configCSS.vs30, True, configCSS.z1pt0, configCSS.z2pt5, i)
        for i in range(0, nSitesObs)])
            
obsResiduals = np.zeros_like(obsLons)
          
            

# Do conditional simulation of residuals
# PGA for correlation model, 1 = nsim, sites are site colletion
# needt modify this function to take in shakemap stddev instead of sampling (0,1) normal dist. so
# we will change the 1 to a value from shakemap
# these are just correlated number of std devs (epsilon) - need to add intraevent uncertainty to this
# here we actually generate a multdimensional array with N simulations. 
# this takes most of hte time.

print "Starting Simulations"

epsilon = csim.conditional_simulation(observedSites, obsResiduals, modelledSites, "PGA", nSims)


####################################
# Loop over results from each simulation and put in directory
#

outputDir = configCSS.scenarioDir + '/' + configCSS.outputDirName

os.system('mkdir %s' % outputDir)

# Write some common files
stationsFile = '%s/stations.xyz' % outputDir
fS = open(stationsFile,'w')
fS.write('lon lat\n')

for i in range(0, len(observedSites.lons)):
    fS.write('%s %s\n' % (observedSites.lons[i], observedSites.lats[i]))
fS.close()

mediansFile = '%s/medians.xyz' % outputDir
fM = open(mediansFile,'w')
fM.write('lon lat GM\n')

for i in range(0, len(gmMedians)):
    fM.write('%s %s %s\n' % (modelledSites.lons[i], modelledSites.lats[i],gmMedians[i]))  
fM.close()

for sim in range(0, nSims):
    
    # epsilon use site specific uncertainty from shakemap here 
    # change the 0 to N here for N simulations.
    simResidualsLog = epsilon[:, sim].A.flatten() * shakeMapUncertainty
    
    # Add these to Shakemap median values
    simGMFLog = np.log(gmMedians) + simResidualsLog
    
    simGMF = np.exp(simGMFLog)
    
    #######################################################################
    # Write simulated GMF to output file
    
    print "Writing Output Files", sim, "of", nSims
    
    outputGridFile = '%s/grid%s.xyz' % (outputDir, sim)
    f = open(outputGridFile,'w' ) 
    f.write('lon lat GM\n')
    
    for i in range(0, len(simGMF)):
        f.write('%s %s %s\n' % (modelledSites.lons[i], modelledSites.lats[i],simGMF[i]))  
    f.close()
    


    if configCSS.generatePlots == True:
    
        ################################################################
        # Do some plotting
        
        print "Plotting", sim, "of", nSims
        
        minLon = min(modelledSites.lons)
        maxLon = max(modelledSites.lons)
        
        minLat = min(modelledSites.lats)
        maxLat = max(modelledSites.lats)
        
        R = '-R%s/%s/%s/%s' % (minLon, maxLon, minLat, maxLat)
        J = '-JM14c'
        plot = '%s/GMF_%s.ps' % (outputDir, sim)
        grid = '%s/GMF_%s.grd' % (outputDir, sim)
        xyzFile = '%s/grid%s.xyz' % (outputDir, sim)
        xyzFileMedians = '%s/medians.xyz' % outputDir
        stationsFile = '%s/stations.xyz' % outputDir
        gridMedians = '%s/GMF_medians.grd' % outputDir
        plotMedians = '%s/GMF_Medians.ps' % outputDir
        
        # FIXME - hardcoded
        I = '0.025'
        
        ############################################
        # GMT 5.1.1 plotting
        os.system('makecpt -Cplotting/Ii.cpt -T0/0.5 -Z > plotting/shake.cpt')
        os.system('xyz2grd %s %s -I%s -h1 -G%s' % (xyzFile, R, I, grid))
        #os.system('grdsample %s -Gtemp.grd -I0.005 -T' % grid)
        os.system('grdimage %s %s %s -Cplotting/shake.cpt -P -K -Y4c > %s' % (grid, R, J, plot))
        os.system('psscale -D7/-1/10/0.5h -Cplotting/shake.cpt -Ba0.1f0.1:"PGA (g)": -O -K -P >> %s' % plot)
        os.system('pscoast %s %s -Dh -Slightblue -Ba1f0.5 -W0.1p -P -O -K >> %s' % (R, J, plot))
        os.system('psxy %s %s %s -St0.4c -Gred -W1p -O -P >> %s' % (R,J,stationsFile,plot))
        os.system('ps2raster -Tg -A %s' % plot)
        
        #os.system('makecpt -CIi.cpt -T0/0.5 -Z > shake.cpt')
        os.system('xyz2grd %s %s -I%s -h1 -G%s' % (xyzFileMedians, R, I, gridMedians))
        #os.system('grdsample %s -Gtemp.grd -I0.005 -T' % grid)
        os.system('grdimage %s %s %s -Cplotting/shake.cpt -P -K -Y4c > %s' % (gridMedians, R, J, plotMedians))
        os.system('psscale -D7/-1/10/0.5h -Cplotting/shake.cpt -Ba0.1f0.1:"PGA (g)": -O -K -P >> %s' % plotMedians)
        os.system('pscoast %s %s -Dh -Slightblue -Ba1f0.5 -W0.1p -P -O -K >> %s' % (R, J, plotMedians))
        os.system('psxy %s %s %s -St0.4c -Gred -W1p -O -P >> %s' % (R,J,stationsFile,plotMedians))
        os.system('ps2raster -Tg -A %s' % plotMedians)
        
        
        #########################################
        # python plotting. Plot epsilon (number of std devs)
        plt.figure()
        plt.scatter(modelledSites.lons, modelledSites.lats, 
                    s=10,
                    marker="s",
                    edgecolor="None",
                    #c=output_resid[:, 0].A.flatten(), 
                    c=epsilon[:, sim].A.flatten())
                    #norm=Normalize(vmin=-3.0, vmax=3.0))
        cbar = plt.colorbar()
        cbar.set_label('Epsilon (No. of Std Devs)', rotation=90)
        plt.scatter(observedSites.lons, observedSites.lats, s=40, c='red')
                    #norm=Normalize(vmin=-3.0, vmax=3.0))
        plt.xlim(min(modelledSites.lons),max(modelledSites.lons))
        plt.ylim(min(modelledSites.lats),max(modelledSites.lats))
        plt.savefig('%s/simEpsilon.png' % outputDir)
        
        # Plot epsilon (number of std devs)
        plt.figure()
        plt.scatter(modelledSites.lons, modelledSites.lats, 
                    s=10,
                    marker="s",
                    edgecolor="None",
                    #c=output_resid[:, 0].A.flatten(), 
                    c=simResidualsLog)
                    #norm=Normalize(vmin=-3.0, vmax=3.0))
        cbar = plt.colorbar()
        cbar.set_label('Log Residuals (simulated)', rotation=90)
        plt.scatter(observedSites.lons, observedSites.lats, s=40, c='red')
                    #norm=Normalize(vmin=-3.0, vmax=3.0))
        plt.xlim(min(modelledSites.lons),max(modelledSites.lons))
        plt.ylim(min(modelledSites.lats),max(modelledSites.lats))
        plt.savefig('%s/simLogresiduals.png' % outputDir)
        
        # Plot simulated ground motion field
        plt.figure()
        plt.scatter(modelledSites.lons, modelledSites.lats, 
                    s=10,
                    marker="s",
                    edgecolor="None",
                    #c=output_resid[:, 0].A.flatten(), 
                    c=simGMF, vmin=0, vmax=0.2)
                    #c=simResiduals[:, 0].A.flatten())
                    #norm=Normalize(vmin=0, vmax=2.0))
        cbar = plt.colorbar()
        cbar.set_label('PGA (g)', rotation=90)
        plt.scatter(observedSites.lons, observedSites.lats, s=40, c='red')
                    #norm=Normalize(vmin=-3.0, vmax=3.0))
        plt.xlim(min(modelledSites.lons),max(modelledSites.lons))
        plt.ylim(min(modelledSites.lats),max(modelledSites.lats))
        plt.savefig('%s/simulatedGMF_%s.png' % (outputDir,sim))
        
        # Plot median ground motions
        plt.figure()
        plt.scatter(modelledSites.lons, modelledSites.lats, 
                    s=10,
                    marker="s",
                    edgecolor="None",
                    #c=output_resid[:, 0].A.flatten(), 
                    c=gmMedians, vmin=0, vmax=0.2)
                    #norm=Normalize(vmin=0, vmax=2.0))
        cbar = plt.colorbar()
        cbar.set_label('PGA (g)', rotation=90)
        plt.scatter(observedSites.lons, observedSites.lats, s=40, c='red')
                    #norm=Normalize(vmin=-3.0, vmax=3.0))
        plt.xlim(min(modelledSites.lons),max(modelledSites.lons))
        plt.ylim(min(modelledSites.lats),max(modelledSites.lats))
        plt.savefig('%s/medianGMF.png' % outputDir)
