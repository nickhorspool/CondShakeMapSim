# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 09:19:01 2015

@author: horspool

Configuration file for running condShakeMapSim.py

"""

# Number of simulations (realizations)
nSims = 3

# scenario directory
scenarioDir = 'data/2015p012816'

# output directory name
outputDirName = 'condSimOutputs'

# do some plotting
generatePlots = True

# Site conditions. Note these are just placeholders for the OQ siteclass 
# object and not used in the calculations
vs30=760.
z1pt0=100.
z2pt5=1.