import numpy as np
import matplotlib.pyplot as plt

### ============================================================================== ###
###                                                                                ###
### This module contains all the functions to plot the dielectric profile          ###
###                                                                                ###
### ============================================================================== ###

def PlotDielectricProfileZ(x,y,zmin,zmax,Nz):
    zarray = np.linspace(zmin,zmax,Nz)


