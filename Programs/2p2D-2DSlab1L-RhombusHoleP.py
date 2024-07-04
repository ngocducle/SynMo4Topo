import numpy as np 
import meep as mp 
from meep import mpb 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab1LRhombusHole 
from DielectricProfile import * 
from BandStructure import * 
from ExportData import * 
from Materials import * 
from LightCone import LightCone 


##### The MAIN program goes here 
def main():
    ##################################################################################
    #                                                                                #
    #           Please define the parameters of the simulation here                  #
    #                                                                                #
    ##################################################################################

    ### The light polarization
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'zeven'
    print('# Polarization = '+polarization)

    ### Resolution 
    resolution = mp.Vector3(8,8,8)  # pixels/a 
    print('# The resolution:'+str(resolution))

    ### Number of bands 
    num_bands = 10 
    print('# The number of bands to simulate: '+str(num_bands))
    
    ### Geometrical parameters 
    h = 0.35        # Thickness of the slab 
    Lz = 5.0        # The height of the unit cell along the z-direction 
    b = 0.35        # The mean of the projections of the diagonal vertices (0 <= b < 0.5)
    
    print('# Thickness of the slab h = '+str(h))
    print('# The edge of the undeformed square hole b = '+str(b))
    print('# The height of the simulation cell Lz = '+str(Lz))

    ### Number of k-points to interpolate between the 2 high-symmetry points 
    Nk = 19 
    print('# The number of points to interpolate the high-symmetry line Nk = '+str(Nk))

    ### The array of synthetic momenta
    N1 = 5  
    k1_array = np.linspace(0.0,0.5,N1)

    N2 = 5 
    k2_array = np.linspace(0.0,0.5,N2)

    # The part of the momentum space to plot the band structure 
    # Choose between: 'BZ', 'Gamma-M', 'M', 'M-vicinity'
    kSpace = 'BZ'
    print('# The part of the momentum space to simulate:'+kSpace)

    ### The set of k-points 
    # kSpace = BZ
    k_points = [
        mp.Vector3(0.0,0.0,0.0),    # Gamma 
        mp.Vector3(0.5,0.0,0.0),    # X 
        mp.Vector3(0.5,0.5,0.0),    # M+
        mp.Vector3(0.0,0.0,0.0),    # Gamma 
        mp.Vector3(-0.5,0.5,0.0)    # M- 
    ]

    # kSpace = M-vicinity 
    #k_points = [
    #    mp.Vector3(0.5,0.5,0.0),    # M 
    #    mp.Vector3(0.45,0.45,0.0),  # Gamma'
    #    mp.Vector3(0.5,0.45,0.0),   # X'
    #    mp.Vector3(0.5,0.5,0.0),    # M 
    #    mp.Vector3(0.55,0.5,0.0)    # Gamma' 
    #] 

    k_points = mp.interpolate(Nk,k_points)

    ### Show the figures (Yes/No)
    show_fig = 'No'
    print('# Show the figures: '+str(show_fig)) 

    ### Materials 
    # The material forming the structure 
    Mater = Si 

    # The environment
    Envir = PMMA 

    ### Calculate the dispersion for the light cone with Envir and k_points 
    lightcone = LightCone(Envir,k_points,resolution)
 
##### Run the MAIN program 
if __name__ == '__main__':
    main()