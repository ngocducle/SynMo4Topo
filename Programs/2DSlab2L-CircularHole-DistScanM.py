import numpy as np 
import meep as mp 
from meep import mpb 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab2LCircularHole
from Materials import * 

##############################################################################
#                                                                            #
#    The program scans the 2DSlab2L with varying interlayer distance         #
#    and calculate the energy eigenvalues at the M-point                     # 
#    The shifts deltax = deltay = 0                                          #
#                                                                            #
##############################################################################


##### The MAIN program goes here 
def main():
    ##########################################################################
    #                                                                        #
    #        Please define the parameters of the simulation here             #
    #                                                                        #
    ##########################################################################
    
    ### The light polarization 
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'zeven'
    print('# Polarization = '+polarization)

    ### We focus on the M-point 
    kSpace = 'M'
    print('# We focus on the M-point')

    ### Resolution 
    resolution = mp.Vector3(32,32,32)   # pixels/a 
    print('# The resolution:'+str(resolution))

    ### Geometrical parameters 
    h = 0.3         # Thickness of the slab 
    radius = 0.4    # Radius of the hole 
    Lz = 10.0       # The height of the unit cell along the z-direction 

    print('# Thickness of the slab h = '+str(h))
    print('# Radius of the circular hole r = '+str(radius))
    print('# the height of the simulation cell Lz = '+str(Lz))

    ### Number of bands 
    num_bands = 30 
    print('# The number of bands to simulate: '+str(num_bands))

    ### Show figure (Yes/No)
    show_fig = 'No'
    print('# Show the figures: '+str(show_fig))

    ### Materials
    # The material forming the structure 
    Mater = Si 

    # The environment
    Envir = PMMA 

    ##########################################################################
    #                                                                        #
    #                   Here the simulation starts                           #
    #                                                                        #
    ##########################################################################

    ### Define the array for the interlayer distance
    Ndist = 31 
    DistArray = np.linspace(0.0,3.0,Ndist)

##### Run the MAIN program 
if __name__ == "__main__":
    main()
