import numpy as np 
import meep as mp 
from meep import mpb

import sys 
sys.path.insert(0,'../src/')
from DielectricProfile import Si,PMMA

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

    ### The part of the momentum spsace to plot the band structure
    # Choose between: 'BZ', 'M-vicinity'
    kSpace = 'BZ'
    print('# The part of the momentum space to simulate:'+kSpace)

    ### The k-point at which we plot the field profile
    k_field = mp.Vector3(0.5,0.5,0.0)   # M-point
    print('# The k-point at which we plot the field profile:'+str(k_field))

    ### Resolution
    resolution = mp.Vector3(8,8,8)  # pixels/a 
    print('# The resolution:'+str(resolution))

    ### Geometrical parameters
    h = 0.3     # Thickness of the slab
    b = 0.9     # Edge length of the square hole
    Lz = 5.0    # The height of the unit cell along the z-direction

    print('# Thickness of the slab h = '+str(h))
    print('# Edge length of the square hole b = '+str(b))
    print('# The height of the simulation cell Lz = '+str(Lz))

    ### Number of bands 
    num_bands = 20 
    print('# The number of bands to simulate: '+str(num_bands))

    ### Number of k-points to interpolate between 2 high-symmetry points
    Nk = 19 
    print('# The number of points to interpolate the high-symmetry line Nk = '+str(Nk))

    ### Show figure (Yes/No)
    show_fig = 'No'
    print('# Show the figure: '+str(show_fig))

    ### Materials 
    # The material forming the structure 
    Mater = Si 

    # The environment
    Envir = PMMA

    ### The parameters to calculate the E- and H-field profiles 
    ### in the planes parallel to Oxy
    Bloch_Phase = 'True'    # Enable the Bloch phase: True/False
    resolution_eps = 81     # Number of pixels per a 
    resolution_field = 81   # Number of pixels per a 
    num_periods = 3         # Number of periods along each direction 

    # The value of z where we take a slice to plot the field 
    zvalue = 0.075 



##### Run the MAIN program
if __name__ == "__main__":
    main()