import numpy as np 
import meep as mp 
from meep import mpb 

import sys 
sys.path.insert(0,'../src/')


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

    

##### Run the MAIN program 
if __name__ == '__main__':
    main()