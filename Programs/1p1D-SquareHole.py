import numpy as np 
import scipy 
import meep as mp 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab2LSquareHole 
from DielectricProfile import * 
from BandStructure import * 
from Materials import * 

##### The MAIN program goes here
def main():
    ##########################################################################
    #                                                                        #
    #        Please define the parameters of the simulation here             #
    #                                                                        #
    ##########################################################################

    ### The light polarization 
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'all'
    print('# Polarization = '+polarization)

    ### The part of the momentum space to plot the band structure
    # Choose between: 'BZ', 'Gamma-M', 'M', 'M-vicinity'
    kSpace = 'Gamma-M' 
    print('# The part of the momentum space to simulate:'+kSpace) 

    ### Resolution 
    resolution = mp.Vector3(16,16,16)   # pixels/a 
    print('# The resolution:'+str(resolution))

##### Run the MAIN program 
if __name__ == "__main__":
    main()