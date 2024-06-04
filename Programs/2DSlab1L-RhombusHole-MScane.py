import numpy as np 
import meep as mp 
from meep import mpb 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab1LRhombusHole


##### The main program goes here 
def main():
    ##################################################################################
    #                                                                                #
    #          Please define the parameters of the simulation here                   #
    #                                                                                #
    ##################################################################################

    ### The light polarization 
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'zeven'
    print('# Polarization = '+polarization)

    ### Resolution 
    resolution = mp.Vector3(16,16,16)   # pixels/a 

    

if __name__ = '__name__':
    main()