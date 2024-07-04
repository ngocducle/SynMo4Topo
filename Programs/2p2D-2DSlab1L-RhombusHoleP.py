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


##### Run the MAIN program 
if __name__ == '__main__':
    main()