import numpy as np 
import meep as mp 
from meep import mpb 

import sys 
sys.path.insert(0,'../src/')
from Materials import * 

##### The main program goes here 
def main():
    ##################################################################################
    #                                                                                # 
    #             Please define the parameters of the simulation here                #
    #                                                                                #
    ##################################################################################

    ### The light polarization 
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'zeven'
    print('# Polarization = '+polarization)

    ### The part of the momentum space to plot the band structure 
    # Choose between: 'BZ', 'Gamma-M', 'M', 'M-vicinity'
    #kSpace = 'BZ'
    #print('# The part of the momentum space to simulate:'+kSpace) 

    ### The k-point at which we plot the field profile 
    k_field = mp.Vector3(0.5,0.5,0.0)   # M-point 
    print('# The k-point at which we plot the field profile:'+str(k_field))

    ### Resolution 
    resolution = mp.Vector3(8,8,8)  # pixels/a 
    print('# The resolution:'+str(resolution))

    ### Geometrical parameters 
    h = 0.35        # Thickness of the slab 
    Lz = 5.0       # The height of the unit cell along the z-direction 
    c = 0.3    # The mean of the projections of the diagonals vertices (0<= mean_d <0.5)
    alpha_d = 0.2   # The anisotropy between the two diagonals (-1<=alpha_d<=1)

    print('# Thickness of the slab h = '+str(h))
    print('# Mean half-axis of the rhombus c = '+str(c))
    print('# The height of the simulation cell Lz = '+str(Lz))

    ### The vertices of the rhombus hole 
    # If alpha > 0 then the diagonal x=y is longer than the diagonal x=-y 
    # If alpha = 0 then the hole is a square 
    # If alpha < 0 then the diagonal x=y is shorter than the diagonal x=-y 
    vertices = [
        mp.Vector3( (1+alpha_d)*c,  (1+alpha_d)*c, 0.0 ),
        mp.Vector3( (1-alpha_d)*c, -(1-alpha_d)*c, 0.0 ),
        mp.Vector3(-(1+alpha_d)*c, -(1+alpha_d)*c, 0.0 ),
        mp.Vector3(-(1-alpha_d)*c,  (1-alpha_d)*c, 0.0 )
    ]

    ### Number of bands 
    num_bands = 20 
    print('# The number of bands to simulateL '+str(num_bands))

    ### Number of k-points to interpolate between the 2 high-symmetry points 
    Nk = 29 
    print('# The number of points to interpolate the high-symmetry line Nk = '+str(Nk))

    ### Show figure (Yes/No)
    show_fig = 'Yes'
    print('# Show the figure: '+str(show_fig))

    ### Materials 
    # The material forming the structure 
    Mater = Si 

    # The envionment 
    Envir = PMMA 





##### Run the MAIN program 
if __name__ == "__main__":
    main() 