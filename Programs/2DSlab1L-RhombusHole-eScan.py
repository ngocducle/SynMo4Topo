import numpy as np 
import meep as mp 
from meep import mpb 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab1LRhombusHole
from Materials import * 
from LightCone import LightCone 

import os 

##### The MAIN program goes here 
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

    ### Resolution 
    resolution = mp.Vector3(8,8,8)   # pixels/a 
    print('# The resolution:'+str(resolution))

    ### Number of bands 
    num_bands = 40 
    print('# The number of bands to simulate: '+str(num_bands))

    ### The thickness of the slab
    h = 0.35 

    ### The edge length of the square hole when the two diagonals are equal 
    b = 0.4  

    ### The height of the unit cell along the z-direction 
    Lz = 5.0

    ### The array of the anisotropy between the two diagonals (-1 <= e <= 1)
    Ne = 6
    e_array = np.linspace(-0.1,0.1,Ne)

    ### Show the figure (Yes/No)
    show_fig = 'No'
    print('# Show the figure: '+str(show_fig))

    ### Materials 
    # The material forming the structure 
    Mater = Si 

    # The environment
    Envir = PMMA 

    ### Number of k-points to interpolate between 2 high-symmetry points 
    Nk = 29 
    print('# The number of points to interpolate the high-symmetry line Nk = '+str(Nk))

    # The set of k-points 
    k_points = [
        mp.Vector3(0.0,0.0,0.0),    # Gamma 
        mp.Vector3(0.5,0.0,0.0),    # X 
        mp.Vector3(0.5,0.5,0.0),    # M 
        mp.Vector3(0.0,0.0,0.0),    # Gamma 
        mp.Vector3(0.0,0.5,0.0),    # Y 
        mp.Vector3(0.5,0.5,0.0)     # M 
    ]

    k_points = mp.interpolate(Nk,k_points)

    ### Calculate the dispersion for the light cone with Envir and k_points 
    lightcone = LightCone(Envir,k_points,resolution)

    ############################################################################
    #                                                                          #
    #                       Here the simulation starts                         #
    #                                                                          #
    ############################################################################
    for e in e_array:
        ### Geometrical parameters 
        print('# Thickness of the slab h = '+str(h))
        print('# The edge length of the parent square is b = '+str(b))
        print('# The anisotropy between the two diagonals e = '+str(e))
        print('# The height of the simulation cell Lz = '+str(Lz))

        ### The vertices of the rhombus hole
        # If e>0 then the diagonal x=y is longer than the diagonal x=-y
        # If e=0 then the hole is a square 
        # If e<0 then the diagonal x=y is shorter than the diagonal x=-y 
        vertices = [
            mp.Vector3(  0.5*(1+e)*b,  0.5*(1+e)*b, 0.0 ),
            mp.Vector3(  0.5*(1-e)*b, -0.5*(1-e)*b, 0.0 ),
            mp.Vector3( -0.5*(1+e)*b, -0.5*(1+e)*b, 0.0 ),
            mp.Vector3( -0.5*(1-e)*b,  0.5*(1-e)*b, 0.0 )
        ]

        ### Define the mode solver 
        ms = _2DSlab1LRhombusHole(h,Lz,b,e,vertices,
                                  num_bands,k_points,resolution,Mater,Envir)
        
        ### Run the simulation 
        if polarization == 'all':
            ms.run()
        elif polarization == 'zeven':
            ms.run_zeven()
        elif polarization == 'zodd':
            ms.run_zodd()
        else: 
            print('ERROR! The polarization does not belong to the allowed list')
            exit()

        ### Extract the frequencies of the modes from the ModeSolver
        freqs = ms.all_freqs 

        ### The number of elements in k_points 
        number = np.arange(len(ms.k_points))

        

##### Run the MAIN program 
if __name__ == "__main__":
    main()