import numpy as np 
import scipy 
import meep as mp 
from meep import mpb 

import sys 
sys.path.insert(0,'.')
from Materials import * 

### ===========================================================================
### FUNCTION: Calculate the light cone dispersion 
### ===========================================================================
def LightCone(Envir,k_points,resolution):
    num_bands = 10 

    # Define the lattice 
    geometry_lattice = mp.Lattice(
        size = mp.Vector3(1,1,0),
        basis1 = mp.Vector3(1,0),
        basis2 = mp.Vector3(0,1)
    )

    geometry = [
        mp.Block(
            center = mp.Vector3(0,0,0),
            size = mp.Vector3(mp.inf,mp.inf,mp.inf),
            material = Envir 
        )
    ]

    # Define the mode solver 
    ms = mpb.ModeSolver(
        geometry = geometry,
        geometry_lattice = geometry_lattice,
        k_points = k_points,
        resolution = resolution,
        num_bands = num_bands
    )

    # Run the simulation 
    ms.run()

    # Extract the frequency 
    freqs = ms.all_freqs
    lightcone = freqs[:,0]

    # Return the result
    return lightcone 