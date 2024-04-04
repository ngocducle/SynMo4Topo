import numpy as np
import meep as mp 
from meep import mpb
import matplotlib.pyplot as plt

# Resolution 
resolution = mp.Vector3(16,16,16)   # pixels/a

# Geometrical parameters
h = 0.40     # Thickness of the slab
radius = 0.4 # Radius of the hole
Lz = 5.0     # The height of the unit cell along the z-direction

# Number of bands
num_bands = 18 

# Number of k-points to interpolate between 2 high-symmetry points
Nk = 19 

