import numpy as np
import meep as mp 
from meep import mpb

import sys
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab1LCircularHole
from PlotBandStructure import *
from ExportData import PrintBandStructure

### The light polarization
# Choose between: 'all', 'zeven', 'zodd' 
polarization = 'all'

### The part of the momentum space to plot the band structure
# Choose between: 'BZ', 'M' 
kSpace = 'BZ'

### Resolution 
resolution = mp.Vector3(8,8,8)   # pixels/a

### Geometrical parameters
h = 0.40     # Thickness of the slab
radius = 0.4 # Radius of the hole
Lz = 5.0     # The height of the unit cell along the z-direction

### Number of bands
num_bands = 8 

### Number of k-points to interpolate between 2 high-symmetry points
Nk = 9 

### Define the mode solver
ms = _2DSlab1LCircularHole(h,Lz,radius,num_bands,Nk,resolution,kSpace)

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

### The title and the name of the files 
namesave = '2DSlab1L-CircularHole-h_'+str(h)+'-r_'+str(radius) \
    +'-'+polarization+'-'+kSpace 

### Print the band structure to file
PrintBandStructure(freqs,number,namesave)

### Plot the band structure
if kSpace == 'BZ':
    PlotBand_BrillouinZone(number,freqs,Nk,namesave)
elif kSpace == 'M':
    PlotBand_M(number,freqs,Nk,namesave)
else:
    print('ERROR! The k-point has not been in the allowed list yet')
    exit()