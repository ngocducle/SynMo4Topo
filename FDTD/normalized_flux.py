import numpy as np 
import scipy 
import meep as mp 

import sys
sys.path.insert(0,'../src')
from Materials import * 

import os 

### Resolution 
resolution = 12 

### PML layer 
dpml = 2.0  # PML thickness 
pml_layers = [mp.PML(dpml)]

### Padding block
pad = 2.0 

### Position of the source and the monitor from the PML layer 
pos = 0.9*pad 

### The materials 
Mater = Si 
Envir = PMMA 
Mater1 = Mater 
Mater2 = Mater 

##### The geometry 
### The number of unit cells along each direction 
Ncellx = 3 
Ncelly = 5 

### The diagonal of one square unit cell 
d = np.sqrt(2.0)

# The layer 1 
h1 = 0.35   # Thickness of the upper layer 
b1 = 0.36   # The edge length of the undeformed square hole
e1 = -0.1    # The anisotropy between the two diagonals of the layer 1 

# The layer 2 
h2 = 0.35   # Thickness of the lower layer
b2 = 0.40   # The edge length of the undeformed square hole 
e2 = -0.1    # The anisotropy between the two diagonals of the layer 2 

# The distance between the 2 layers 
dist = 0.1 

# The total thickness of the bilayer 
htotal = h1 + dist + h2 

# The height of the unit cell along the z-direction 
Lz = 1.5*htotal 

### The size of the unit cell 
sx = 2*(dpml+pad)
sy = Ncelly*d + 2*(dpml+pad)
sz = Lz + 2*dpml 

### The simulation cell 
cell = mp.Vector3(sx,sy,sz)

##### Source 
fcen = 0.30     # pulse center frequency 
df = 0.15       # pulse width 
nfreq = 101     # number of frequencies 
comp = mp.Ez    # the component 

sources = [
    mp.Source(
        mp.GaussianSource(fcen,fwidth=df),
        component = comp,
        center = mp.Vector3(-0.5*sx+dpml+pos,0,0),
        size = mp.Vector3(0,sy,htotal)
    )
]

##### The position of the point we measure the decay rate 
pt = mp.Vector3(0.5*sx-dpml-0.5*pad,0,0)

##### The array of synthetic momenta 
Nq = 6 
q_array = np.linspace(0,0.5,Nq)

##### The geometry 
geometry = [
    mp.Block(
        center = mp.Vector3(0,0,0),
        size = mp.Vector3(mp.inf,mp.inf,mp.inf),
        material = Envir 
    ),

    mp.Block(
        center = mp.Vector3(-0.25*sx,0,0.5*(htotal-h1)),
        size = mp.Vector3(0.5*sx,sy,h1),
        material = Mater1 
    ),

    mp.Block(
        center = mp.Vector3(-0.25*sx,0,0.5*(-htotal+h2)),
        size = mp.Vector3(0.5*sx,sy,h2),
        material = Mater2 
    ),

    mp.Block(
        center = mp.Vector3(0.25*sx,0,0.5*(htotal-h2)),
        size = mp.Vector3(0.5*sx,sy,h2),
        material = Mater2 
    ),

    mp.Block(
        center = mp.Vector3(0.25*sx,0,0.5*(-htotal+h1)),
        size = mp.Vector3(0.5*sx,sy,h1),
        material = Mater1 
    )
]

##### The simulation 
sim = mp.Simulation(cell_size = cell,
                    geometry = geometry,
                    sources = sources,
                    boundary_layers = pml_layers,
                    resolution = resolution)

##### Transmitted flux 
tran_fr = mp.FluxRegion(
    center = mp.Vector3(0.5*sx-dpml-pos,0,0),
    size = mp.Vector3(0,sy,htotal)
)

tran = sim.add_flux(fcen,df,nfreq,tran_fr)

##### Run the simulation 
sim.run(
    until_after_sources = mp.stop_when_fields_decayed(
        1000,comp,pt,1e-4
    )
)

##### The transmission flux 
tran_norm_flux = np.array(mp.get_fluxes(tran))

##### The array of data 
#q_array_print = q*np.ones(nfreq)
freq_array = np.linspace(fcen-df,fcen+df,nfreq)

dataexport = np.column_stack((freq_array,tran_norm_flux))

##### Write the transmission flux to file 
with open("norm_tran_flux.txt","w") as file:
    np.savetxt(file,dataexport,fmt='%.8f') 