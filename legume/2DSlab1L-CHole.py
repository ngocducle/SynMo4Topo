import numpy as np
import scipy 
import legume
import matplotlib.pyplot as plt

# 2D photonic crystal slab monolayer parameters
a = 1          # period 
h = 0.3        # thickness of the slab
radius = 0.4   # radius of a hole
epss = 3.54**2 # Si 
epse = 1.46**2 # PMMA

# Initialize the lattice and the photonic crystal slab
lattice = legume.Lattice('square')
phc = legume.PhotCryst(lattice,eps_l=epse,eps_u=epse)

# Add the layer
phc.add_layer(d=h,eps_b=epss)

# Add the circular hole
phc.layers[-1].add_shape(legume.Circle(eps=epse,r=radius))

# Define the GME object
gme = legume.GuidedModeExp(phc,gmax=15)

# Define the PWE object
#pwe = legume.PlaneWaveExp(phc,gmax=10)

# Visualize the structure
legume.viz.structure(phc,xz=True,yz=True,figsize=3)

# The path in the Brillouin zone along which we plot the band structure
path = lattice.bz_path(['G','X','M','G'],[39,39,39]) 
gme.run(
        kpoints = path['kpoints'],
        gmode_inds = [1,2,3,4],
        numeig = 10,
        verbose = False
        )

#pwe.run(
#        kpoints = path['kpoints'],
#        pol = 'te'
#        )

#freqs_te = pwe.freqs 

# Plot the figure
fig,ax = plt.subplots(1,figsize=(7,5))
legume.viz.bands(gme,figsize=(5,5),Q=True,ax=ax)
#plt.plot(freqs_te)
ax.set_xticks(path['indexes'])
ax.set_xticklabels(path['labels'])
ax.xaxis.grid('True')

# Show the figure
plt.show()
