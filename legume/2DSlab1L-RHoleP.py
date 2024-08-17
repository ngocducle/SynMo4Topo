import numpy as np
import scipy 
import legume
import matplotlib.pyplot as plt

# 2D photonic crystal slab monolayer parameters
a = 1           # period
h = 0.35        # thickness of the slab
b = 0.38        # edge of a undeformed hole
e = -0.1        # deformation factor 
f = (1+e)/(1-e) 
epss = 3.54**2  # Si 
epse = 1.46**2  # PMMA 

# Initialize the lattice and the photonic crystal slab
lattice = legume.Lattice('square')
phc = legume.PhotCryst(lattice,eps_l=epse,eps_u=epse)

# Add the layer 
phc.add_layer(d=h,eps_b=epss)

# Add the rhombus hole
#phc.layers[-1].add_shape(legume.Circle(eps=epse,r=radius))
rhombus = legume.Poly(x_edges = [0.5*b*f,-0.5*b/f,-0.5*b*f,0.5*b/f],
                      y_edges = [0.5*b*f,0.5*b/f,-0.5*b*f,-0.5*b/f],
                      eps = epse)

phc.layers[-1].add_shape(rhombus)

# Define the GME object
gme = legume.GuidedModeExp(phc,gmax=4)

# Visualize the structure 
legume.viz.structure(phc,xz=True,yz=True,figsize=3)

# The path in the Brillouin zone along which we plot the band structure
G1 = [0.45*np.pi,0.45*np.pi]
X1 = [0.5*np.pi,0.45*np.pi]
#G2 = [0.55*np.pi,0.45*np.pi]

path = lattice.bz_path(['G','X','M','G'],[39,39,39])
gme.run(
    kpoints = path['kpoints'],
    gmode_inds = [0,2,4,6],
    numeig = 10,
    verbose = False
)

# The file name 
namesave = '2DSlab1L-RHoleP-h_'+str(h)+'-b_'+str(b)+'-e_'+str(e)

# Plot the figure
fig,ax = plt.subplots(1,figsize=(7,5))
legume.viz.bands(gme,figsize=(5,5),Q=True,ax=ax)
#plt.plot(freqs_te)
ax.set_xticks(path['indexes'])
ax.set_xticklabels(path['labels'])
ax.xaxis.grid('True')
plt.savefig(namesave+'.png')

plt.show()
