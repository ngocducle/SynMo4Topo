import numpy as np
import meep as mp 
from meep import mpb
import matplotlib.pyplot as plt

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
num_bands = 18 

### Number of k-points to interpolate between 2 high-symmetry points
Nk = 9 

##### FUNCTION: ModeSolver 2D photonic crystal slab monolayer with circular holes
def _2DSlab1LCircularHole(h,Lz,radius,num_bands,Nk,resolution,kSpace):
    ### Define the k-points
    if kSpace == 'BZ':  # The whole Brillouin zone
        k_points = [
            mp.Vector3(0.0,0.0,0.0),    # Gamma
            mp.Vector3(0.5,0.0,0.0),    # X
            mp.Vector3(0.5,0.5,0.0),    # M 
            mp.Vector3(0.0,0.0,0.0)     # Gamma 
        ]
    elif kSpace == 'M': # In the vicinity of the M-point
        k_points = [
            mp.Vector3(0.3,0.3,0.0),    # Gamma'
            mp.Vector3(0.5,0.3,0.0),    # X
            mp.Vector3(0.5,0.5,0.0),    # M 
            mp.Vector3(0.3,0.3,0.0)     # Gamma 
        ]
    else:
        print('ERROR! The k-point has not been in the allowed list yet')
        exit()

    # Interpolate to get the points on the lines
    k_points = mp.interpolate(Nk,k_points)

    ### Define the materials
    Si = mp.Medium(index = 3.54)
    SiO2 = mp.Medium(index = 1.46)
    Dielectric = mp.Medium(epsilon = 12.0)
    Air = mp.Medium(epsilon = 1.0)

    Environment = Air 

    ### Define the lattice
    geometry_lattice = mp.Lattice(
        size = mp.Vector3(1.0,1.0,Lz),
        basis1 = mp.Vector3(1.0,0.0),
        basis2 = mp.Vector3(0.0,1.0)
    )

    ### Define the geometry
    geometry = [
        mp.Block(
            center = mp.Vector3(0.0,0.0,0.0),
            size = mp.Vector3(mp.inf,mp.inf,mp.inf),
            material = Environment
        ),

        mp.Block(
            center = mp.Vector3(0.0,0.0,0.0),
            size = mp.Vector3(1.0,1.0,h),
            material = Si
        ),

        mp.Cylinder(
            center = mp.Vector3(0.0,0.0,0.0),
            radius = radius,
            height = h,
            axis = mp.Vector3(0,0,1),
            material = Environment
        )
    ]

    ### The ModeSolver
    ModeSolver = mpb.ModeSolver(
        geometry = geometry,
        geometry_lattice = geometry_lattice,
        k_points = k_points,
        resolution = resolution,
        num_bands = num_bands
    )

    ### Return the mode solver
    return ModeSolver

##### FUNCTION: plot the band structure
def PlotBand_BrillouinZone(number,freqs,namesave):
    fig, ax = plt.subplots()
    ax.plot(number, freqs)
    plt.vlines(Nk+1,0.0,1.0,linestyle='dashed',color='black')
    plt.vlines(2*(Nk+1),0.0,1.0,linestyle='dashed',color='black')
    plt.xlim(0,3*(Nk+1))
    plt.ylim(0,0.5)
    tick_locs = [i*(Nk+1) for i in range(4)]
    tick_labs = [r'$\Gamma$','X','M',r'$\Gamma$']
    ax.set_xticks(tick_locs)
    ax.set_xticklabels(tick_labs,size=16)
    ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14)
    plt.title(namesave,fontsize=14)
    plt.savefig(namesave+'.png')
    plt.show()

##### The MAIN program continues here
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
namesave = '2DSlab1L-CircularHole-h_'+str(h)+'-r_'+str(radius)+'-'+polarization 
print(namesave)

### Plot the band structure
PlotBand_BrillouinZone(number,freqs,namesave)