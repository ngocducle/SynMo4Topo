import numpy as np
import meep as mp 
from meep import mpb
import matplotlib.pyplot as plt

### ====================================================================
### Function: calculate the band structure and parity for 2D PhC slab bilayer
### with circular hole at the center

def Parity_2DPhCSlabMonolayerCircularHole_M(h,radius,dist,num_bands):
    # Resolution 
    resolution = mp.Vector3(16,16,16)

    # The height of a unit
    Lz = 5.0

    # Define the materials
    Si = mp.Medium(index = 3.54)
    SiO2 = mp.Medium(index = 1.46)
    PMMA = mp.Medium(index = 1.46)
    Dielectric = mp.Medium(epsilon = 12.0)
    Air = mp.Medium(epsilon = 1.0)
    
    Environment = PMMA 

    # Define the lattice 
    geometry_lattice = mp.Lattice(size = mp.Vector3(1.0,1.0,Lz),
                                  basis1 = mp.Vector3(1.0,0.0),
                                  basis2 = mp.Vector3(0.0,1.0))
    
    # Define the geometry 
    geometry = [ mp.Block(center = mp.Vector3(0.0,0.0,0.0), 
                          size = mp.Vector3(mp.inf,mp.inf,mp.inf), 
                          material = Environment),

                 mp.Block(center = mp.Vector3(0.0,0.0,dist), 
                          size = mp.Vector3(1.0,1.0,h), 
                          material = Si),

                 mp.Cylinder(center = mp.Vector3(0.0,0.0,dist), 
                             radius = radius, 
                             height = h, 
                             axis = mp.Vector3(0.0,0.0,1.0), 
                             material = Environment)
                ] 

    # Define the k-point
    k_points = [mp.Vector3(0.5,0.5,0.0)]

    # Run the simulation 
    ms = mpb.ModeSolver(geometry = geometry,
                        geometry_lattice = geometry_lattice, 
                        k_points = k_points,
                        resolution = resolution,
                        num_bands = num_bands 
                        )
    
    ms.run()

    freqs = ms.all_freqs 

    parity = ms.compute_yparities()

    # Return the results
    return freqs, parity

### =================================================================================
### Run the simulation

# Thickness of one layer
h = 0.3

# Radius of a hole
radius = 0.4

# Number of bands 
num_bands = 10  

# Array of distance between the two layers
Ndist = 21
Dist = np.linspace(0.0,2.0,Ndist)

### ==================================================================================
### Scan over the distarray and find the frequencies and parities of the bands
Bands = np.zeros((Ndist,num_bands))
Parity = np.zeros((Ndist,num_bands))

print('# Shape of Parity_array:'+str(np.shape(Parity)))

for i in np.arange(Ndist):
    dist = Dist[i]
    print('# Distance = '+str(dist))

    freqs,parity = Parity_2DPhCSlabMonolayerCircularHole_M(h,radius,dist,num_bands)

    Bands[i,:] = freqs 
    Parity[i,:] = parity 

### =================================================================================
### Print the frequencies and parities of the modes into a file
### Column 0: distarray
### Column j (1 <= j <= num_bands): Band j 

    
with open('2Dbilayer-freqs.txt','w') as file:
    for i in range(Ndist):
        file.write('%.8f        ' % Dist[i])

        for j in range(num_bands):
            file.writelines('%.8f        ' % Bands[i,j])

        file.write('\n') 

with open('2Dbilayer-parity.txt','w') as file:
    for i in range(Ndist):
        file.write('%.8f        ' % Dist[i])

        for j in range(num_bands):
            file.writelines('%.8f       ' % Parity[i,j])

        file.write('\n')

### ==================================================================================
### Plot the figure of all the bands
#fig, ax = plt.subplots(figsize=(12,12))
#plt.scatter(Dist,Bands[:,0],s=5,c=Parity[:,0],cmap='rainbow')
#plt.scatter(Dist,Bands[:,1],s=5,c=Parity[:,1],cmap='rainbow')
#plt.scatter(Dist,Bands[:,2],s=5,c=Parity[:,2],cmap='rainbow')
#plt.scatter(Dist,Bands[:,3],s=5,c=Parity[:,3],cmap='rainbow')
#plt.scatter(Dist,Bands[:,4],s=5,c=Parity[:,4],cmap='rainbow')
#plt.scatter(Dist,Bands[:,5],s=5,c=Parity[:,5],cmap='rainbow')
#plt.scatter(Dist,Bands[:,6],s=5,c=Parity[:,6],cmap='rainbow')
#plt.scatter(Dist,Bands[:,7],s=5,c=Parity[:,7],cmap='rainbow')
#cbar = plt.colorbar()
#cbar.set_ticks(ticks=[-1.0,1.0],labels=['Odd','Even'],fontsize=14)
#cbar.set_label('Parity',fontsize=14)
#plt.xticks(fontsize=14)
#plt.yticks(fontsize=14)
#plt.xlabel('d/a',fontsize=14)
#plt.ylabel(r'$\omega a/(2\pi c)$',fontsize=14)
#plt.savefig('PhCMonolayer-AllBands-Parity.png')
#plt.show()