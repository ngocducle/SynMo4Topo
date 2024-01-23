import numpy as np
import scipy 
import meep as mp
from meep import mpb 
import matplotlib.pyplot as plt 

# Resolution 
resolution = 32 # pixels / a

# Number of bands 
num_bands = 5   

# Number of k-points to interpolate between 2 high symmetry points
N_k = 19 

# Define the materials 
Air  = mp.Medium(epsilon = 1.0)
Si   = mp.Medium(index = 3.54)
SiO2 = mp.Medium(index = 1.46)

# Define the lattice 
# Set the size Lz of the cell to be 0 because the structure 
# is translationally invariant along the z-direction 
geometry_lattice = mp.Lattice(size = mp.Vector3(1,1,0), 
                              basis1 = mp.Vector3(1,0),
                              basis2 = mp.Vector3(0,1)    
                             )

geometry = [ mp.Block(center = mp.Vector3(0,0,0),
                      size   = mp.Vector3(mp.inf, mp.inf, mp.inf),
                      material = SiO2)
           ]

# Define the k-points 
k_points = [mp.Vector3(0,0,0),      # Gamma 
            mp.Vector3(0.5,0,0),    # X 
            mp.Vector3(0.5,0.5,0),  # M 
            mp.Vector3(0,0,0),      # Gamma 
           ]

k_points = mp.interpolate(N_k,k_points) 

# Run the simulation 
ms = mpb.ModeSolver(geometry   = geometry,
                    geometry_lattice = geometry_lattice, 
                    k_points   = k_points,
                    resolution = resolution,
                    num_bands  = num_bands 
                    )

ms.run_te()
te_freqs = ms.all_freqs 

ms.run_tm() 
tm_freqs = ms.all_freqs 

ms.run_zeven()
zeven_freqs = ms.all_freqs 

ms.run_zodd()
zodd_freqs  = ms.all_freqs 

# Print the photonic band of band TE to a file
with open('SiO2-TE.txt','w') as file:

    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in te_freqs[n])
        file.write('\n') 

number = np.arange(len(k_points)) 

# Plot the TE bands 
fig, ax = plt.subplots() 
ax.plot(number, te_freqs)      
#plt.vlines(0,0,0.2,linestyle = 'dashed',color='black')
plt.vlines(N_k+1,0,1.0,linestyle = 'dashed',color='black') 
plt.vlines(2*(N_k+1),0,1.0,linestyle = 'dashed',color='black') 
#plt.vlines(3*(N_k+1),0,0.2,linestyle = 'dashed',color='black') 
plt.xlim(0,3*(N_k+1)) 
plt.ylim(0,0.8)   
#tick_locs = [ 0, N_k+1, 2*(N_k+1), 3*(N_k+1)]   
tick_locs = [i*(N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$',fontsize = 14) 
plt.title('TE band',fontsize = 14) 
plt.savefig('SiO2-TE.png') 
plt.show() 

# Print the photonic band of band TM to a file
with open('SiO2-TM.txt','w') as file:

    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in tm_freqs[n])
        file.write('\n') 

# Plot the TM bands 
fig, ax = plt.subplots() 
ax.plot(number, tm_freqs)       
#plt.vlines(0,0,0.2,linestyle = 'dashed',color='black')
plt.vlines(N_k+1,0,1.0,linestyle = 'dashed',color='black') 
plt.vlines(2*(N_k+1),0,1.0,linestyle = 'dashed',color='black') 
#plt.vlines(3*(N_k+1),0,0.2,linestyle = 'dashed',color='black') 
plt.xlim(0,3*(N_k+1)) 
plt.ylim(0,0.8)   
#tick_locs = [ 0, N_k+1, 2*(N_k+1), 3*(N_k+1)]   
tick_locs = [i*(N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$',fontsize = 14) 
plt.title('TM band',fontsize = 14) 
plt.savefig('SiO2-TM.png') 
plt.show() 

# Print the photonic band of band even to a file
with open('SiO2-even.txt','w') as file:

    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in zeven_freqs[n])
        file.write('\n') 

# Plot the z-even bands 
fig, ax = plt.subplots() 
ax.plot(number, zeven_freqs)       
#plt.vlines(0,0,0.2,linestyle = 'dashed',color='black')
plt.vlines(N_k+1,0,1.0,linestyle = 'dashed',color='black') 
plt.vlines(2*(N_k+1),0,1.0,linestyle = 'dashed',color='black') 
#plt.vlines(3*(N_k+1),0,0.2,linestyle = 'dashed',color='black') 
plt.xlim(0,3*(N_k+1)) 
plt.ylim(0,0.8)   
#tick_locs = [ 0, N_k+1, 2*(N_k+1), 3*(N_k+1)]   
tick_locs = [i*(N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$',fontsize = 14) 
plt.title('z-even band',fontsize = 14) 
plt.savefig('SiO2-zeven.png') 
plt.show() 

# Print the photonic band of band odd to a file
with open('SiO2-zodd.txt','w') as file:

    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in zodd_freqs[n])
        file.write('\n') 

# Plot the z-odd bands 
fig, ax = plt.subplots() 
ax.plot(number, zodd_freqs)       
#plt.vlines(0,0,0.2,linestyle = 'dashed',color='black')
plt.vlines(N_k+1,0,1.0,linestyle = 'dashed',color='black') 
plt.vlines(2*(N_k+1),0,1.0,linestyle = 'dashed',color='black') 
#plt.vlines(3*(N_k+1),0,0.2,linestyle = 'dashed',color='black') 
plt.xlim(0,3*(N_k+1)) 
plt.ylim(0,0.8)    
#tick_locs = [ 0, N_k+1, 2*(N_k+1), 3*(N_k+1)]   
tick_locs = [i*(N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$',fontsize = 14) 
plt.title('z-odd band',fontsize = 14) 
plt.savefig('SiO2-zodd.png') 
plt.show() 