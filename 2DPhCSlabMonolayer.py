import numpy as np 
import meep as mp 
from meep import mpb 
import matplotlib.pyplot as plt 

# Resolution
resolution = mp.Vector3(16,16,49) # pixels / a 

# Geometrical parameters 
h = 0.50       # Thickness of one layer 
radius = 0.10  # Radius of the hole 
Lz = 5         # Length of the unit cell along the z-direction 

# Number of bands 
num_bands = 10

# Number of k-points to interpolate between 2 high symmetry points 
N_k = 19        

# Define the materials 
Si = mp.Medium(index = 3.54)
SiO2 = mp.Medium(index = 1.46)
Dielectric = mp.Medium(epsilon = 12.0) 
Air = mp.Medium(epsilon = 1.0)

Environnement = Air 

# Define the lattice
geometry_lattice = mp.Lattice(size = mp.Vector3(1, 1, Lz),  
                              basis1 = mp.Vector3(1, 0),
                              basis2 = mp.Vector3(0, 1))

# Define the geometry 
geometry = [ mp.Block(center = mp.Vector3(0, 0, 0), 
                      size = mp.Vector3(mp.inf, mp.inf, mp.inf), 
                      material = Environnement),

             mp.Block(center = mp.Vector3(0, 0, 0), 
                      size = mp.Vector3(1, 1, h), 
                      material = Si),

             mp.Cylinder(center = mp.Vector3(0, 0, 0), 
                         radius = radius, 
                         height = h, 
                         axis = mp.Vector3(0, 0, 1), 
                         material = Environnement)
           ] 

# Define the k-points 
k_points = [mp.Vector3(0, 0, 0),      # Gamma 
            mp.Vector3(0.5, 0, 0),    # X 
            mp.Vector3(0.5, 0.5, 0),  # M 
            mp.Vector3(0, 0, 0),      # Gamma 
           ]

k_points = mp.interpolate(N_k, k_points) 
print('# Number of points to interpolate between 2 high symmetry points:'+str(N_k))
print('# Length of the array k_points:'+str(len(k_points)))   

# Run the simulation 
ms = mpb.ModeSolver(geometry = geometry,
                    geometry_lattice = geometry_lattice, 
                    k_points = k_points,
                    resolution = resolution,
                    num_bands = num_bands 
                    )

ms.run() 

ms.run_te()
te_freqs = ms.all_freqs 

ms.run_tm() 
tm_freqs = ms.all_freqs 

ms.run_zeven()
zeven_freqs = ms.all_freqs 

ms.run_zodd()
zodd_freqs  = ms.all_freqs 

# Define the grid for the unit cell 
Nx = 300  
Ny = 300 

gridx = np.linspace(-2.5, 2.5, Nx)
gridy = np.linspace(-2.5, 2.5, Ny)
gridz = np.linspace(-h, h, 9) 

xplot = np.zeros((Nx, Ny))
for j in range(Ny):
    xplot[:, j] = gridx 

yplot = np.zeros((Nx, Ny))
for i in range(Nx):
    yplot[i, :] = gridy 

for k in range(9): 
    refractiveindex = np.zeros((Nx, Ny))

    for i in range(Nx):
        for j in range(Ny):
            refractiveindex[i,j] = ms.get_epsilon_point \
                (mp.Vector3(gridx[i], gridy[j], gridz[k])) 

    fig, ax = plt.subplots()
    plt.pcolormesh(xplot,yplot,refractiveindex,cmap='RdBu',vmin=0,vmax=12)  
    plt.xlabel('x', fontsize = 14)
    plt.ylabel('y', fontsize = 14)
    plt.title('z = ' + str(gridz[k]), fontsize=14) 
    plt.colorbar() 
    ax.set_aspect('equal')  
    plt.savefig('2Dmonolayer-h_' + str(h) + '-r_' + str(radius) \
        + '-z_' + str(gridz[k]) + '.png')
    plt.show()  

### Show the dielectric profile with respect to the z-direction at 
# the point (x,y,z)
x = 0.27
y = 0.36
Nz = 500  
zarray = np.linspace(-0.5*Lz, 0.5*Lz, Nz) # Default 
#zarray = np.linspace(-0.3, -0.2, Nz)
epsilonarray = np.zeros(Nz)

for i in range(Nz):
    epsilonarray[i] = ms.get_epsilon_point(mp.Vector3(x,y,zarray[i]))

with open('2Dmonolayer-h_' + str(h) + '-r_' + str(radius) \
    + '-epsilon-z.txt', 'w') as file:

    for i in range(Nz):
        file.write('%.8f     ' % x)
        file.writelines('%.8f     ' % y)
        file.writelines('%.8f      ' % zarray[i])
        file.writelines('%.8f      ' % epsilonarray[i])
        file.write('\n') 

fig = plt.figure() 
plt.plot(zarray, epsilonarray)
#plt.xlim(-0.32,-0.18) 
plt.xlabel('z', fontsize = 14)
plt.ylabel(r'$\epsilon$', fontsize = 14) 
plt.title('x = '+str(x)+', y = '+str(y), fontsize = 14) 
plt.savefig('2Dmonolayer-h_' + str(h) + '-r_' + str(radius) \
    + '-epsilon-z.png') 
plt.show() 


# Print the photonic band of band TE to a file
with open('2Dmonolayer-h_' + str(h) + '-r_' + str(radius) \
    + '-TE.txt', 'w') as file:

    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in te_freqs[n])
        file.write('\n') 

number = np.arange(len(k_points)) 

# Plot the TE bands 
fig, ax = plt.subplots() 
ax.plot(number, te_freqs)      
#plt.vlines(0, 0, 0.2, linestyle = 'dashed', color = 'black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color='black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color='black') 
#plt.vlines(3 * (N_k+1), 0, 0.2, linestyle = 'dashed', color='black') 
plt.xlim(0, 3 * (N_k+1)) 
plt.ylim(0, 0.5)   
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$',fontsize = 14) 
plt.title('TE band',fontsize = 14) 
plt.savefig('2Dmonolayer-h_' + str(h) + '-r_' + str(radius) \
    + '-TE.png') 
plt.show() 

# Print the photonic band of band TM to a file
with open('2Dmonolayer-h_' + str(h) + '-r_' + str(radius) \
    +'-TM.txt', 'w') as file:

    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in tm_freqs[n])
        file.write('\n') 

# Plot the TM bands 
fig, ax = plt.subplots() 
ax.plot(number, tm_freqs)       
#plt.vlines(0, 0, 0.2, linestyle = 'dashed', color = 'black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color = 'black') 
#plt.vlines(3 * (N_k+1), 0, 0.2, linestyle = 'dashed', color = 'black') 
plt.xlim(0,3 * (N_k+1)) 
plt.ylim(0, 0.5)   
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('TM band', fontsize = 14) 
plt.savefig('2Dmonolayer-h_' + str(h) + '-r_' + str(radius) \
    + '-TM.png') 
plt.show() 

# Print the photonic band of band z-even to a file
with open('2Dmonolayer-h_' + str(h) + '-r_' + str(radius) \
    + '-zeven.txt', 'w') as file:

    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in zeven_freqs[n])
        file.write('\n') 

# Plot the z-even bands 
fig, ax = plt.subplots() 
ax.plot(number, zeven_freqs)       
#plt.vlines(0, 0, 0.2, linestyle = 'dashed', color = 'black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color = 'black') 
#plt.vlines(3 * (N_k+1), 0, 0.2, linestyle = 'dashed', color = 'black') 
plt.xlim(0,3 * (N_k+1)) 
plt.ylim(0, 0.5)     
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('z-even band', fontsize = 14) 
plt.savefig('2Dmonolayer-h_' + str(h) + '-r_' + str(radius) \
    + '-zeven.png') 
plt.show() 

# Print the photonic band of band z-odd to a file
with open('2Dmonolayer-h_' + str(h) + '-r_' + str(radius) \
    + '-zodd.txt', 'w') as file:

    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in zodd_freqs[n])
        file.write('\n') 

# Plot the z-odd bands 
fig, ax = plt.subplots() 
ax.plot(number, zodd_freqs)       
#plt.vlines(0, 0, 0.2, linestyle = 'dashed', color = 'black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color = 'black') 
#plt.vlines(3 * (N_k+1), 0, 0.2, linestyle = 'dashed', color = 'black') 
plt.xlim(0, 3 * (N_k+1)) 
plt.ylim(0, 0.5)     
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('z-odd band', fontsize = 14) 
plt.savefig('2Dmonolayer-h_' + str(h) + '-r_' + str(radius) \
    + '-zodd.png') 
plt.show() 
