import numpy as np 
import meep as mp 
from meep import mpb 
import matplotlib.pyplot as plt 

# Resolution
resolution = mp.Vector3(16,16,32) # pixels / a 

# Geometrical parameters  
h = 0.5              # Thickness of one layer 
radius = 0.4         # Radius of the hole 
dist = 0.5           # Distance between two layers 
delta1 = 0.0         # Relative displacement along the x-direction 
delta2 = 0.0         # Relative displacement along the y-direction 
Lz = 10              # Size of the cell along the z-direction

# Number of bands 
num_bands = 20 

# Number of k-points to interpolate between 2 high symmetry points 
N_k = 19     

# Define the materials 
Si   = mp.Medium(index = 3.54)
SiO2 = mp.Medium(index = 1.46)
PMMA = mp.Medium(index = 1.46)
Dielectric = mp.Medium(epsilon = 12.0) 
Air  = mp.Medium(epsilon = 1.0)

Environment = PMMA 

# Define the lattice
geometry_lattice = mp.Lattice(size = mp.Vector3(1,1,Lz), 
                              basis1 = mp.Vector3(1,0),
                              basis2 = mp.Vector3(0,1))

# Define the geometry 
geometry = [ mp.Block(center = mp.Vector3(0,0,0),  
                      size = mp.Vector3(mp.inf, mp.inf, mp.inf), 
                      material = Environment),

             mp.Block(center = mp.Vector3(0,0,0.5*(h+dist)), 
                      size = mp.Vector3(1,1,h), 
                      material = Si), 

             mp.Block(center = mp.Vector3(0,0,-0.5*(h+dist)),
                      size = mp.Vector3(1,1,h),
                      material = Si), 

             mp.Cylinder(center = (0.5*delta1,0.5*delta2,0.5*(dist+h)), 
                         radius = radius,  
                         height = h, 
                         axis = mp.Vector3(0,0,1), 
                         material = Environment),

             mp.Cylinder(center = (-0.5*delta1,-0.5*delta2,-0.5*(dist+h)), 
                         radius = radius,  
                         height = h, 
                         axis = mp.Vector3(0,0,1), 
                         material = Environment)    
           ] 

# Define the k-points 
k_points = [mp.Vector3(0, 0, 0),      # Gamma 
            mp.Vector3(0.5, 0, 0),    # X 
            mp.Vector3(0.5, 0.5, 0),  # M 
            mp.Vector3(0, 0, 0),      # Gamma 
           ]

# Set of k-points to fit 
#k_points = [ mp.Vector3(0.45, 0.45, 0),  # Gamma'  
#             mp.Vector3(0.5, 0.45, 0),   # X'  
#             mp.Vector3(0.5, 0.5, 0),    # M  
#             mp.Vector3(0.45, 0.45,0)    # Gamma'   
#           ]

k_points = mp.interpolate(N_k, k_points) 
#print('# Number of points to interpolate between 2 high symmetry points:'+str(N_k))
#print('# Length of the array k_points:'+str(len(k_points)))   

# Print the k-points to file 
print(np.array(k_points)) 

with open('2Dmonolayer-h_' + str(h) + '-r_' + str(radius) \
    + '-kpoints.txt', 'w') as file:

    for i in range(3*(N_k+1)+1): 
        kvalue = np.array(k_points[i]) 
        file.write('%.8f     ' % i) 
        file.write('%.8f     ' % kvalue[0])
        file.writelines('%.8f     ' % kvalue[1])
        file.writelines('%.8f      ' % kvalue[2])
        file.write('\n') 

# Run the simulation 
ms = mpb.ModeSolver(geometry = geometry,
                    geometry_lattice = geometry_lattice, 
                    k_points = k_points,
                    resolution = resolution,
                    num_bands  = num_bands 
                    )

ms.run()
freqs = ms.all_freqs 

ms.run_zeven()
zeven_freqs = ms.all_freqs 

ms.run_zodd() 
zodd_freqs = ms.all_freqs 

# Check the geometry
# Define the grid for the unit cell
Nx = 300 
Ny = 300 

gridx = np.linspace(-2.5, 2.5, Nx)
gridy = np.linspace(-2.5, 2.5, Ny)
gridz = [ -2*h, \
          -(h+dist) - 0.1, -(h+dist) + 0.1, \
          -0.5 * (h+dist), \
          -0.5*dist - 0.1, -0.5*dist + 0.1, \
           0.0, \
           0.5*dist - 0.1, 0.5*dist +  0.1, \
           0.5 * (h+dist), \
           h + dist - 0.1, h + dist + 0.1, \
           2*h ]  

xplot = np.zeros((Nx, Ny))
for j in range(Ny):
    xplot[:,j] = gridx 

yplot = np.zeros((Nx,Ny))
for i in range(Nx):
    yplot[i,:] = gridy 

for k in gridz:
    refractiveindex = np.zeros((Nx,Ny))

    for i in range(Nx):
        for j in range(Ny):
            refractiveindex[i, j] = ms.get_epsilon_point \
                (mp.Vector3(gridx[i], gridy[j], k)) 

    fig, ax = plt.subplots()
    plt.pcolormesh(xplot, yplot, refractiveindex, \
         cmap = 'RdBu', vmin = 0, vmax = 12)
    plt.xlabel('x',fontsize = 14)
    plt.ylabel('y',fontsize = 14)
    plt.title('z = '+str(k), fontsize = 14)
    plt.colorbar()
    ax.set_aspect('equal')
    plt.savefig('2Dbilayer-h_'+str(h)+'-r_'+str(radius) \
        +'-delta1_'+str(delta1)+'-delta2_'+str(delta2) \
        +'-z_'+str(k)+'.png')
    plt.show() 

### Show the dielectric profile with respect to the z-direction at 
# the point (x,y,z)
x = 0.47
y = 0.49
Nz = 500  
zarray = np.linspace(-0.5*Lz, 0.5*Lz, Nz) # Default 
#zarray = np.linspace(-1.0,1.0, Nz)
epsilonarray = np.zeros(Nz)

for i in range(Nz):
    epsilonarray[i] = ms.get_epsilon_point(mp.Vector3(x,y,zarray[i]))

with open('2Dbilayer-h_' + str(h) + '-r_' + str(radius) \
    + '-epsilon-z.txt', 'w') as file:

    for i in range(Nz):
        file.write('%.8f     ' % x)
        file.writelines('%.8f     ' % y)
        file.writelines('%.8f      ' % zarray[i])
        file.writelines('%.8f      ' % epsilonarray[i])
        file.write('\n') 

fig = plt.figure() 
plt.plot(zarray, epsilonarray)
plt.xlim(-0.5*Lz,0.5*Lz)
plt.xlabel('z', fontsize = 14)
plt.ylabel(r'$\epsilon$', fontsize = 14) 
plt.title('x = '+str(x)+', y = '+str(y), fontsize = 14) 
plt.savefig('2Dbilayer-h_' + str(h) + '-r_' + str(radius) \
    +'-epsilon-z.png') 
plt.show() 

# Print the photonic band of all bands to a file
with open('2Dbilayer-h_'+str(h)+'-r_'+str(radius)+'-delta1_' \
    +str(delta1)+'-delta2_'+str(delta2) \
    +'-all.txt','w') as file:
    
    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in freqs[n])
        file.write('\n') 

number = np.arange(len(k_points)) 

# Plot all the bands 
fig, ax = plt.subplots() 
ax.plot(number, freqs)      
#plt.vlines(0, 0, 0.2, linestyle = 'dashed', color='black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color='black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color='black') 
#plt.vlines(3 * (N_k+1), 0, 0.2, linestyle = 'dashed', color='black') 
plt.xlim(0,3 * (N_k+1))  
plt.ylim(0, 0.5)   
#tick_locs = [ 0, N_k+1, 2 * (N_k+1), 3 * (N_k+1)]   
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
#tick_labs = [r'$\Gamma \prime$', r'X $\prime$', 'M', r'$\Gamma \prime$']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$',fontsize = 14) 
plt.title('All bands',fontsize = 14) 
plt.savefig('2Dbilayer-h_'+str(h)+'-r_'+str(radius) \
    +'-delta1_'+str(delta1)+'-delta2_'+str(delta2) \
    +'-all.png') 
plt.show() 

# Print the photonic band of zeven-bands to a file
with open('2Dbilayer-h_'+str(h)+'-r_'+str(radius)+'-delta1_' \
    +str(delta1)+'-delta2_'+str(delta2) \
    +'-zeven.txt','w') as file:
    
    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in zeven_freqs[n])
        file.write('\n') 

number = np.arange(len(k_points)) 

# Plot the z-even bands 
fig, ax = plt.subplots() 
ax.plot(number, zeven_freqs)      
#plt.vlines(0, 0, 0.2, linestyle = 'dashed', color='black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color='black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color='black') 
#plt.vlines(3 * (N_k+1), 0, 0.2, linestyle = 'dashed', color='black') 
plt.xlim(0,3 * (N_k+1))  
plt.ylim(0, 0.5)   
#tick_locs = [ 0, N_k+1, 2 * (N_k+1), 3 * (N_k+1)]   
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
#tick_labs = [r'$\Gamma \prime$', r'X $\prime$', 'M', r'$\Gamma \prime$']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$',fontsize = 14) 
plt.title('z-even band',fontsize = 14) 
plt.savefig('2Dbilayer-h_'+str(h)+'-r_'+str(radius) \
    +'-delta1_'+str(delta1)+'-delta2_'+str(delta2) \
    +'-zeven.png') 
plt.show() 

# Print the photonic band of zodd-bands to a file
with open('2Dbilayer-h_'+str(h)+'-r_'+str(radius) \
    +'-delta1_'+str(delta1)+'-delta2_'+str(delta2) \
    +'-zodd.txt','w') as file:
    
    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in zodd_freqs[n])
        file.write('\n') 

# Plot the zodd-bands 
fig, ax = plt.subplots() 
ax.plot(number, zodd_freqs)       
#plt.vlines(0, 0, 0.2, linestyle = 'dashed', color='black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color='black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color='black') 
#plt.vlines(3 * (N_k+1), 0, 0.2, linestyle = 'dashed', color='black') 
plt.xlim(0,3 * (N_k+1))  
plt.ylim(0, 0.5)    
#tick_locs = [ 0, N_k + 1, 2 * (N_k+1), 3 * (N_k+1)]   
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
#tick_labs = [r'$\Gamma \prime$', r'X $\prime$', 'M', r'$\Gamma \prime$']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$',fontsize = 14) 
plt.title('z-odd band',fontsize = 14) 
plt.savefig('2Dbilayer-h_'+str(h)+'-r_'+str(radius) \
    +'-delta1_'+str(delta1)+'-delta2_'+str(delta2) \
    +'-zodd.png')  
plt.show() 