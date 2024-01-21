import numpy as np 
import meep as mp 
from meep import mpb 
import matplotlib.pyplot as plt 

# Resolution
resolution = 16 # pixels / a 

# Geometrical parameters 
h = 0.25            # Thickness of one layer 
radius = 0.2             # Radius of the hole 
dist = 0.00            # Distance between two layers 
delta1 = 0.3             # Relative displacement along the x-direction 
delta2 = 0.4             # Relative displacement along the y-direction 
Lz = 2*h + dist + 5  # Size of the cell along the z-direction

# Number of bands 
num_bands = 10  

# Number of k-points to interpolate between 2 high symmetry points 
N_k = 19    

# Define the materials 
Si   = mp.Medium(index = 3.54)
SiO2 = mp.Medium(index = 1.46)
Dielectric = mp.Medium(epsilon = 12.0) 
Air  = mp.Medium(epsilon = 1.0)

# Define the lattice
geometry_lattice = mp.Lattice(size = mp.Vector3(1, 1, Lz), 
                              basis1 = mp.Vector3(1, 0),
                              basis2 = mp.Vector3(0, 1))

# Define the geometry 
geometry = [ mp.Block(center = mp.Vector3(0, 0, 0), 
                      size = mp.Vector3(mp.inf, mp.inf, mp.inf), 
                      material = SiO2),

             mp.Block(center = mp.Vector3(0, 0, 0), 
                      size = mp.Vector3(1, 1, 2*h + dist), 
                      material = Si), 

             mp.Cylinder(center = (0.5*delta1, 0.5*delta2, 0.5*(dist + h)), 
                         radius = radius,  
                         height = h, 
                         axis = mp.Vector3(0, 0, 1), 
                         material = SiO2),

             mp.Cylinder(center = (-0.5*delta1, -0.5*delta2, -0.5*(dist + h)), 
                         radius = radius,  
                         height = h, 
                         axis = mp.Vector3(0, 0, 1), 
                         material = SiO2)   
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
                    num_bands  = num_bands 
                    )

ms.run_te()
te_freqs = ms.all_freqs 

ms.run_tm() 
tm_freqs = ms.all_freqs 

# Check the geometry
# This part doesn't work for 3D unit cell 
#md  = mpb.MPBData(rectify = True, periods = 5, resolution = 32)
#eps = ms.get_epsilon() 
#converted_eps = md.convert(eps) 

#plt.imshow(converted_eps.T, interpolation = 'spline36', cmap = 'binary')
#plt.axis('off')
#plt.show() 

# Plot the photonic band
#fig, ax = plt.subplots()
#plt.plot(k_points, te_freqs[:,0])
#plt.show() 

# Define the grid for the unit cell
Nx = 300 
Ny = 300 

gridx = np.linspace(-2.5, 2.5, Nx)
gridy = np.linspace(-2.5, 2.5, Ny)
gridz = [ -2*h, \
          -(h+dist) - 1e-4, -(h+dist) + 1e-4, \
          -0.5 * (h+dist), \
          -0.5*dist - 1e-4, -0.5*dist + 1e-4, \
           0.0, \
           0.5*dist - 1e-4, 0.5*dist + 1e-4, \
           0.5 * (h+dist), \
           h + dist - 1e-4, h + dist + 1e-4, \
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
    plt.savefig('2Dbilayer-h_'+str(h)+'-r_'+str(radius)+ \
        '-delta1_'+str(delta1)+'-delta2_'+str(delta2)+ \
        '-res_'+str(resolution)+'-z_'+str(k)+'.png')
    plt.show() 

# Print the photonic band of band TE to a file
with open('2Dbilayer-h_'+str(h)+'-r_'+str(radius)+'-delta1_' \
    +str(delta1)+'-delta2_'+str(delta2)+'-res'+str(resolution) \
    +'-TE.txt','w') as file:
    
    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in te_freqs[n])
        file.write('\n') 

number = np.arange(len(k_points)) 

# Plot the TE bands 
fig, ax = plt.subplots() 
ax.plot(number, te_freqs)      
#plt.vlines(0, 0, 0.2, linestyle = 'dashed', color='black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color='black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color='black') 
#plt.vlines(3 * (N_k+1), 0, 0.2, linestyle = 'dashed', color='black') 
plt.xlim(0,3 * (N_k+1))  
plt.ylim(0, 0.5)   
#tick_locs = [ 0, N_k+1, 2 * (N_k+1), 3 * (N_k+1)]   
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$',fontsize = 14) 
plt.title('TE band',fontsize = 14) 
plt.savefig('2Dbilayer-h_'+str(h)+'-r_'+str(radius)+ \
    '-delta1_'+str(delta1)+'-delta2_'+str(delta2)+ \
    '-res_'+str(resolution)+'-TE.png') 
plt.show() 

# Print the photonic band of band TM to a file
with open('2Dbilayer-h_'+str(h)+'-r_'+str(radius)+ \
    '-delta1_'+str(delta1)+'-delta2_'+str(delta2)+ \
    '-res_'+str(resolution)+'-TM.txt','w') as file:
    
    for n in range(len(k_points)):
        file.write('%.8f   ' % n)    
        file.writelines('%.8f    ' % w for w in tm_freqs[n])
        file.write('\n') 

# Plot the TM bands 
fig, ax = plt.subplots() 
ax.plot(number, tm_freqs)       
#plt.vlines(0, 0, 0.2, linestyle = 'dashed', color='black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color='black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color='black') 
#plt.vlines(3 * (N_k+1), 0, 0.2, linestyle = 'dashed', color='black') 
plt.xlim(0,3 * (N_k+1))  
plt.ylim(0, 0.5)    
#tick_locs = [ 0, N_k + 1, 2 * (N_k+1), 3 * (N_k+1)]   
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$',fontsize = 14) 
plt.title('TM band',fontsize = 14) 
plt.savefig('2Dbilayer-h_'+str(h)+'-r_'+str(radius)+ \
    '-delta1_'+str(delta1)+'-delta2_'+str(delta2)+ \
    '-res'+str(resolution)+'-TM.png')  
plt.show() 