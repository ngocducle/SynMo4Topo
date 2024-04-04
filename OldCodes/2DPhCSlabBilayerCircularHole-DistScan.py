import numpy as np
import meep as mp
from meep import mpb 
import matplotlib.pyplot as plt 

### =======================================================================
### Function: calculate the band structure for 2D PhC slab bilayer
### with circular hole at the center 
### The varying parameter is the distance between the two layers 

def _2DPhCSlabBilayerCircularHole_M(h,radius,dist,delta1,delta2):
    # Resolution 
    resolution = mp.Vector3(16,16,32) 
    
    # The height of a unit cell 
    Lz = 10

    # Number of bands
    num_bands = 20 

    # Define the materials 
    Si = mp.Medium(index = 3.54)
    SiO2 = mp.Medium(index = 1.46)
    PMMA = mp.Medium(index = 1.46)
    Dielectric = mp.Medium(epsilon = 12.0)
    Air = mp.Medium(epsilon = 1.0)

    Environment = PMMA 

    # Define the lattice 
    geometry_lattice = mp.Lattice(size = mp.Vector3(1,1,Lz),
                                  basis1 = mp.Vector3(1,0),
                                  basis2 = mp.Vector3(0,1))

    # Define the geometry 
    geometry = [ mp.Block(center = mp.Vector3(0, 0, 0), 
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
                            axis = mp.Vector3(0, 0, 1), 
                            material = Environment)    
                ] 

    # Define the k-point
    k_points = [mp.Vector3(0.5,0.5,0)] # M 

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

    # Return the results 
    return freqs, zeven_freqs, zodd_freqs 

### =======================================================================
### Run the simulation 

# Thickness of one layer 
h = 0.3 

# Radius of a hole 
radius = 0.4 

# Array of distance between the two layers
Ndist = 20   
distarray = np.linspace(2.05,3.0,Ndist)

# Initialize the arrays of bands
Band1 = np.zeros(Ndist)
Band2 = np.zeros(Ndist)
Band3 = np.zeros(Ndist)
Band4 = np.zeros(Ndist)
Band5 = np.zeros(Ndist)
Band6 = np.zeros(Ndist)
Band7 = np.zeros(Ndist)
Band8 = np.zeros(Ndist) 
Band9 = np.zeros(Ndist)
Band10 = np.zeros(Ndist)
Band11 = np.zeros(Ndist)
Band12 = np.zeros(Ndist)
Band13 = np.zeros(Ndist)
Band14 = np.zeros(Ndist)
Band15 = np.zeros(Ndist)
Band16 = np.zeros(Ndist)
Band17 = np.zeros(Ndist)
Band18 = np.zeros(Ndist)
Band19 = np.zeros(Ndist)
Band20 = np.zeros(Ndist)

Bandeven1 = np.zeros(Ndist)
Bandeven2 = np.zeros(Ndist)
Bandeven3 = np.zeros(Ndist)
Bandeven4 = np.zeros(Ndist)
Bandeven5 = np.zeros(Ndist)
Bandeven6 = np.zeros(Ndist)
Bandeven7 = np.zeros(Ndist)
Bandeven8 = np.zeros(Ndist)
Bandeven9 = np.zeros(Ndist)
Bandeven10 = np.zeros(Ndist)
Bandeven11 = np.zeros(Ndist)
Bandeven12 = np.zeros(Ndist)
Bandeven13 = np.zeros(Ndist)
Bandeven14 = np.zeros(Ndist)
Bandeven15 = np.zeros(Ndist)
Bandeven16 = np.zeros(Ndist)
Bandeven17 = np.zeros(Ndist)
Bandeven18 = np.zeros(Ndist)
Bandeven19 = np.zeros(Ndist)
Bandeven20 = np.zeros(Ndist)

Bandodd1 = np.zeros(Ndist)
Bandodd2 = np.zeros(Ndist)
Bandodd3 = np.zeros(Ndist)
Bandodd4 = np.zeros(Ndist)
Bandodd5 = np.zeros(Ndist)
Bandodd6 = np.zeros(Ndist)
Bandodd7 = np.zeros(Ndist)
Bandodd8 = np.zeros(Ndist) 
Bandodd9 = np.zeros(Ndist)
Bandodd10 = np.zeros(Ndist)
Bandodd11 = np.zeros(Ndist)
Bandodd12 = np.zeros(Ndist)
Bandodd13 = np.zeros(Ndist)
Bandodd14 = np.zeros(Ndist)
Bandodd15 = np.zeros(Ndist)
Bandodd16 = np.zeros(Ndist)
Bandodd17 = np.zeros(Ndist)
Bandodd18 = np.zeros(Ndist)
Bandodd19 = np.zeros(Ndist)
Bandodd20 = np.zeros(Ndist)

# Scan over distarray and find the frequencies of the bands 
for i in np.arange(Ndist):
    dist = distarray[i]

    freqs, zeven_freqs, zodd_freqs = _2DPhCSlabBilayerCircularHole_M(\
        h,radius,dist,0,0) 

    #print(zeven_freqs) 

    Band1[i] = freqs[0,0]
    Band2[i] = freqs[0,1]
    Band3[i] = freqs[0,2]
    Band4[i] = freqs[0,3]
    Band5[i] = freqs[0,4]
    Band6[i] = freqs[0,5]
    Band7[i] = freqs[0,6]
    Band8[i] = freqs[0,7] 
    Band9[i] = freqs[0,8]
    Band10[i] = freqs[0,9]
    Band11[i] = freqs[0,10]
    Band12[i] = freqs[0,11]
    Band13[i] = freqs[0,12]
    Band14[i] = freqs[0,13]
    Band15[i] = freqs[0,14]
    Band16[i] = freqs[0,15] 
    Band17[i] = freqs[0,16]
    Band18[i] = freqs[0,17]
    Band19[i] = freqs[0,18]
    Band20[i] = freqs[0,19] 

    Bandeven1[i] = zeven_freqs[0,0]
    Bandeven2[i] = zeven_freqs[0,1]
    Bandeven3[i] = zeven_freqs[0,2]
    Bandeven4[i] = zeven_freqs[0,3]
    Bandeven5[i] = zeven_freqs[0,4]
    Bandeven6[i] = zeven_freqs[0,5]
    Bandeven7[i] = zeven_freqs[0,6]
    Bandeven8[i] = zeven_freqs[0,7] 
    Bandeven9[i] = zeven_freqs[0,8]
    Bandeven10[i] = zeven_freqs[0,9]
    Bandeven11[i] = zeven_freqs[0,10]
    Bandeven12[i] = zeven_freqs[0,11]
    Bandeven13[i] = zeven_freqs[0,12]
    Bandeven14[i] = zeven_freqs[0,13]
    Bandeven15[i] = zeven_freqs[0,14]
    Bandeven16[i] = zeven_freqs[0,15] 
    Bandeven17[i] = zeven_freqs[0,16]
    Bandeven18[i] = zeven_freqs[0,17]
    Bandeven19[i] = zeven_freqs[0,18]
    Bandeven20[i] = zeven_freqs[0,19] 

    Bandodd1[i] = zodd_freqs[0,0]
    Bandodd2[i] = zodd_freqs[0,1]
    Bandodd3[i] = zodd_freqs[0,2]
    Bandodd4[i] = zodd_freqs[0,3]
    Bandodd5[i] = zodd_freqs[0,4]
    Bandodd6[i] = zodd_freqs[0,5]
    Bandodd7[i] = zodd_freqs[0,6]
    Bandodd8[i] = zodd_freqs[0,7]
    Bandodd9[i] = zodd_freqs[0,8]
    Bandodd10[i] = zodd_freqs[0,9]
    Bandodd11[i] = zodd_freqs[0,10]
    Bandodd12[i] = zodd_freqs[0,11]
    Bandodd13[i] = zodd_freqs[0,12]
    Bandodd14[i] = zodd_freqs[0,13]
    Bandodd15[i] = zodd_freqs[0,14]
    Bandodd16[i] = zodd_freqs[0,15] 
    Bandodd17[i] = zodd_freqs[0,16]
    Bandodd18[i] = zodd_freqs[0,17]
    Bandodd19[i] = zodd_freqs[0,18]
    Bandodd20[i] = zodd_freqs[0,19] 

# Print the photonic bands into a file 
# Column 0: distarray 
# Column j (1 <= j <= 8): Bandj 
with open('2Dbilayer-alls.txt','w') as file:  

    for i in range(Ndist):
        file.write('%.8f    ' % distarray[i]) 
        file.writelines('%.8f    ' % Band1[i])
        file.writelines('%.8f    ' % Band2[i])
        file.writelines('%.8f    ' % Band3[i])
        file.writelines('%.8f    ' % Band4[i])
        file.writelines('%.8f    ' % Band5[i])
        file.writelines('%.8f    ' % Band6[i])
        file.writelines('%.8f    ' % Band7[i])
        file.writelines('%.8f    ' % Band8[i])
        file.writelines('%.8f    ' % Band9[i])
        file.writelines('%.8f    ' % Band10[i])
        file.writelines('%.8f    ' % Band11[i])
        file.writelines('%.8f    ' % Band12[i])
        file.writelines('%.8f    ' % Band13[i])
        file.writelines('%.8f    ' % Band14[i])
        file.writelines('%.8f    ' % Band15[i])
        file.writelines('%.8f    ' % Band16[i])
        file.writelines('%.8f    ' % Band17[i])
        file.writelines('%.8f    ' % Band18[i])
        file.writelines('%.8f    ' % Band19[i])
        file.writelines('%.8f    ' % Band20[i])
        file.write('\n')

with open('2Dbilayer-even.txt','w') as file:  

    for i in range(Ndist):
        file.write('%.8f    ' % distarray[i]) 
        file.writelines('%.8f    ' % Bandeven1[i])
        file.writelines('%.8f    ' % Bandeven2[i])
        file.writelines('%.8f    ' % Bandeven3[i])
        file.writelines('%.8f    ' % Bandeven4[i])
        file.writelines('%.8f    ' % Bandeven5[i])
        file.writelines('%.8f    ' % Bandeven6[i])
        file.writelines('%.8f    ' % Bandeven7[i])
        file.writelines('%.8f    ' % Bandeven8[i])
        file.writelines('%.8f    ' % Bandeven9[i])
        file.writelines('%.8f    ' % Bandeven10[i])
        file.writelines('%.8f    ' % Bandeven11[i])
        file.writelines('%.8f    ' % Bandeven12[i])
        file.writelines('%.8f    ' % Bandeven13[i])
        file.writelines('%.8f    ' % Bandeven14[i])
        file.writelines('%.8f    ' % Bandeven15[i])
        file.writelines('%.8f    ' % Bandeven16[i])
        file.writelines('%.8f    ' % Bandeven17[i])
        file.writelines('%.8f    ' % Bandeven18[i])
        file.writelines('%.8f    ' % Bandeven19[i])
        file.writelines('%.8f    ' % Bandeven20[i])
        file.write('\n')

with open('2Dbilayer-odd.txt','w') as file:  

    for i in range(Ndist):
        file.write('%.8f    ' % distarray[i]) 
        file.writelines('%.8f    ' % Bandodd1[i])
        file.writelines('%.8f    ' % Bandodd2[i])
        file.writelines('%.8f    ' % Bandodd3[i])
        file.writelines('%.8f    ' % Bandodd4[i])
        file.writelines('%.8f    ' % Bandodd5[i])
        file.writelines('%.8f    ' % Bandodd6[i])
        file.writelines('%.8f    ' % Bandodd7[i])
        file.writelines('%.8f    ' % Bandodd8[i])
        file.writelines('%.8f    ' % Bandodd9[i])
        file.writelines('%.8f    ' % Bandodd10[i])
        file.writelines('%.8f    ' % Bandodd11[i])
        file.writelines('%.8f    ' % Bandodd12[i])
        file.writelines('%.8f    ' % Bandodd13[i])
        file.writelines('%.8f    ' % Bandodd14[i])
        file.writelines('%.8f    ' % Bandodd15[i])
        file.writelines('%.8f    ' % Bandodd16[i])
        file.writelines('%.8f    ' % Bandodd17[i])
        file.writelines('%.8f    ' % Bandodd18[i])
        file.writelines('%.8f    ' % Bandodd19[i])
        file.writelines('%.8f    ' % Bandodd20[i])
        file.write('\n')

# Plot the figure for all bands 
fig, ax = plt.subplots() 
plt.plot(distarray, Band1,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Band2,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Band3,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Band4,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Band5,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Band6,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Band7,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Band8,'o',markerfacecolor='green',markeredgecolor='green')
ax.set_xlabel('d/a', fontsize = 14)
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14)
plt.title('All bands',fontsize = 14) 
plt.savefig('2Dbilayer-all.png') 
plt.show() 

# Plot the figure fpr zeven bands 
fig, ax = plt.subplots() 
plt.plot(distarray, Bandeven1,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandeven2,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandeven3,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandeven4,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandeven5,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandeven6,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandeven7,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandeven8,'o',markerfacecolor='green',markeredgecolor='green')
ax.set_xlabel('d/a', fontsize = 14)
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14)
plt.title('z-even band',fontsize = 14) 
plt.savefig('2Dbilayer-zeven.png') 
plt.show() 

# Plot the figure fpr zodd bands 
fig, ax = plt.subplots() 
plt.plot(distarray, Bandodd1,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandodd2,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandodd3,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandodd4,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandodd5,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandodd6,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandodd7,'o',markerfacecolor='green',markeredgecolor='green')
plt.plot(distarray, Bandodd8,'o',markerfacecolor='green',markeredgecolor='green')
ax.set_xlabel('d/a', fontsize = 14)
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14)
plt.title('z-odd band',fontsize = 14) 
plt.savefig('2Dbilayer-zodd.png') 
plt.show() 