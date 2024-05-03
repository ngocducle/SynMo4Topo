import numpy as np 
import matplotlib.pyplot as plt

### Values of omega at M 
omega_M = 0.48431971

### Effective group velocity 
v = omega_M * np.sqrt(2.0)

### The arrays of kx and ky 
Nk = 500
kx = np.linspace(0.4,0.6,Nk)
ky = np.linspace(0.0,1.0,Nk)

### The array of the energy of the light cone 
LightCone = np.zeros((Nk,Nk))

### The positions of the Gamma points (in unit: 2pi/a)
### (0,0), (1,0), (0,1), (1,1)

for i in range(Nk):
    for j in range(Nk):
        kxval = kx[i]
        kyval = ky[j]

        if (kxval > 0.5):
            kxval = kxval - 1.0 
        
        if (kyval > 0.5):
            kyval = kyval - 1.0 
                
        LightCone[i,j] = v*np.sqrt(kxval**2+kyval**2)


X, Y = np.meshgrid(kx,ky)

### Plot the light cone 
ax = plt.figure(figsize=(12,10)).add_subplot(projection='3d')
lightcone = ax.plot_surface(X,Y,LightCone)
ax.set_xlabel(r'$k_x a / (2 \pi)$',fontsize=14)
ax.set_ylabel(r'$k_y a / (2 \pi)$',fontsize=14)
plt.show()