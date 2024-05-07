import numpy as np 
import scipy 
import scipy.linalg as sla
import cmath 
import matplotlib.pyplot as plt 
from matplotlib.colors import LightSource 

import sys 
sys.path.insert(0,'../src/')
from EffectiveModels import EffModel_2DSlab2L_M 

### ================================================================================= ###
#                                                                                       #
#          The effective model for 2D photonic crystal slab bilayer                     #
#          (1+1)D with:                                                                 #
#               kx = ky = k                                                             #
#               delta_x = delta_y = delta                                               #
#                                                                                       #
### ================================================================================= ###

### The frequency at the M point
omegaM = 0.0 
print('omegaM = '+str(omegaM))

### Group velocity 
v = 0.4
print('Group velocity = '+str(v))

### Interlayer coupling strength 
W = 0.010 
print('Interlayer coupling strength W = '+str(W))

### Counter-propagating coupling strength 
UW = -2.5
U = UW * W 
print('Counter-propagating coupling strength U = '+str(U))

### Interlayer coupling strength 
VW = 0.0
V = VW * W 
print('Interlayer coupling strength V = '+str(V))

### The array of intrinsic momenta
Nk = 201 
k_array = np.linspace(0.0,0.05,Nk)

### The array of synthetic momenta 
Ndelta = 201 
delta_array = np.linspace(0.0,1.0,Ndelta)

### The frequency of the light-cone at the X point
kxLightCone_X = 0.5
kyLightCone_X = 0.0 
LightCone_X = 0.34246575

### Interpolate the array of light-cone
kxLightCone_array = k_array + 0.5
kyLightCone_array = np.zeros(Nk)
LightCone_array = np.zeros(Nk)

for i in np.arange(Nk):
    kxval = kxLightCone_array[i]
    kyval = 0.5

    if (kxval > 0.5):
        kxval = kxval - 1.0 

    if (kyval > 0.5):
        kyval = kyval - 1.0 

    LightCone_array[i] = LightCone_X/kxLightCone_X * np.sqrt(kxval**2 + kyval**2)

# Transfer to 2D 
LightCone = np.zeros((Nk,Ndelta))

for j in np.arange(Ndelta):
    LightCone[:,j] = LightCone_array 

### Initiate the arrays for frequencies
E1 = np.empty((Nk,Ndelta),dtype=complex)
E2 = np.empty((Nk,Ndelta),dtype=complex)
E3 = np.empty((Nk,Ndelta),dtype=complex)
E4 = np.empty((Nk,Ndelta),dtype=complex)
E5 = np.empty((Nk,Ndelta),dtype=complex)
E6 = np.empty((Nk,Ndelta),dtype=complex)
E7 = np.empty((Nk,Ndelta),dtype=complex)
E8 = np.empty((Nk,Ndelta),dtype=complex)

### Calculate the energy bands of bands arising from 1L modes
for i in range(Nk):
    for j in range(Ndelta):
        # Extract the genuine momentum k and the synthetic momentum delta
        k = k_array[i]
        delta = delta_array[j]

        evalues, evectors = EffModel_2DSlab2L_M(k,k,delta,delta,omegaM,v,U,W,V)

        E1[i,j] = evalues[0]
        E2[i,j] = evalues[1]
        E3[i,j] = evalues[2]
        E4[i,j] = evalues[3]
        E5[i,j] = evalues[4]
        E6[i,j] = evalues[5]
        E7[i,j] = evalues[6]
        E8[i,j] = evalues[7]

for i in range(Nk):
    for j in range(Ndelta):
        E1[i,j] = min(E1[i,j],LightCone[i,j])
        E2[i,j] = min(E2[i,j],LightCone[i,j])
        E3[i,j] = min(E3[i,j],LightCone[i,j])
        E4[i,j] = min(E4[i,j],LightCone[i,j])
        E5[i,j] = min(E5[i,j],LightCone[i,j])
        E6[i,j] = min(E6[i,j],LightCone[i,j])
        E7[i,j] = min(E7[i,j],LightCone[i,j])
        E8[i,j] = min(E8[i,j],LightCone[i,j])

### Plot the dispersion surface
X,Y = np.meshgrid(k_array,delta_array)

X = X.transpose()+0.5
Y = Y.transpose()

print('Shape of X: ')
print(np.shape(X))

print('Shape of Y:')
print(np.shape(Y))

#lightcolor = np.array([0.0, 2.0, 0.0])
rgb1 = np.tile([2.0,0.0,0.0],(E1.shape[0],E1.shape[1],1))
rgb2 = np.tile([2.0,2.0,0.0],(E2.shape[0],E2.shape[1],1))
rgb3 = np.tile([0.0,2.5,0.0],(E3.shape[0],E3.shape[1],1))
rgb4 = np.tile([0.0,0.0,3.0],(E4.shape[0],E4.shape[1],1))
rgb5 = np.tile([0.8,0.9,0.0],(E5.shape[0],E5.shape[1],1))
rgb6 = np.tile([0.1,1.2,1.5],(E6.shape[0],E6.shape[1],1))
rgb7 = np.tile([1.0,0.1,0.1],(E7.shape[0],E7.shape[1],1))
rgb8 = np.tile([0.1,1.2,0.1],(E8.shape[0],E8.shape[1],1))

#Parameters
#    ----------
#    azdeg : float, default: 315 degrees (from the northwest)
#       The azimuth (0-360, degrees clockwise from North) of the light
#        source.
#    altdeg : float, default: 45 degrees
#        The altitude (0-90, degrees up from horizontal) of the light
#        source.
ls = LightSource(azdeg=225,altdeg=30)
illuminated_surface1 = ls.shade_rgb(rgb1,E1)
illuminated_surface2 = ls.shade_rgb(rgb2,E2)
illuminated_surface3 = ls.shade_rgb(rgb3,E3)
illuminated_surface4 = ls.shade_rgb(rgb4,E4)
illuminated_surface5 = ls.shade_rgb(rgb5,E5)
illuminated_surface6 = ls.shade_rgb(rgb6,E6)
illuminated_surface7 = ls.shade_rgb(rgb5,E7)
illuminated_surface8 = ls.shade_rgb(rgb6,E8)

ax1 = plt.figure(figsize=(12,10)).add_subplot(projection='3d')
surf11 = ax1.plot_surface(X,Y,E1,
                          antialiased='True',
                          linewidth=2,
                          facecolors=illuminated_surface1)
surf12 = ax1.plot_surface(X,Y,E2,
                          antialiased='True',
                          linewidth=2,
                          facecolors=illuminated_surface2)
surf13 = ax1.plot_surface(X,Y,E3,
                          antialiased='True',
                          linewidth=2,
                          facecolors=illuminated_surface3)
surf14 = ax1.plot_surface(X,Y,E4,
                          antialiased='True',
                          linewidth=2,
                          facecolors=illuminated_surface4)
surf15 = ax1.plot_surface(X,Y,E5,
                          antialiased='True',
                          linewidth=2,
                          facecolors=illuminated_surface5)
surf16 = ax1.plot_surface(X,Y,E6,
                          antialiased='True',
                          linewidth=2,
                          facecolors=illuminated_surface6)
surf17 = ax1.plot_surface(X,Y,E7,
                          antialiased='True',
                          linewidth=2,
                          facecolors=illuminated_surface7)
surf18 = ax1.plot_surface(X,Y,E8,
                          antialiased='True',
                          linewidth=2,
                          facecolors=illuminated_surface8)

ax1.set_xticks([0.55,0.525,0.50,0.475,0.45])
ax1.set_yticks([0.0,0.25,0.5,0.75,1.0])
ax1.set_xlabel(r'$k a / (2 \pi)$',fontsize=14)
ax1.set_ylabel(r'$\delta$',fontsize=14)
ax1.set_title('v = '+str(v)+', W = '+str(W)+', U/W = '+str(UW)+', V/W = '+str(VW),
             fontsize=14)
ax1.view_init(elev=10,azim=-150,roll=0)
plt.savefig('W_{:.2f}-UsW_{:.2f}-VsW_{:.2f}-v_{:.2f}-4bands.png'.format(W,UW,VW,v))

"""
ax2 = plt.figure(figsize=(12,10)).add_subplot(projection='3d')
surf21 = ax2.plot_surface(X,Y,E1,cmap='winter')
surf22 = ax2.plot_surface(X,Y,E2,cmap='spring')
surf23 = ax2.plot_surface(X,Y,E3,cmap='autumn')
surf24 = ax2.plot_surface(X,Y,E4,cmap='summer')
surf25 = ax2.plot_surface(X,Y,E5,cmap='summer')
surf26 = ax2.plot_surface(X,Y,E6,cmap='autumn')
surf27 = ax2.plot_surface(X,Y,E7,cmap='spring')
surf28 = ax2.plot_surface(X,Y,E8,cmap='winter')
ax2.set_xlabel(r'$k a / (2 \pi)$',fontsize=14)
ax2.set_ylabel(r'$\delta$',fontsize=14)
ax2.set_title('v = '+str(v)+', W = '+str(W)+', U/W = '+str(UW)+', V/W = '+str(VW),
             fontsize=14)
ax2.view_init(elev=5,azim=5,roll=0)
plt.savefig('W_{:.2f}-UsW_{:.2f}-VsW_{:.2f}-v_{:.2f}-8bands.png'.format(W,UW,VW,v))
"""

plt.show()