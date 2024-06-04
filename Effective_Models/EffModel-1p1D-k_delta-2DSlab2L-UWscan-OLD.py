import numpy as np 
import scipy 
import scipy.linalg as sla 
import cmath 
import matplotlib.pyplot as plt 

import sys 
sys.path.insert(0,'../src/')
from EffectiveModels import EffModel_2DSlab2L_M 

### =============================================================================== ###
#                                                                                     #
#       Calculate the 2D band structure for 2D photonic crystal slab bilayer          #
#       in the vicinity of the M point (0.5,0.5) with kx = ky = k,                    #
#       delta_x = delta_y = delta                                                     #
#                                                                                     #
### =============================================================================== ### 

### The array of U/W 
N_UW = 21
#UW = np.linspace(-2,2,N_UW)
UW = [0.0,0.5,1.0,1.5,2.0,5.0]

### The array of V/W 
N_VW = 21 
#VW = np.linspace(-2,2,N_VW)
VW = [0.0,0.5,1.0,1.5,2.0,5.0]

### The array of intrinsic momenta 
Nk = 51
k_array = np.linspace(-0.02,0.02,Nk)

### The array of synthetic momenta 
Ndelta = 51
delta_array = np.linspace(-0.5,0.5,Ndelta)

##### The light cone in the vicinity of the M point 
### The frequency of the light-cone at Gamma point
kLightCone_Gamma = 0.0 
LightCone_Gamma = 0.0 

### The frequency of the ligth cone at the X point 
kxLightCone_X = 0.5
kyLightCone_X = 0.0 
LightCone_X = 0.34246575

### Interpolate the array of light-cone
kxLightCone_array = k_array + 0.5
kyLightCone_array = k_array + 0.5
LightCone_array = np.zeros(Nk)

for i in np.arange(Nk):
    kxval = kxLightCone_array[i]
    kyval = kyLightCone_array[i]

    if (kxval > 0.5):
        kxval = kxval - 1.0 

    if (kyval > 0.5):
        kyval = kyval - 1.0 

    LightCone_array[i] = LightCone_X/kxLightCone_X * np.sqrt(kxval**2+kyval**2)

# Transfer to 2D
LightCone = np.zeros((Nk,Ndelta))

for j in np.arange(Ndelta):
    LightCone[:,j] = LightCone_array 

##### The parameters to the effective model 
### Frequency at the M point 
omegaM = 0.00
print('omegaM = '+str(omegaM))

### Group velocity 
v = 0.35 
print('Group velocity v = '+str(v))

### The value of W 
W = 0.01 

##### We scan over the arrays of UW and VW 
for i_UW in range(N_UW):
    for i_VW in range(N_VW):
        U = UW[i_UW]*W 
        V = VW[i_VW]*W  

        ### Initiate the arrays for frequencies 
        Eeven1 = np.empty((Nk,Ndelta),dtype=complex)
        Eeven2 = np.empty((Nk,Ndelta),dtype=complex)
        Eeven3 = np.empty((Nk,Ndelta),dtype=complex)
        Eeven4 = np.empty((Nk,Ndelta),dtype=complex)
        Eeven5 = np.empty((Nk,Ndelta),dtype=complex)
        Eeven6 = np.empty((Nk,Ndelta),dtype=complex)
        Eeven7 = np.empty((Nk,Ndelta),dtype=complex)
        Eeven8 = np.empty((Nk,Ndelta),dtype=complex)

        for i in range(Nk):
            for j in range(Ndelta):
                # Extract the intrinsic momentum k and the synthetum momentum delta
                k = k_array[i]
                delta = delta_array[j]

                evalues,evectors = EffModel_2DSlab2L_M(k,0.0,delta,0.0,omegaM,v,U,W,V)

                Eeven1[i,j] = evalues[0]
                Eeven2[i,j] = evalues[1]
                Eeven3[i,j] = evalues[2]
                Eeven4[i,j] = evalues[3]
                Eeven5[i,j] = evalues[4]
                Eeven6[i,j] = evalues[5]
                Eeven7[i,j] = evalues[6]
                Eeven8[i,j] = evalues[7]

                ### Plot the dispersion surface 
                X,Y = np.meshgrid(k_array,delta_array)

                X = X.transpose()+0.5
                Y = Y.transpose()

        ax = plt.figure(figsize=(12,10)).add_subplot(projection='3d')
        #light_surface = ax.plot_surface(X,Y,LightCone,cmap='gray')
        surfe1 = ax.plot_surface(X,Y,Eeven1,cmap='winter')
        surfe2 = ax.plot_surface(X,Y,Eeven2,cmap='spring')
        surfe3 = ax.plot_surface(X,Y,Eeven3,cmap='autumn')
        surfe4 = ax.plot_surface(X,Y,Eeven4,cmap='summer')
        #surfe5 = ax.plot_surface(X,Y,Eeven5,cmap='summer')
        #surfe6 = ax.plot_surface(X,Y,Eeven6,cmap='autumn')
        #surfe7 = ax.plot_surface(X,Y,Eeven7,cmap='spring')
        #surfe8 = ax.plot_surface(X,Y,Eeven8,cmap='winter')

        #ax.set_xticks([0.55,0.525,0.50,0.475,0.45])
        #ax.set_yticks([0.0,0.25,0.5,0.75,1.0])
        ax.set_xlabel(r'$k_x a / (2 \pi)$',fontsize=14)
        ax.set_ylabel(r'$\delta_x$',fontsize=14)
        ax.set_title('U/W = '+str(UW[i_UW])+', V/W = '+str(VW[i_VW]),fontsize=14)
        ax.view_init(elev=5,azim=75,roll=0)  # View from delta side
        plt.savefig('UW_'+str(UW[i_UW])+'_VW_'+str(VW[i_VW])+'.png')
        #plt.show()