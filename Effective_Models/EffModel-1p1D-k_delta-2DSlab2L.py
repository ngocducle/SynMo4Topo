import numpy as np
import scipy 
import scipy.linalg as sla
import cmath 
import matplotlib.pyplot as plt 

import sys 
sys.path.insert(0,'../src/')
from EffectiveModels import EffModel_2DSlab2L_M

### ================================================================================ ###
#                                                                                      #
#           Fit the effective model for 2D photonic crystal slab bilayer               #
#                                                                                      #
### ================================================================================ ###

### Calculate the modes arising from monolayer EVEN bands

# Frequency at the M point
omegaMe = 0.300464408
print('omegaMe = '+str(omegaMe))

# Group velocity 
ve = 0.35
print('Group velocity = '+str(ve))

# Counter-propagating coupling strength 
Ue = -0.016615673
print('Ue = '+str(Ue))

# Orthogonally-propagating coupling strength
We = 0.001744918
print('We = '+str(We))

# Interlayer coupling strength at zero interlayer distance
V0e = 0.044
print('V0e = '+str(V0e))

# Characteristic length for interlayer coupling
d0e = 0.35 
print('d0e = '+str(d0e))

# The interlayer distance
dist = 0.1
print('Interlayer distance dist = '+str(dist))

Ve = V0e*np.exp(-dist/d0e) 

# The array of intrinsic momenta
Nk = 201 
k_array = np.linspace(-0.1,0.1,Nk)
#k_array = [0.0,0.001]

# The array of synthetic momenta 
Ndelta = 201
delta_array = np.linspace(-0.5,0.5,Ndelta)
#delta_array = [0.0,0.01]

# Initiate the arrays for frequencies
Eeven1 = np.empty((Nk,Ndelta),dtype=complex)
Eeven2 = np.empty((Nk,Ndelta),dtype=complex)
Eeven3 = np.empty((Nk,Ndelta),dtype=complex)
Eeven4 = np.empty((Nk,Ndelta),dtype=complex)
Eeven5 = np.empty((Nk,Ndelta),dtype=complex)
Eeven6 = np.empty((Nk,Ndelta),dtype=complex)
Eeven7 = np.empty((Nk,Ndelta),dtype=complex)
Eeven8 = np.empty((Nk,Ndelta),dtype=complex)

### Calculate the energy bands of bands arsing from 1L even modes 
for i in range(Nk):
    for j in range(Ndelta):
        # Extract the intrinsic momentum k and the synthetic momentum delta 
        k = k_array[i]
        delta = delta_array[j]

        evalues,evectors = EffModel_2DSlab2L_M(k,k,delta,delta,omegaMe,ve,Ue,We,Ve)

        Eeven1[i,j] = evalues[0]
        Eeven2[i,j] = evalues[1]
        Eeven3[i,j] = evalues[2]
        Eeven4[i,j] = evalues[3]
        Eeven5[i,j] = evalues[4]
        Eeven6[i,j] = evalues[5]
        Eeven7[i,j] = evalues[6]
        Eeven8[i,j] = evalues[7]

#print(np.shape(Eeven1))

### Calculate the modes arising from monolayer ODD bands

# Frequency at the M point
omegaMo = 0.380151458
print('omegaMo = '+str(omegaMo))

# Group velocity 
vo = 0.41 
print('Group velocity = '+str(vo))

# Counter-propagating coupling strength 
Uo = 0.001794338
print('Uo = '+str(Uo))

# Orthogonally-propagating coupling strength
Wo = 0.006309538
print('Wo = '+str(Wo)) 

# Interlayer coupling strength amplitude
V0o = 0.049
print('V0o = '+str(V0o))

# Characteristic length for interlayer coupling
d0o = 0.42
print('d0o = '+str(d0o)) 
Vo = V0o*np.exp(-dist/d0o)

# Initiate the arrays for frequencies
Eodd1 = np.empty((Nk,Ndelta),dtype=complex)
Eodd2 = np.empty((Nk,Ndelta),dtype=complex)
Eodd3 = np.empty((Nk,Ndelta),dtype=complex)
Eodd4 = np.empty((Nk,Ndelta),dtype=complex)
Eodd5 = np.empty((Nk,Ndelta),dtype=complex)
Eodd6 = np.empty((Nk,Ndelta),dtype=complex)
Eodd7 = np.empty((Nk,Ndelta),dtype=complex)
Eodd8 = np.empty((Nk,Ndelta),dtype=complex)

"""
### Calculate the energy bands of bands arsing from 1L even modes 
for i in range(Nk):
    for j in range(Ndelta):
        # Extract the intrinsic momentum k and the synthetic momentum delta 
        k = k_array[i]
        delta = delta_array[j]

        evalues,evectors = EffModel_2DSlab2L_M(k,k,delta,delta,omegaMo,vo,Uo,Wo,Vo)

        Eodd1[i,j] = evalues[0]
        Eodd2[i,j] = evalues[1]
        Eodd3[i,j] = evalues[2]
        Eodd4[i,j] = evalues[3]
        Eodd5[i,j] = evalues[4]
        Eodd6[i,j] = evalues[5]
        Eodd7[i,j] = evalues[6]
        Eodd8[i,j] = evalues[7]
"""


### Plot the dispersion surface
X, Y = np.meshgrid(k_array,delta_array)

X = X.transpose()+0.5
Y = Y.transpose()

print(np.shape(X))
print(X)

print(np.shape(Y))
print(Y)

ax = plt.figure(figsize=(12,10)).add_subplot(projection='3d')
surfe1 = ax.plot_surface(X,Y,Eeven1,cmap='winter')
surfe2 = ax.plot_surface(X,Y,Eeven2,cmap='spring')
surfe3 = ax.plot_surface(X,Y,Eeven3,cmap='autumn')
surfe4 = ax.plot_surface(X,Y,Eeven4,cmap='summer')
#surfe5 = ax.plot_surface(X,Y,Eeven5,cmap='summer')
#surfe6 = ax.plot_surface(X,Y,Eeven6,cmap='autumn')
#surfe7 = ax.plot_surface(X,Y,Eeven7,cmap='spring')
#surfe8 = ax.plot_surface(X,Y,Eeven8,cmap='winter')

#surfo1 = ax.plot_surface(X,Y,Eodd1,cmap='Blues')
#surfo2 = ax.plot_surface(X,Y,Eodd2,cmap='Blues')
#surfo3 = ax.plot_surface(X,Y,Eodd3,cmap='Blues')
#surfo4 = ax.plot_surface(X,Y,Eodd4,cmap='Blues')

#ax.set_zlim(0.345,0.36)
ax.set_xticks([0.4,0.45,0.5,0.55,0.6])
ax.set_yticks([-0.50,-0.25,0,0.25,0.50])
ax.set_xlabel(r'$k a / (2 \pi)$',fontsize=14)
ax.set_ylabel(r'$\delta$',fontsize=14)
ax.set_title('dist = '+str(dist),fontsize=14)
ax.view_init(elev=5,azim=75,roll=0) # View from delta side
#ax.view_init(elev=5,azim=15,roll=0) # View from q side 
plt.savefig('dist_'+str(dist)+'.png')
plt.show()