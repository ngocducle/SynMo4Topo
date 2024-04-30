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
dist = 0.8881
print('Interlayer distance dist = '+str(dist))

Ve = V0e*np.exp(-dist/d0e) 

# The array of intrinsic momenta 
k = 0.0 

# The array of synthetic momenta 
Ndelta = 1001
delta_array = np.linspace(-0.5,0.5,Ndelta)

# Initiate the arrays for frequencies
Eeven1 = np.empty(Ndelta,dtype=complex)
Eeven2 = np.empty(Ndelta,dtype=complex)
Eeven3 = np.empty(Ndelta,dtype=complex)
Eeven4 = np.empty(Ndelta,dtype=complex)
Eeven5 = np.empty(Ndelta,dtype=complex)
Eeven6 = np.empty(Ndelta,dtype=complex)
Eeven7 = np.empty(Ndelta,dtype=complex)
Eeven8 = np.empty(Ndelta,dtype=complex)

### Calculate the energy bands of bands arsing from 1L even modes 
for j in range(Ndelta):
    # Extract the synthetic momentum delta 
    delta = delta_array[j]

    evalues,evectors = EffModel_2DSlab2L_M(k,k,delta,delta,omegaMe,ve,Ue,We,Ve)

    Eeven1[j] = evalues[0]
    Eeven2[j] = evalues[1]
    Eeven3[j] = evalues[2]
    Eeven4[j] = evalues[3]
    Eeven5[j] = evalues[4]
    Eeven6[j] = evalues[5]
    Eeven7[j] = evalues[6]
    Eeven8[j] = evalues[7]

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
Eodd1 = np.empty(Ndelta,dtype=complex)
Eodd2 = np.empty(Ndelta,dtype=complex)
Eodd3 = np.empty(Ndelta,dtype=complex)
Eodd4 = np.empty(Ndelta,dtype=complex)
Eodd5 = np.empty(Ndelta,dtype=complex)
Eodd6 = np.empty(Ndelta,dtype=complex)
Eodd7 = np.empty(Ndelta,dtype=complex)
Eodd8 = np.empty(Ndelta,dtype=complex)

### Calculate the energy bands of bands arsing from 1L even modes 
for j in range(Ndelta):
    # Extract the intrinsic momentum k and the synthetic momentum delta 
    delta = delta_array[j]

    evalues,evectors = EffModel_2DSlab2L_M(k,k,delta,delta,omegaMo,vo,Uo,Wo,Vo)

    Eodd1[j] = evalues[0]
    Eodd2[j] = evalues[1]
    Eodd3[j] = evalues[2]
    Eodd4[j] = evalues[3]
    Eodd5[j] = evalues[4]
    Eodd6[j] = evalues[5]
    Eodd7[j] = evalues[6]
    Eodd8[j] = evalues[7]

fig,ax = plt.subplots(figsize=(12,10)) 

#plt.plot(delta_array,Eeven1,color='blue')
plt.plot(delta_array,Eeven2,color='purple')
plt.plot(delta_array,Eeven3,color='red')
#plt.plot(delta_array,Eeven4,color='green')
#plt.plot(delta_array,Eeven5,color='green')
#plt.plot(delta_array,Eeven6,color='red')
#plt.plot(delta_array,Eeven7,color='purple')
#plt.plot(delta_array,Eeven8,color='blue')

#plt.plot(delta_array,Eodd1,color='blue')
#plt.plot(delta_array,Eodd2,color='purple')
#plt.plot(delta_array,Eodd3,color='red')
#plt.plot(delta_array,Eodd4,color='green')
#plt.plot(delta_array,Eodd5,color='green')
#plt.plot(delta_array,Eodd6,color='red')
#plt.plot(delta_array,Eodd7,color='purple')
#plt.plot(delta_array,Eodd8,color='blue')

ax.set_xticks([-0.50,-0.25,0,0.25,0.50])
ax.set_xlabel(r'$\delta$',fontsize=14)
ax.set_title('k = 0.5, dist = '+str(dist),fontsize=14)
plt.savefig('k_050-dist_'+str(dist)+'.png')
plt.show()