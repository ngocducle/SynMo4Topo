import numpy as np 
import scipy 
import scipy.linalg as sla
import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'../src/')
from EffectiveModels import EffModel_2DSlab2L_M

### ================================================================================ ###
#                                                                                      #
#           Fit the effective model for 2D photonic crystal slab bilayer               #
#                                                                                      #
### ================================================================================ ###

# The input parameters to the model
# Energy at point M
omegaM = 0.37773038

# The group velocity 
v = 0.35

# The counter-propagating coupling strength
U = -0.02365418

# The orthogonally-propagating coupling strength
W = 0.00721296

# Shift along the x-direction
deltax = 0.0

# Shift along the y-direction
deltay = 0.0

# Interlayer distance
d = 0.05

# Interlayer coupling strength at zero interlayer distance 
V0 = 0.044

# Characteristic length for interlayer distance
d0 = 0.35

# Interlayer coupling strength
V = V0*np.exp(-d/d0)

# The array of momenta
# Gamma' = (-Kmax,-Kmax)
# X' = (-Kmax,0)
# M  = (0,0)

Nk = 29 # number of momenta
Kmax = 0.05 # maximum value of q from M-point

qxarray1 = np.linspace(-Kmax,0,Nk+2)  # Gamma' -> X' 
qxarray2 = np.linspace(0,0,Nk+2)      # X' -> M 
qxarray3 = np.linspace(0,-Kmax,Nk+2)  # M -> Gamma 
qxarray = np.concatenate( (qxarray1,
                           qxarray2[1:Nk+2],
                           qxarray3[1:Nk+2]),
                           axis = 0 ) 

qyarray1 = np.linspace(-Kmax,-Kmax,Nk+2) # Gamma' -> X' 
qyarray2 = np.linspace(-Kmax,0,Nk+2)     # X' -> M 
qyarray3 = np.linspace(0,-Kmax,Nk+2)     # M -> Gamma' 
qyarray = np.concatenate( (qyarray1,
                           qyarray2[1:Nk+2],
                           qyarray3[1:Nk+2]),
                           axis = 0 ) 

# Diagonalize the Hamiltonian
E1 = np.empty(3*Nk+4,dtype=complex)
E2 = np.empty(3*Nk+4,dtype=complex)
E3 = np.empty(3*Nk+4,dtype=complex)
E4 = np.empty(3*Nk+4,dtype=complex)
E5 = np.empty(3*Nk+4,dtype=complex)
E6 = np.empty(3*Nk+4,dtype=complex)
E7 = np.empty(3*Nk+4,dtype=complex)
E8 = np.empty(3*Nk+4,dtype=complex)

index = np.arange(3*Nk+4)
 
# We filter the following bands:
# Band1-Band8, Band2-Band9, Band6-Band7, Band12-Band13
for i in range(3*Nk+4):
    qx = qxarray[i]
    qy = qyarray[i]

    evalues,evectors = EffModel_2DSlab2L_M(qx,qy,deltax,deltay,omegaM,v,U,W,V)

    E1[i] = evalues[0]
    E2[i] = evalues[1]
    E3[i] = evalues[2]
    E4[i] = evalues[3]
    E5[i] = evalues[4]
    E6[i] = evalues[5]
    E7[i] = evalues[6]
    E8[i] = evalues[7]

# Load the band structure calculated from MPB
PhotonicBands = np.loadtxt('2DSlab2L-CircularHole-h_0.3-r_0.4-all-M-vicinity-Band.txt')

# Number of k-points
number = PhotonicBands[:,0]

# Load the energy of the photonic bands
Band1 = PhotonicBands[:,1]
Band2 = PhotonicBands[:,2]
Band3 = PhotonicBands[:,6]
Band4 = PhotonicBands[:,7]
Band5 = PhotonicBands[:,8]
Band6 = PhotonicBands[:,9]
Band7 = PhotonicBands[:,12]
Band8 = PhotonicBands[:,13]

### ===========================================================================
### Load the light-cone
LightCone = np.loadtxt('PMMA-zeven-Mvicinity-Nk29.txt')
LightConeMomenta = LightCone[:,0]
LightConeEnergy = LightCone[:,1]

### ============================================================================
### Plot the band structure
fig,ax = plt.subplots()
plt.plot(index,E1,color='blue')
plt.plot(index,E2,color='blue')
plt.plot(index,E3,color='blue')
plt.plot(index,E4,color='blue')
plt.plot(index,E5,color='blue')
plt.plot(index,E6,color='blue')
plt.plot(index,E7,color='blue')
plt.plot(index,E8,color='blue')

plt.plot(number,Band1,'o',markersize=2,
         markerfacecolor='blue',markeredgecolor='blue')
plt.plot(number,Band2,'o',markersize=2,
         markerfacecolor='blue',markeredgecolor='blue')
plt.plot(number,Band3,'o',markersize=2,
         markerfacecolor='blue',markeredgecolor='blue')
plt.plot(number,Band4,'o',markersize=2,
         markerfacecolor='blue',markeredgecolor='blue')
plt.plot(number,Band5,'o',markersize=2,
         markerfacecolor='blue',markeredgecolor='blue')
plt.plot(number,Band6,'o',markersize=2,
         markerfacecolor='blue',markeredgecolor='blue')
plt.plot(number,Band7,'o',markersize=2,
         markerfacecolor='blue',markeredgecolor='blue')
plt.plot(number,Band8,'o',markersize=2,
         markerfacecolor='blue',markeredgecolor='blue')

plt.plot(LightConeMomenta, LightConeEnergy, linewidth = 3, color = 'black')
ax.fill_between(LightConeMomenta, LightConeEnergy, 0.5, facecolor = 'green') 

plt.vlines(Nk+1, 0, 1, linestyle = 'dashed', color = 'black') 
plt.vlines(2*(Nk+1),0,1, linestyle = 'dashed', color = 'black') 
plt.xlim(0,3*(Nk+1))
plt.ylim(0.275,0.475)
tick_locs = [0, Nk+1, 2*(Nk+1), 3*(Nk+1)]  
tick_labs = [ r'$\Gamma \prime$', r'$X\prime$', 'M', r'$\Gamma \prime$' ]
ax.set_xticks(tick_locs) 
ax.set_xticklabels(tick_labs, size = 16) 
#ax.set_xlabel(r'$k a / (2 \pi)$', fontsize = 14) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.savefig('FitResults-2DSlab2L.png')  
plt.show()