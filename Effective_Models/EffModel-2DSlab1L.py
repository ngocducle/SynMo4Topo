import numpy as np 
import scipy 
import scipy.linalg as sla
import matplotlib.pyplot as plt

import sys
sys.path.insert(0,'../src/')
from EffectiveModels import EffModel_2DSlab1L_M

### ================================================================================ ###
#                                                                                      #
#           Fit the effective model for 2D photonic crystal slab monolayer             #
#                                                                                      #
### ================================================================================ ###

### ====================================================================================
### EVEN MODES: 
# Values of the bands 
E1calc_e = 0.40138452
E2calc_e = 0.40156082
E3calc_e = 0.33956213
E4calc_e = 0.36841401

# Energy at point M 
omegaM_e = (E1calc_e + E2calc_e + E3calc_e + E4calc_e) / 4    
print('omegaM_e = '+str(omegaM_e))

# The group velocity 
ve = 0.35      
print('Group velocity ve = '+str(ve)) 

# The counter-propagating coupling strength 
Ue = omegaM_e - E1calc_e 
print('Ue = '+str(Ue)) 

# The orthogonally-propagating coupling strength 
We = (E4calc_e - E3calc_e) / 4
print('We = '+str(We))

# The array of momenta 
# Gamma' = (-Kmax, -Kmax)
# X' = (-Kmax, 0)
# M  = (0, 0) 

Nk = 19 # number of momenta 
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
E1_e = np.empty(3*Nk+4, dtype = complex) 
E2_e = np.empty(3*Nk+4, dtype = complex)
E3_e = np.empty(3*Nk+4, dtype = complex)
E4_e = np.empty(3*Nk+4, dtype = complex)  

index = np.arange(3*Nk+4) 

for i in range(3*Nk+4): 
    qx = qxarray[i]  
    qy = qyarray[i]

    evalues, evectors = EffModel_2DSlab1L_M(qx, qy, omegaM_e, ve, Ue, We) 

    E1_e[i] = evalues[0]
    E2_e[i] = evalues[1]
    E3_e[i] = evalues[2]
    E4_e[i] = evalues[3] 

# Load the band structure calculated from MPB
PhotonicBands_e = np.loadtxt('2DSlab1L-CircularHole-h_0.3-r_0.4-zeven-M-Band.txt') 

# Number of k-points 
number_e = PhotonicBands_e[:,0]

# Load the Energy of the photonic bands 
Band1_e = PhotonicBands_e[:,1]
Band2_e = PhotonicBands_e[:,2]
Band3_e = PhotonicBands_e[:,3]
Band4_e = PhotonicBands_e[:,4] 

### ==========================================================================
### ODD MODES
# Values of the bands 
E1calc_o = 0.42742564
E2calc_o = 0.42744911
E3calc_o = 0.41345118
E4calc_o = 0.47543789

# Energy at point M 
omegaM_o = (E1calc_o + E2calc_o + E3calc_o + E4calc_o) / 4    
print('omegaM_o = '+str(omegaM_o)) 

# The group velocity 
vo = 0.41       
print('Group velocity vo = '+str(vo)) 

# The counter-propagating coupling strength 
Uo = omegaM_o - E1calc_o 
print('Uo = '+str(Uo)) 

# The orthogonally-propagating coupling strength 
Wo = (E4calc_o - E3calc_o) / 4
print('Wo = '+str(Wo))

# The array of momenta 
# Gamma' = (-0.2, -0.2)
# X' = (-0.2, 0)
# M  = (0, 0) 

Nk = 19 # number of momenta 
Kmax = -0.05 # maximum value of q from M-point 

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
E1_o = np.empty(3*Nk+4, dtype = complex) 
E2_o = np.empty(3*Nk+4, dtype = complex)
E3_o = np.empty(3*Nk+4, dtype = complex)
E4_o = np.empty(3*Nk+4, dtype = complex)  

index = np.arange(3*Nk+4) 

for i in range(3*Nk+4): 
    qx = qxarray[i]  
    qy = qyarray[i]

    evalues, evectors = EffModel_2DSlab1L_M(qx, qy, omegaM_o, vo, Uo, Wo) 

    E1_o[i] = evalues[0]
    E2_o[i] = evalues[1]
    E3_o[i] = evalues[2]
    E4_o[i] = evalues[3] 

# Load the band structure calculated from MPB
PhotonicBands_o = np.loadtxt('2DSlab1L-CircularHole-h_0.3-r_0.4-zodd-M-Band.txt') 

# Number of k-points 
number_o = PhotonicBands_o[:,0]

# Load the Energy of the photonic bands 
Band1_o = PhotonicBands_o[:,1]
Band2_o = PhotonicBands_o[:,2]
Band3_o = PhotonicBands_o[:,3]
Band4_o = PhotonicBands_o[:,4] 

### ===========================================================================
### Load the light-cone
# Load the light cone 
LightCone = np.loadtxt('PMMA-zeven-M-Kmax005.txt')
LightConeMomenta = LightCone[:,0]
LightConeEnergy = LightCone[:,1]

### ===========================================================================
# Plot the band structure 
fig, ax = plt.subplots() 
plt.plot(index, E1_e, color = 'blue')  
plt.plot(index, E2_e, color = 'blue') 
plt.plot(index, E3_e, color = 'blue')
plt.plot(index, E4_e, color = 'blue') 
plt.plot(number_e, Band1_e, 'o', markersize = 2, 
    markerfacecolor = 'blue', markeredgecolor = 'blue')
plt.plot(number_e, Band2_e, 'o', markersize = 2, 
    markerfacecolor = 'blue', markeredgecolor = 'blue')
plt.plot(number_e, Band3_e, 'o', markersize = 2, 
    markerfacecolor = 'blue', markeredgecolor = 'blue')
plt.plot(number_e, Band4_e, 'o', markersize = 2, 
    markerfacecolor = 'blue', markeredgecolor = 'blue') 
plt.plot(index, E1_o, color = 'red')  
plt.plot(index, E2_o, color = 'red') 
plt.plot(index, E3_o, color = 'red')
plt.plot(index[34:44], E4_o[34:44], color = 'red')
plt.plot(number_o, Band1_o, 'o', markersize = 2, 
    markerfacecolor = 'red', markeredgecolor = 'red')
plt.plot(number_o, Band2_o, 'o', markersize = 2, 
    markerfacecolor = 'red', markeredgecolor = 'red')
plt.plot(number_o, Band3_o, 'o', markersize = 2, 
    markerfacecolor = 'red', markeredgecolor = 'red')
plt.plot(number_o[34:44], Band4_o[34:44], 'o', markersize = 2, 
    markerfacecolor = 'red', markeredgecolor = 'red') 
plt.plot(LightConeMomenta, LightConeEnergy, linewidth = 3, color = 'black')
ax.fill_between(LightConeMomenta, LightConeEnergy, 0.5, facecolor = 'green') 
plt.vlines(Nk+1, 0, 1, linestyle = 'dashed', color = 'black') 
plt.vlines(2*(Nk+1),0,1, linestyle = 'dashed', color = 'black') 
plt.xlim(0,3*(Nk+1))
plt.ylim(0.3,0.5)         
tick_locs = [0, Nk+1, 2*(Nk+1), 3*(Nk+1)]  
tick_labs = [ r'$\Gamma \prime$', r'$X\prime$', 'M', r'$\Gamma \prime$' ]
ax.set_xticks(tick_locs) 
ax.set_xticklabels(tick_labs, size = 16) 
#ax.set_xlabel(r'$k a / (2 \pi)$', fontsize = 14) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.savefig('FitResults-2DSlab1L.png')  
plt.show() 