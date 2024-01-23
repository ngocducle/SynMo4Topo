import numpy as np 
import meep as mp 
import matplotlib.pyplot as plt 

# Load the light cone 
LightCone = np.loadtxt('Air-TE.txt')
LightConeMomenta = LightCone[:,0]
LightConeEnergy = LightCone[:,1]

# Number of k-points
N_k = 19 

# Load the Energy of the photonic bands
PhotonicBands = np.loadtxt('2Dmonolayer-h_0.5-r_0.2-res_16-TE.txt')
number = PhotonicBands[:,0]
Bands = PhotonicBands[:,1:18] 

# Plot the figure
fig, ax = plt.subplots() 
plt.plot(number, Bands) 
plt.plot(LightConeMomenta, LightConeEnergy, linewidth = 2, color = 'black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.xlim(0,3 * (N_k+1)) 
plt.ylim(0, 0.8)    
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('TE-bands', fontsize = 14) 
plt.savefig('2Dmonolayer-h_0.5-r_0.2-res_16-TE.png') 
plt.show() 

