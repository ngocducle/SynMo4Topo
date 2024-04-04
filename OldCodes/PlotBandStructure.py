import numpy as np 
import meep as mp 
import matplotlib.pyplot as plt 

# Name save 
namesave = '2Dmonolayer-h_0.5-r_0.1'   

# Load the light cone 
LightCone = np.loadtxt('SiO2-TE.txt')
LightConeMomenta = LightCone[:,0]
LightConeEnergy = LightCone[:,1]

# Number of k-points
N_k = 19 

# Load the Energy of the photonic bands
PhotonicBands = np.loadtxt(namesave+'-TE.txt')
number = PhotonicBands[:,0]
Bands = PhotonicBands[:,1:10] 

# Plot the figure
fig, ax = plt.subplots() 
plt.plot(number, Bands, linewidth = 3, color = 'green')  
plt.plot(LightConeMomenta, LightConeEnergy, linewidth = 3, color = 'black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.xlim(0,3 * (N_k+1)) 
plt.ylim(0, 0.5)     
ax.fill_between(LightConeMomenta, LightConeEnergy, 0.5, facecolor = 'green')  
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('TE band', fontsize = 14) 
plt.savefig(namesave+'-TE.png')
plt.show() 

# Load the Energy of the photonic bands
PhotonicBands = np.loadtxt(namesave+'-TM.txt')
number = PhotonicBands[:,0]
Bands = PhotonicBands[:,1:10] 

# Plot the figure
fig, ax = plt.subplots() 
plt.plot(number, Bands, linewidth = 3, color = 'green')  
plt.plot(LightConeMomenta, LightConeEnergy, linewidth = 3, color = 'black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.xlim(0,3 * (N_k+1)) 
plt.ylim(0, 0.5)     
ax.fill_between(LightConeMomenta, LightConeEnergy, 0.5, facecolor = 'green')  
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('TM band', fontsize = 14) 
plt.savefig(namesave+'-TM.png')
plt.show() 

# Load the Energy of the photonic bands
PhotonicBands = np.loadtxt(namesave+'-zeven.txt')
number = PhotonicBands[:,0]
Bands = PhotonicBands[:,1:10] 

# Plot the figure
fig, ax = plt.subplots() 
plt.plot(number, Bands, linewidth = 3, color = 'green')  
plt.plot(LightConeMomenta, LightConeEnergy, linewidth = 3, color = 'black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.xlim(0,3 * (N_k+1)) 
plt.ylim(0, 0.5)     
ax.fill_between(LightConeMomenta, LightConeEnergy, 0.5, facecolor = 'green')  
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('z-even band', fontsize = 14) 
plt.savefig(namesave+'-zeven.png')
plt.show() 

# Load the Energy of the photonic bands
PhotonicBands = np.loadtxt(namesave+'-zodd.txt')
number = PhotonicBands[:,0]
Bands = PhotonicBands[:,1:10] 

# Plot the figure
fig, ax = plt.subplots() 
plt.plot(number, Bands, linewidth = 3, color = 'green')  
plt.plot(LightConeMomenta, LightConeEnergy, linewidth = 3, color = 'black')
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.xlim(0,3 * (N_k+1)) 
plt.ylim(0, 0.5)     
ax.fill_between(LightConeMomenta, LightConeEnergy, 0.5, facecolor = 'green')  
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('z-odd band', fontsize = 14) 
plt.savefig(namesave+'-zodd.png')
plt.show() 
