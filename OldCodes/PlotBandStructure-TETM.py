import numpy as np 
import meep as mp 
import matplotlib.pyplot as plt 

# Name save 
namesave = '2Dmonolayer-h_0.4-r_0.4' 

# Load the light cone
LightCone = np.loadtxt('SiO2-TE.txt')
LightConeMomenta = LightCone[:,0]
LightConeEnergy = LightCone[:,1]

# Number of k-points 
N_k = 19 

# Load the Energy of the TE bands 
TEBandsdata = np.loadtxt(namesave+'-TE.txt')
TEnumber = TEBandsdata[:,0]
TEBands = TEBandsdata[:,1:10]

# Load the Energy of the TM bands
TMBandsdata = np.loadtxt(namesave+'-TM.txt')
TMnumber = TMBandsdata[:,0]
TMBands = TMBandsdata[:,1:10]

# Plot the figure 
fig, ax = plt.subplots()
plt.plot(TEnumber, TEBands, linewidth=1, color='green')
plt.plot(TMnumber, TMBands, linewidth=1, linestyle='dashed',color='red')
plt.plot(LightConeMomenta, LightConeEnergy, linewidth = 3, color = 'black') 
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.xlim(0,3 * (N_k+1)) 
plt.ylim(0, 0.6)     
#ax.fill_between(LightConeMomenta, LightConeEnergy, 0.5, facecolor = 'green')  
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('h = 0.4, r = 0.4; TE: green, TM: red', fontsize=14) 
plt.savefig(namesave+'-TETM.png') 
plt.show() 


# Load the Energy of the TE bands 
AllBandsdata = np.loadtxt(namesave+'-All.txt')
Allnumber = AllBandsdata[:,0]
AllBands = AllBandsdata[:,1:10]

# Plot the figure 
fig, ax = plt.subplots()
plt.plot(Allnumber, AllBands, linewidth=1, color='green')
plt.plot(LightConeMomenta, LightConeEnergy, linewidth = 3, color = 'black') 
plt.vlines(N_k+1, 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.vlines(2 * (N_k+1), 0, 1.0, linestyle = 'dashed', color = 'black') 
plt.xlim(0,3 * (N_k+1)) 
plt.ylim(0, 0.6)     
#ax.fill_between(LightConeMomenta, LightConeEnergy, 0.5, facecolor = 'green')  
tick_locs = [i * (N_k+1) for i in range(4)] 
tick_labs = [r'$\Gamma$', 'X', 'M', r'$\Gamma$'] 
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs, size = 16) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('h = 0.4, r = 0.4; All bands', fontsize=14) 
plt.savefig(namesave+'-All.png') 
plt.show() 