import numpy as np
import scipy 
import matplotlib.pyplot as plt

### Load the data of the simulation  
namesave = '2DSlab1L-CircularHole-h_0.3-r_0.4-zeven-BZ-Band'
data = np.loadtxt(namesave+'.txt')

shape = np.shape(data)

Nk = shape[1]-2

number = data[:,0]
freqs = data[:,1:Nk]

### Load the light cone 
LightCone = np.loadtxt('PMMA-zeven-BZ-Nk19.txt')
LightConeMomenta = LightCone[:,0]
LightConeEnergy = LightCone[:,1]

fig, ax = plt.subplots()
ax.plot(number, freqs, linewidth=1,color='green')
plt.plot(LightConeMomenta, LightConeEnergy, linewidth = 3, color = 'black')
plt.vlines(Nk+1,0.0,1.0,linestyle='dashed',color='black')
plt.vlines(2*(Nk+1),0.0,1.0,linestyle='dashed',color='black')
plt.xlim(0,3*(Nk+1))
plt.ylim(0,0.5)
ax.fill_between(LightConeMomenta, LightConeEnergy, 0.8, facecolor = 'green') 
tick_locs = [i*(Nk+1) for i in range(4)]
tick_labs = [r'$\Gamma$','X','M',r'$\Gamma$']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs,size=16)
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14)
plt.title(namesave,fontsize=14)
plt.savefig(namesave+'.png')
plt.show()