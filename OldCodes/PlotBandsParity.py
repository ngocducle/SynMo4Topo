import numpy as np 
import scipy 
import matplotlib.pyplot as plt 
from matplotlib import cm, colors 

### Load the data files
Bands = np.loadtxt('2Dbilayer-freqs.txt')
Parity = np.loadtxt('2Dbilayer-parity.txt')

Dist = Bands[:,0]   # Array of the position of the center ofthe monolayer

### Plot the figure of all the bands
cmap = 'rainbow'

fig, ax = plt.subplots(figsize=(12,12))
plt.scatter(Dist,Bands[:,1],s=5,c=Parity[:,1],cmap=cmap)
plt.scatter(Dist,Bands[:,2],s=5,c=Parity[:,2],cmap=cmap)
plt.scatter(Dist,Bands[:,3],s=5,c=Parity[:,3],cmap=cmap)
plt.scatter(Dist,Bands[:,4],s=5,c=Parity[:,4],cmap=cmap)
plt.scatter(Dist,Bands[:,5],s=5,c=Parity[:,5],cmap=cmap)
plt.scatter(Dist,Bands[:,6],s=5,c=Parity[:,6],cmap=cmap)
plt.scatter(Dist,Bands[:,7],s=5,c=Parity[:,7],cmap=cmap)
plt.scatter(Dist,Bands[:,8],s=5,c=Parity[:,8],cmap=cmap)

Parity = Parity[:,1:8]
vmin = -abs(Parity).max() 
vmax = abs(Parity).max()
norm = colors.Normalize(vmin=vmin, vmax=vmax) 
cbar=fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                  orientation='vertical',  
                  shrink=1.0, 
                  ax=ax)  
cbar.set_ticks(ticks=[vmin,0.0,vmax],
               labels=['Odd-like','0.0','Even-like'],fontsize=16)
#cbar.set_label('Parity',fontsize=16)

plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('d/a',fontsize=16)
plt.ylabel(r'$\omega a/(2\pi c)$',fontsize=16)
plt.title('y-parity',fontsize=16)
plt.show()