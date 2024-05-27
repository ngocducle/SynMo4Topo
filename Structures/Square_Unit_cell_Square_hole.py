import numpy as np 
import scipy 
import matplotlib.pyplot as plt 

### ==================================================================
###
### Plot the unit cell for the case of diagonal symmetry-breaking 
###
### ================================================================== 

### The diagonal of the square 
d = 0.25

### Number of points in a line 
Nk = 100 

### The diagonal 
xDiag = np.linspace(-0.5,0.5,Nk)
yDiag = np.linspace(-0.5,0.5,Nk)

### Plot the unit cell
fig,ax = plt.subplots()

plt.vlines(-0.5,-0.5,0.5,color='black',linewidth=2)
plt.vlines(0.5,-0.5,0.5,color='black',linewidth=2)
plt.hlines(-0.5,-0.5,0.5,color='black',linewidth=2)
plt.hlines(0.5,-0.5,0.5,color='black',linewidth=2)

plt.vlines(-d,-d,d,color='black',linewidth=2)
plt.vlines(d,-d,d,color='black',linewidth=2)
plt.hlines(-d,-d,d,color='black',linewidth=2)
plt.hlines(d,-d,d,color='black',linewidth=2)

plt.plot(xDiag,yDiag,color='black',linewidth=2,linestyle='dashed')
plt.plot(xDiag,-yDiag,color='black',linewidth=2,linestyle='dashed')

ax.fill_between(np.linspace(-0.5,-d,Nk),-0.5,0.5,color='yellow')
ax.fill_between(np.linspace(-d,d,Nk),-0.5,-d,color='yellow')
ax.fill_between(np.linspace(-d,d,Nk),0.5,d,color='yellow')
ax.fill_between(np.linspace(d,0.5,Nk),-0.5,0.5,color='yellow')

plt.xlim(-0.6,0.6)
plt.ylim(-0.6,0.6)
plt.xticks([-0.5,-0.25,0.0,0.25,0.5],fontsize=14)
plt.yticks([-0.5,-0.25,0.0,0.25,0.5],fontsize=14)
plt.xlabel('x/a',fontsize=14)
plt.ylabel('y/a',fontsize=14)
plt.title(r'$\alpha = 0$',fontsize=18)
ax.set_aspect('equal')
plt.savefig('Square_Unit_cell_Square_hole.png')

### Show the figure
plt.show()