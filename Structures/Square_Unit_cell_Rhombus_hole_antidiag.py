import numpy as np 
import scipy 
import matplotlib.pyplot as plt 

### ==================================================================
###
### Plot the unit cell for the case of diagonal symmetry-breaking 
###
### ================================================================== 

### The diagonal of the rhombus (0< d < sqrt(2)/2)
d = 0.25

### The anisotropic parameter (-1 <= alpha < 0)
alpha = -0.2 

### The length of the diagonal x=y
d1 = d*(1+alpha)

### The length of the diagonal x=-y 
d2 = d*(1-alpha) 

### The vertices of the rhombus 
xA,yA = (1+alpha)*d, (1+alpha)*d 
xB,yB = (1-alpha)*d, -(1-alpha)*d 
xC,yC = -(1+alpha)*d,-(1+alpha)*d 
xD,yD = -(1-alpha)*d,(1-alpha)*d 

#print(min(xA,yA))

### Number of points in a line 
Nk = 100 

### Line AB
xAB = np.linspace(xA,xB,Nk)
yAB = np.linspace(yA,yB,Nk)

### Line BC
xBC = np.linspace(xB,xC,Nk)
yBC = np.linspace(yB,yC,Nk)

### Line CD
xCD = np.linspace(xC,xD,Nk)
yCD = np.linspace(yC,yD,Nk)

### Line DA
xDA = np.linspace(xD,xA,Nk)
yDA = np.linspace(yD,yA,Nk) 

### Domain LD
xLD = np.linspace(-0.5,xD,Nk)

### Domain DC 
xDC = np.linspace(xD,xC,Nk)
yDC_Lower = yD + (yD-yC)/(xD-xC)*(xDC - xD)
yDC_Upper = yD + (yD-yA)/(xD-xA)*(xDC - xD)

### Domain CA
xCA = np.linspace(xC,xA,Nk)
yCA_Lower = yC + (yB-yC)/(xB-xC)*(xCA - xC)
yCA_Upper = yA + (yA-yD)/(xA-xD)*(xCA - xA)

### Domain AB 
yAB_Lower = yB + (yC-yB)/(xC-xB)*(xAB - xB)
yAB_Upper = yA + (yA-yB)/(xA-xB)*(xAB - xA) 

### Domain BR
xBR = np.linspace(xB,0.5,Nk)

### The diagonal 
xDiag = np.linspace(-0.5,0.5,Nk)
yDiag = np.linspace(-0.5,0.5,Nk)

### Plot the unit cell
fig,ax = plt.subplots()

plt.vlines(-0.5,-0.5,0.5,color='black',linewidth=2)
plt.vlines(0.5,-0.5,0.5,color='black',linewidth=2)
plt.hlines(-0.5,-0.5,0.5,color='black',linewidth=2)
plt.hlines(0.5,-0.5,0.5,color='black',linewidth=2)

plt.plot(xAB,yAB,color='black',linewidth=2)
plt.plot(xBC,yBC,color='black',linewidth=2)
plt.plot(xCD,yCD,color='black',linewidth=2)
plt.plot(xDA,yDA,color='black',linewidth=2) 

plt.plot(xDiag,yDiag,color='black',linewidth=2,linestyle='dashed')
plt.plot(xDiag,-yDiag,color='black',linewidth=2,linestyle='dashed')

ax.fill_between(xLD,-0.5,0.5,color='yellow')
ax.fill_between(xDC,-0.5,yDC_Lower,color='yellow')
ax.fill_between(xDC,yDC_Upper,0.5,color='yellow')
ax.fill_between(xCA,-0.5,yCA_Lower,color='yellow')
ax.fill_between(xCA,yCA_Upper,0.5,color='yellow')
ax.fill_between(xAB,-0.5,yAB_Lower,color='yellow')
ax.fill_between(xAB,yAB_Upper,0.5,color='yellow')
ax.fill_between(xBR,-0.5,0.5,color='yellow')

"""plt.arrow(-0.025,0.025,d1,d1,head_width=0.03,length_includes_head=True)
plt.arrow(-0.025,0.025,-d1,-d1,head_width=0.03,length_includes_head=True)

plt.arrow(0.025,0.025,d2,-d2,head_width=0.03,length_includes_head=True)
plt.arrow(0.025,0.025,-d2,d2,head_width=0.03,length_includes_head=True)"""

plt.xlim(-0.6,0.6)
plt.ylim(-0.6,0.6)
plt.xticks([-0.5,-0.25,0.0,0.25,0.5],fontsize=14)
plt.yticks([-0.5,-0.25,0.0,0.25,0.5],fontsize=14)
plt.xlabel('x/a',fontsize=14)
plt.ylabel('y/a',fontsize=14)
plt.title(r'$\alpha < 0$',fontsize=18)
ax.set_aspect('equal')
plt.savefig('Square_Unit_cell_Rhombus_hole_antidiag.png')

### Show the figure
plt.show()
