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

### The anisotropic parameter (-1 <= alpha <= 1)
alpha = 0.2 

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
xAB = np.linspace(min(xA,xB),max(xA,xB),Nk)
yAB = np.linspace(min(yA,yB),max(yA,yB),Nk)

### Line BC
xBC = np.linspace(min(xB,xC),max(xB,xC),Nk)
yBC = np.linspace(min(yB,yC),max(yB,yC),Nk)

### Line CD
xCD = np.linspace(min(xC,xD),max(xC,xD),Nk)
yCD = np.linspace(min(yC,yD),max(yC,yD),Nk)

### Line DA
xDA = np.linspace(min(xD,xA),max(xD,xA),Nk)
yDA = np.linspace(min(yD,yA),max(yD,yA),Nk) 

### Domain LC
xLC = np.linspace(-0.5,xC,Nk)

### Domain CD 
yCD_Lower = yC + (yB-yC)/(xB-xC)*(xCD - xC)
yCD_Upper = yC + (yD-yC)/(xD-xC)*(xCD - xC)

### Domain DB
xDB = np.linspace(xD,xB,Nk)
yDB_Lower = yB + (yB-yC)/(xB-xC)*(xDB - xB)
yDB_Upper = yD + (yA-yD)/(xA-xD)*(xDB - xD)

### Domain BA 
yAB_Lower = yB + (yA-yB)/(xA-xB)*(xAB - xB)
yAB_Upper = yA + (yA-yD)/(xA-xD)*(xAB - xA) 

### Domain AR
xAR = np.linspace(xA,0.5,Nk)

### The diagonal 
xDiag = np.linspace(-0.5,0.5,Nk)
yDiag = np.linspace(-0.5,0.5,Nk)

### Plot the unit cell
fig,ax = plt.subplots(figsize=(10,8))

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

ax.fill_between(xLC,-0.5,0.5,color='yellow')
ax.fill_between(xCD,-0.5,yCD_Lower,color='yellow')
ax.fill_between(xCD,yCD_Upper,0.5,color='yellow')
ax.fill_between(xDB,-0.5,yDB_Lower,color='yellow')
ax.fill_between(xDB,yDB_Upper,0.5,color='yellow')
ax.fill_between(xAB,-0.5,yAB_Lower,color='yellow')
ax.fill_between(xAB,yAB_Upper,0.5,color='yellow')
ax.fill_between(xAR,-0.5,0.5,color='yellow')

"""plt.arrow(-0.025,0.025,d1,d1,head_width=0.03,length_includes_head=True)
plt.arrow(-0.025,0.025,-d1,-d1,head_width=0.03,length_includes_head=True)

plt.arrow(0.025,0.025,d2,-d2,head_width=0.03,length_includes_head=True)
plt.arrow(0.025,0.025,-d2,d2,head_width=0.03,length_includes_head=True)"""

plt.xlim(-0.6,0.6)
plt.ylim(-0.6,0.6)
plt.xticks([-0.5,-0.25,0.0,0.25,0.5],fontsize=18)
plt.yticks([-0.5,-0.25,0.0,0.25,0.5],fontsize=18)
plt.xlabel('x/a',fontsize=18)
plt.ylabel('y/a',fontsize=18)
plt.title(r'$\alpha > 0$',fontsize=18)
ax.set_aspect('equal')
plt.savefig('Square_Unit_cell_Rhombus_hole_diag.png')

### Show the figure
plt.show()
