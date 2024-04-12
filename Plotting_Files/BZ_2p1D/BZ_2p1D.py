import numpy as np 
import scipy 
import matplotlib.pyplot as plt 

ax = plt.figure().add_subplot(projection='3d')

# Line MA
x1,y1,z1 = np.linspace(0,0.1,100),np.linspace(0,0,100),np.linspace(0,0,100)
plt.plot(x1,y1,z1,linestyle='dashed')

# Line AB
x2,y2,z2 = np.linspace(0.1,0.1,100),np.linspace(0.0,0.1,100),np.linspace(0,0,100) 
plt.plot(x2,y2,z2,linestyle='dashed')

# Line BC
x3,y3,z3 = np.linspace(0,0.1,100),np.linspace(0.1,0.1,100),np.linspace(0,0,100)
plt.plot(x3,y3,z3,linestyle='dashed') 

# Line CM
x4,y4,z4 = np.linspace(0,0,100),np.linspace(0,0.1,100),np.linspace(0,0,100)
plt.plot(x4,y4,z4,linestyle='dashed')

# Line MB
x5,y5,z5 = np.linspace(0,0.1,100),np.linspace(0,0.1,100),np.linspace(0,0,100)
plt.plot(x5,y5,z5,linestyle='dashed')

# Line RD
x6,y6,z6 = np.linspace(0,0.1,100),np.linspace(0,0,100),np.linspace(0.5,0.5,100)
plt.plot(x6,y6,z6,linestyle='dashed')

# Line DE
x7,y7,z7 = np.linspace(0.1,0.1,100),np.linspace(0.0,0.1,100),np.linspace(0.5,0.5,100) 
plt.plot(x7,y7,z7,linestyle='dashed')

# Line EF
x8,y8,z8 = np.linspace(0,0.1,100),np.linspace(0.1,0.1,100),np.linspace(0.5,0.5,100)
plt.plot(x8,y8,z8,linestyle='dashed') 

# Line FR
x9,y9,z9 = np.linspace(0,0,100),np.linspace(0,0.1,100),np.linspace(0.5,0.5,100)
plt.plot(x9,y9,z9,linestyle='dashed')

# Line RE
x10,y10,z10 = np.linspace(0,0.1,100),np.linspace(0,0.1,100),np.linspace(0.5,0.5,100)
plt.plot(x10,y10,z10,linestyle='dashed')

# Line MR
x11,y11,z11 = np.linspace(0,0,100),np.linspace(0,0,100),np.linspace(0,0.5,100)
plt.plot(x11,y11,z11,linestyle='dashed')

# Line AD
x12,y12,z12 = np.linspace(0.1,0.1,100),np.linspace(0,0,100),np.linspace(0,0.5,100)
plt.plot(x12,y12,z12,linestyle='dashed')

# Line BE
x13,y13,z13 = np.linspace(0.1,0.1,100),np.linspace(0.1,0.1,100),np.linspace(0,0.5,100)
plt.plot(x13,y13,z13,linestyle='dashed')

# Line CF
x14,y14,z14 = np.linspace(0,0,100),np.linspace(0.1,0.1,100),np.linspace(0,0.5,100)
plt.plot(x14,y14,z14,linestyle='dashed')

# Line ME
x15,y15,z15 = np.linspace(0,0.1,100),np.linspace(0,0.1,100),np.linspace(0,0.5,100)
plt.plot(x15,y15,z15,linestyle='dashed')

### Label the point
ax.text(-0.02,-0.02,0,"M",fontsize=14)
ax.text(0.1,-0.02,0.0,"A",fontsize=14)
ax.text(0.12,0.12,0.0,"B",fontsize=14)
ax.text(-0.02,0.1,0.0,"C",fontsize=14)
ax.text(-0.02,-0.02,0.5,"R",fontsize=14)
ax.text(0.1,-0.02,0.5,"D",fontsize=14)
ax.text(0.12,0.12,0.5,"E",fontsize=14)
ax.text(-0.02,0.1,0.5,"F",fontsize=14)

ax.set_xlim(-0.01,0.11)
ax.set_ylim(-0.01,0.11)
ax.set_zlim(-0.01,0.51)
ax.set_xticks([0,0.05,0.1])
ax.set_yticks([0,0.05,0.1])
ax.set_zticks([0,0.1,0.2,0.3,0.4,0.5])
ax.set_xlabel(r'$q_x$',fontsize=14)
ax.set_ylabel(r'$q_y$',fontsize=14)
ax.set_zlabel('$q$',fontsize=14)

ax.view_init(elev=20., azim=35, roll=0)
plt.savefig('BZ_2p1D.png')
plt.show()