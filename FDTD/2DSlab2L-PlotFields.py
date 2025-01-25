import numpy as np 
import scipy 
import matplotlib.pyplot as plt 

##### =============================================================================
##### Plot the field profile of 2D photonic crystal slab 2L
##### Analyze the data taken from the program 2DSlab2L-RHoleP-hj-FieldProfile
##### =============================================================================

### Load the data 
x_array = np.loadtxt('x_array.txt')
y_array = np.loadtxt('y_array.txt')
z_array = np.loadtxt('z_array.txt')

Nz = len(z_array)
print('Nz = '+str(Nz))

# Load the dielectric profile 
dielectric = []
for i in range(Nz):
    data = np.loadtxt(f'dielectric-z_{i:d}.txt')
    dielectric.append(data)

# Rearrange the axes of dielectric
dielectric = np.transpose(np.array(dielectric),(1,2,0))
print('Shape of dielectric: ')
print(np.shape(dielectric))

# Load the Ey field profile 
Field_Ey = []
for i in range(Nz):
    data = np.loadtxt(f'Field_Ey-z_{i:d}.txt',dtype=complex)
    Field_Ey.append(data)

Field_Ey = np.transpose(np.array(Field_Ey),(1,2,0))
print('Shape of Field_Ey: ')
print(np.shape(Field_Ey))

### The z coordinates of the planes to take the field
dist = 0.1
h = 0.35

z1 = 0.5*(dist+h)
z2 = 0
z3 = -z1 

# The corresponding indices
def find_nearest(array, values):
    indices = np.abs(np.subtract.outer(array, values)).argmin(0)
    return indices

ind1 = find_nearest(z_array,z1)
print(ind1)
ind2 = find_nearest(z_array,z2)
print(ind2)
ind3 = find_nearest(z_array,z3)
print(ind3)

### Plot the field profile 
X,Y = np.meshgrid(x_array,y_array)

fig,ax = plt.subplots()
ax.pcolormesh(X,Y,np.real(Field_Ey[:,:,ind1].T),shading='gouraud',cmap='RdBu')
ax.set_xlabel('x',fontsize=14)
ax.set_ylabel('y',fontsize=14)
ax.set_title(f'z = {z1:.4f}',fontsize=14)
ax.set_aspect('equal')
plt.show()