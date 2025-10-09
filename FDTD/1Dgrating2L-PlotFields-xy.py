import numpy as np
import scipy
import matplotlib.pyplot as plt

# Import data
x_array = np.loadtxt('x_array.txt')
y_array = np.loadtxt('y_array.txt')
dielectric = np.loadtxt('dielectric.txt')
Ez = np.loadtxt('dft_field_xy_Ez.txt')

print(np.shape(dielectric))
print(np.shape(Ez))

X,Y = np.meshgrid(x_array,y_array)

print(np.shape(X))
print(np.shape(Y))

fig,ax = plt.subplots(figsize=(12,8))
ax.contour(X,Y,dielectric.T,cmap='binary')
ax.pcolormesh(X,Y,Ez.T,shading='gouraud',cmap='bwr')
ax.set_aspect('equal')
ax.set_xlabel('x')
ax.set_ylabel('y')
plt.show()
