import numpy as np
import scipy 
import matplotlib.pyplot as plt

##### Plot the field profile of 1Dgrating2L

### Load the data 
data = np.loadtxt('dft_field_along_x_Re.txt')
x_array = data[:,0]
Efield = data[:,1]

### Plot the figure
fig,ax = plt.subplots()
ax.plot(x_array,Efield)
plt.show()