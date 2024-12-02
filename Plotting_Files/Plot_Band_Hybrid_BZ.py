import numpy as np 
import scipy 
import matplotlib.pyplot as plt 

##### ===================================================================================
#####                           
#####       Plot the band structure in the hybrid Brillouin zone 
#####       calculated using the program 
#####       Programs/1p1D-2DSlab2L-RhombusHoleP-kq-LinehybridBZ.py
#####
##### ===================================================================================

##### Load the data 
namesave = '1p1D-2DSlab2L-RHoleP-h1_0.35-b1_0.38-e1_-0.1-h2_0.35-b2_0.38-e2_-0.1-d_0.1_kq-Band.txt'
data = np.loadtxt(namesave)

### Number of bands 
num_bands = 12 

### Array of hybrid-momenta 
shape = np.shape(data)
length = shape[0]
Number = np.arange(length)

##### Plot the figure 
fig,ax = plt.subplots()
ax.plot(Number,data[:,2:num_bands])
ax.set_ylabel(r'$\omega a / (2\pi c)$',fontsize=14)
plt.show()