import numpy as np 
import scipy 
import matplotlib.pyplot as plt

##### Load the data 
data = np.loadtxt('transmission_flux.txt')

freq_array = data[:,0]
trans_flux = data[:,1]

##### Load the normalized flux
data0 = np.loadtxt('normalized_flux.txt')

norm_flux = data0[:,1]

##### Calculate the transmission
trans = np.divide(trans_flux,norm_flux)

##### Plot the transmitted flux
fig,ax = plt.subplots()
ax.plot(freq_array,trans)
ax.set_xlabel(r'$\omega a/(2 \pi c)$',fontsize=14)
ax.set_ylabel('Transmission (a.u.)',fontsize=14)
plt.show()