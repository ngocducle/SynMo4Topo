import numpy as np 
import scipy 
import matplotlib.pyplot as plt 

### Load the data 
data = np.loadtxt('data.txt') 
#print('data = ') 
#print(data)

data_shift = data[:,0]
data_freq = data[:,1] 
data_trans = data[:,2] 

### The array of shift 
Nshift = 51
minshift = 0.0
maxshift = 0.5
shift_array = np.linspace(minshift,maxshift,Nshift) 

### The array of frequencies 
fcen = 0.2387 
df = 0.01
nfreq = 500 
freq_array = np.linspace(fcen-df,fcen+df,nfreq)

### Reshape the data 
data_shift = np.reshape(data_shift,(Nshift,nfreq)) 
print('data_shift = ')
print(data_shift) 

data_freq = np.reshape(data_freq,(Nshift,nfreq))
print('data_freq = ') 
print(data_freq) 

data_trans = np.reshape(data_trans,(Nshift,nfreq)) 
print('data_trans = ') 
print(data_trans) 

### Plot the transmission flux 
fig,ax = plt.subplots() 
plt.pcolormesh(data_shift,data_freq,data_trans,shading='gouraud',cmap='hot') 
plt.xlabel('q') 
plt.ylabel(r'$\omega a /(2 \pi c)$') 
plt.show() 


