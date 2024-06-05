import numpy as np 
import scipy 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

# The upper layer 
h1 = 0.35 # Thickness of the upper layer 
b1 = 0.4  # The edge length of the undeformed square hole 
d1 = b1/np.sqrt(2) # The projection of mean diagonal on the x and y axes 
e1 = 0.05  # The anisotropy between the two diagonals of the upper layer 

# The lower layer 
h2 = 0.35 # Thickness of the lower layer 
b2 = 0.4  # The edge length of the undeformed square hole 
d2 = b2/np.sqrt(2) # The projection of mean diagonal on the x and y axes 
e2 = 0.05  # The anisotropy between the two diagonals of the lower layer

# The interlayer parameters
dist = 0.40   # Distance between the two layers 
Lz = 10.0     # The height of the unit cell along the z-direction 

### The number of k-points (genuine momenta) to interpolate 
Nk = 29 

### The number of q-points (synthetic momenta)
Nq = 16

### The title and the name of the files 
namesave = '1p1D-2DSlab2L-RHole-h1_'+str(h1)+'-b1_'+str(b1)+'-e1_'+str(e1) \
    + '-h2_'+str(h2)+'-b2_'+str(b2)+'-e2_'+str(e2) \
    + '-d_'+str(dist)+'_kq' 


##### Load the data file 
data = np.loadtxt(namesave+'-Band.txt')

q_array = data[:,0].reshape((Nq,Nk),order='C')
k_array = data[:,1].reshape((Nq,Nk),order='C')

#print(q_array)
#print(k_array)

Band1 = data[:,3].reshape((Nq,Nk),order='C') 
Band2 = data[:,4].reshape((Nq,Nk),order='C')
Band3 = data[:,5].reshape((Nq,Nk),order='C')
Band4 = data[:,6].reshape((Nq,Nk),order='C')
Band5 = data[:,7].reshape((Nq,Nk),order='C')


##### Plot the band structure 
fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                      figsize=(12,10))
ax.plot_surface(q_array,k_array,Band1)
ax.plot_surface(q_array,k_array,Band2)
ax.plot_surface(q_array,k_array,Band3)
ax.plot_surface(q_array,k_array,Band4)
ax.plot_surface(q_array,k_array,Band5)
ax.set_xlabel('q',fontsize=14)
ax.set_ylabel('k',fontsize=14)
plt.show()

