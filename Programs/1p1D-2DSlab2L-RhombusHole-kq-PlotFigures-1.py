import numpy as np 
import scipy 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

# The upper layer 
h1 = 0.35 # Thickness of the upper layer 
b1 = 0.38  # The edge length of the undeformed square hole 
e1 = -0.1 # The anisotropy between the two diagonals of the upper layer 

# The lower layer 
h2 = 0.35 # Thickness of the lower layer 
b2 = 0.38  # The edge length of the undeformed square hole  
e2 = -0.1  # The anisotropy between the two diagonals of the lower layer

# The interlayer parameters
dist = 0.10   # Distance between the two layers 
Lz = 10.0     # The height of the unit cell along the z-direction 

### The number of k-points (genuine momenta) to interpolate 
Nk = 19 

### The number of q-points (synthetic momenta)
Nq = 51

### The array of synthetic momenta 
q_array = []

### The array of genuine momenta 
k_array = []

### Bands 
Band1 = []
Band2 = []
Band3 = []
Band4 = []
Band5 = []
Band6 = []
Band7 = []
Band8 = []

q_array_1D = np.linspace(0.0,0.50,Nq)

##### We scan over the range of q points
for iq in range(Nq):
    q = q_array_1D[iq]
    ### The title and the name of the file 
    #namesave = '1p1D-2DSlab2L-RHole-h1_0.35-b1_0.4-e1_0.05-h2_0.35-b2_0.4-e2_0.05-d_0.4_kq-Band-No_'+str(iq)+'.txt'
    
    #namesave = '1p1D-2DSlab2L-RHole-h1_'+str(h1)+'-b1_'+str(b1)+'-e1_'+str(e1) \
    #    + '-h2_'+str(h2)+'-b2_'+str(b2)+'-e2_'+str(e2) \
    #    + '-d_'+str(dist)+'_kq-Band-No_'+str(iq)+'.txt'

    namesave = '1p1D-2DSlab2L-RHoleP-h1_'+str(h1)+'-b1_'+str(b1)+'-e1_'+str(e1) \
        + '-h2_'+str(h2)+'-b2_'+str(b2)+'-e2_'+str(e2) \
        + '-d_'+str(dist)+'_kq'
    namesavetxt = namesave+'-Band'
    namesavetxt = namesavetxt + "-q__{0:.4f}".format(q) 
    namesavetxt = namesavetxt + '.txt'

    ### Load the data file 
    data = np.loadtxt(namesavetxt)

    ### Add the data 
    q_array.append(data[:,0])
    k_array.append(data[:,1])
    Band1.append(data[:,2])
    Band2.append(data[:,3])
    Band3.append(data[:,4])
    Band4.append(data[:,5])
    Band5.append(data[:,6])
    Band6.append(data[:,7])
    Band7.append(data[:,8])
    Band8.append(data[:,9])

q_array = 2*np.array(q_array)
k_array = np.array(k_array)-0.5
Band1 = np.array(Band1)
Band2 = np.array(Band2)
Band3 = np.array(Band3)
Band4 = np.array(Band4)
Band5 = np.array(Band5)
Band6 = np.array(Band6)
Band7 = np.array(Band7)
Band8 = np.array(Band8)

##### Reshape the arrays of data
q_array = q_array.reshape((Nq,Nk+2),order='C')
k_array = k_array.reshape((Nq,Nk+2),order='C')

#print(q_array)
#print(k_array)

Band1 = Band1.reshape((Nq,Nk+2),order='C') 
Band2 = Band2.reshape((Nq,Nk+2),order='C')
Band3 = Band3.reshape((Nq,Nk+2),order='C')
Band4 = Band4.reshape((Nq,Nk+2),order='C')
Band5 = Band5.reshape((Nq,Nk+2),order='C')
Band6 = Band6.reshape((Nq,Nk+2),order='C') 
Band7 = Band7.reshape((Nq,Nk+2),order='C')
Band8 = Band8.reshape((Nq,Nk+2),order='C')

##### Plot the band structure 
fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                      figsize=(12,10))
ax.plot_surface(q_array,k_array,Band1,antialiased=True)
ax.plot_surface(q_array,k_array,Band2,antialiased=True)
ax.plot_surface(q_array,k_array,Band3,antialiased=True)
ax.plot_surface(q_array,k_array,Band4,antialiased=True)
#ax.plot_surface(q_array,k_array,Band5)
#ax.plot_surface(q_array,k_array,Band6)
#ax.plot_surface(q_array,k_array,Band7)
#ax.plot_surface(q_array,k_array,Band8)
ax.set_xlabel('q',fontsize=14)
ax.set_ylabel(r'$k_x=k_y$',fontsize=14)
ax.set_title('e = '+str(e1),fontsize=14)
ax.view_init(elev=10,azim=60,roll=0)
plt.savefig(namesave+'-Band.png')
plt.show()
