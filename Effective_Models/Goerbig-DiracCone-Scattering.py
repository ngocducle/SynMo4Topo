import numpy as np 
import scipy
import scipy.linalg as sla  
import matplotlib.pyplot as plt 

### The parameters 
m = 2
c = 0.1
Delta = -0.1

### Array of momenta 
Nx = 501 
qx_array = np.linspace(-1,1,Nx)
Ny = 501 
qy_array = np.linspace(-1,1,Ny)

### The arrays of energy 
E_array = np.zeros((Nx,Ny,2))

### The arrays of states
States = np.zeros((Nx,Ny,2,2),dtype=complex)

### We scan over qx and qy 
for ix in range(Nx):
    for iy in range(Ny):
        qx = qx_array[ix]
        qy = qy_array[iy]

        H = np.zeros((2,2),dtype=complex)

        H[0,1] = -(qx+1j*qy)*(qx+1j*qy)/(2*m) + c*(qx-1j*qy) + Delta 
        H[1,0] = -(qx-1j*qy)*(qx-1j*qy)/(2*m) + c*(qx+1j*qy) + np.conjugate(Delta)

        #print('Hamiltonian = ')
        #print(H)

        # Diagonalize the Hamiltonian
        evals,evecs = sla.eigh(H)

        E_array[ix,iy,:] = evals 
        States[ix,iy,:,:] = evecs 

### Plot the band structure 
X,Y = np.meshgrid(qx_array,qy_array)

fig,ax = plt.subplots(subplot_kw={'projection':'3d'})
ax.plot_surface(X,Y,E_array[:,:,0].T)
ax.plot_surface(X,Y,E_array[:,:,1].T)
ax.set_xlabel('qx',fontsize=14)
ax.set_ylabel('qy',fontsize=14)
plt.show()

##### ========================================================================
#####       Calculate the Berry winding
##### ========================================================================
### The center of the circular contour 
center_x = 0.7
center_y = 0.0 

### The radius of the contour
contour_xaxis = 0.2
contour_yaxis = 0.2

### The array of angles 
Ntheta = 200 
dtheta = 2*np.pi/Ntheta 
theta_array = np.arange(0,2*np.pi,dtheta)
#print('theta_array = ')
#print(theta_array)

### The array of q and p of the contour 
contour_x_array = center_x + contour_xaxis * np.cos(theta_array)
contour_y_array = center_y + contour_yaxis * np.sin(theta_array)

### The list of eigenstate
U_list = []

### Scan around the contour 
for j in range(Ntheta):
    # The genuine momentum 
    qx = contour_x_array[j]

    # The synthetic momentum 
    qy = contour_y_array[j]

    H = np.zeros((2,2),dtype=complex)

    H[0,1] = -(qx+1j*qy)*(qx+1j*qy)/(2*m) + c*(qx-1j*qy) + Delta 
    H[1,0] = -(qx-1j*qy)*(qx-1j*qy)/(2*m) + c*(qx+1j*qy) + np.conjugate(Delta)

    ### Diagonalize the Hamiltonian
    E,states = sla.eigh(H)

    ### The eigenstate of the lowest band
    U_list.append(states[:,0])

### Calculate the Berry winding
Berry_winding = 0

for j in range(Ntheta-1):
    Berry_winding = Berry_winding + np.angle( np.vdot(U_list[j],U_list[j+1]) )

Berry_winding = Berry_winding + np.angle( np.vdot(U_list[Ntheta-1],U_list[0]) )
Berry_winding = Berry_winding / np.pi 

print('Berry winding = ')
print(Berry_winding)

### Plot the contour 
fig,ax = plt.subplots()
ax.contour(X,Y,E_array[:,:,0].T)
for j in range(Ntheta):
    ax.plot(contour_x_array[j],contour_y_array[j],'o',markerfacecolor='red',markeredgecolor='red',markersize=1)
ax.set_xlabel('qx',fontsize=14)
ax.set_ylabel('qy',fontsize=14)
ax.set_title(f'Delta = {Delta:.4f}, w = {Berry_winding:.4f}',fontsize=14)
plt.savefig('Berry_winding_qx_0.5.png')
plt.show()