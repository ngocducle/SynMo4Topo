import numpy as np 
import scipy
import scipy.linalg as sla  
import matplotlib.pyplot as plt 
from matplotlib import colors 
from matplotlib.patches import Circle 

### The parameters 
m = 2
c = 0
Delta = 0.5
Mass = 0.0001

### Array of momenta 
Nx = 501 
qx_array = np.linspace(-1.5,1.5,Nx)
Ny = 501 
qy_array = np.linspace(-1.5,1.5,Ny)

### The arrays of energy 
E_array = np.zeros((Nx,Ny,2))

### The arrays of states
States = np.zeros((Nx,Ny,2,2),dtype=complex)

### The array of Berry curvature 
F_array = np.zeros((Nx,Ny,2))

### We scan over qx and qy 
for ix in range(Nx):
    for iy in range(Ny):
        qx = qx_array[ix]
        qy = qy_array[iy]

        H = np.zeros((2,2),dtype=complex)

        H[0,0] = Mass 
        H[0,1] = -(qx+1j*qy)*(qx+1j*qy)/(2*m) + c*(qx-1j*qy) + Delta 
        H[1,0] = -(qx-1j*qy)*(qx-1j*qy)/(2*m) + c*(qx+1j*qy) + np.conjugate(Delta)
        H[1,1] = -Mass 

        #print('Hamiltonian = ')
        #print(H)

        dHx = np.zeros((2,2),dtype=complex)

        dHx[0,1] = -(qx+1j*qy)/m + c 
        dHx[1,0] = -(qx-1j*qy)/m + c

        dHy = np.zeros((2,2),dtype=complex)

        dHy[0,1] = -1j*(qx+1j*qy)/m - 1j*c 
        dHy[1,0] = 1j*(qx-1j*qy)/m + 1j*c 

        # Diagonalize the Hamiltonian
        evals,evecs = sla.eigh(H)

        E_array[ix,iy,:] = evals 
        States[ix,iy,:,:] = evecs 

        dHxe = np.matmul((evecs.conjugate()).transpose(),np.matmul(dHx,evecs))
        dHye = np.matmul((evecs.conjugate()).transpose(),np.matmul(dHy,evecs))

        for n in range(2):
                for m in range(2):
                    if (m != n):
                        val = -2.0*np.imag(dHxe[n,m]*dHye[m,n]) / (evals[n]-evals[m])**2
                        F_array[ix,iy,n] = F_array[ix,iy,n] + val

### Plot the band structure 
X,Y = np.meshgrid(qx_array,qy_array)

ls = colors.LightSource(azdeg=225, altdeg=5, 
                        hsv_min_val=0, hsv_max_val=1, 
                        hsv_min_sat=1, hsv_max_sat=0)

#fcolors1 = ls.shade(E_array[:,:,0],cmap='binary',vert_exag=0.1,blend_mode='soft')
#fcolors2 = ls.shade(E_array[:,:,1],cmap='binary',vert_exag=0.1,blend_mode='soft')

fig,ax = plt.subplots(subplot_kw={'projection':'3d'},figsize=(12,10))

x_axis_x = np.linspace(-1.8,1.8,Nx)
x_axis_y = np.linspace(0,0,Nx)
x_axis_z = np.linspace(0,0,Nx)
ax.plot(x_axis_x,x_axis_y,x_axis_z,color='black',linewidth=2)
ax.text(1.5,0.25,-0.32,r'$k_x$',fontsize=24)

y_axis_x = np.linspace(0,0,Ny)
y_axis_y = np.linspace(-1.8,1.8,Ny)
y_axis_z = np.linspace(0,0,Ny)
ax.plot(y_axis_x,y_axis_y,y_axis_z,color='black',linewidth=2)
ax.text(0,1.5,-0.4,r'$k_y$',fontsize=24)

ax.plot_surface(X,Y,E_array[:,:,0].T)
ax.plot_surface(X,Y,E_array[:,:,1].T)
ax.set_xlabel('qx',fontsize=14)
ax.set_ylabel('qy',fontsize=14)
ax.view_init(elev=5, azim=30, roll=0)
plt.axis('off')
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
plt.savefig('Berry_winding.png')
plt.show()

##### ========================================================================
###         Plot the dispersion surfaces with Berry curvature 
##### ========================================================================
##### For band 1 and band 2 only 
F_array_3D = F_array 
maxabs = abs(F_array_3D[:,:,0]).max()
vmin, vmax = -maxabs,maxabs
cmap = 'coolwarm'
norm = colors.Normalize(vmin=vmin,vmax=vmax)
scamap = plt.cm.ScalarMappable(norm=norm,cmap=cmap)
linewidth = 0 

fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                      figsize=(12,10))

fcolors = scamap.to_rgba(F_array_3D[:,:,0].T)
ax.plot_surface(X,Y,E_array[:,:,0].T,
                linewidth=linewidth,
                antialiased='False',
                rstride=5,
                cstride=5,
                facecolors=fcolors,
                cmap=cmap)

fcolors = scamap.to_rgba(F_array_3D[:,:,1].T)
ax.plot_surface(X,Y,E_array[:,:,1].T,
                linewidth=linewidth,
                antialiased='False',
                rstride=5,
                cstride=5,
                facecolors=fcolors,
                cmap=cmap)

x_axis_x = np.linspace(-1.8,1.8,Nx)
x_axis_y = np.linspace(0,0,Nx)
x_axis_z = np.linspace(0,0,Nx)
ax.plot(x_axis_x,x_axis_y,x_axis_z,color='black',linewidth=2)
ax.text(1.5,0.25,-0.32,r'$k_x$',fontsize=24)

y_axis_x = np.linspace(0,0,Ny)
y_axis_y = np.linspace(-1.8,1.8,Ny)
y_axis_z = np.linspace(0,0,Ny)
ax.plot(y_axis_x,y_axis_y,y_axis_z,color='black',linewidth=2)
ax.text(0,1.5,-0.4,r'$k_y$',fontsize=24)

ax.set_xlabel(r'$q_x$',fontsize=14)
ax.set_ylabel(r'$q_y$',fontsize=14)
ax.set_zlabel(r'$\omega a / (2 \pi c)$',fontsize=14)
"""fig.colorbar(scamap,
             orientation='vertical',
             shrink=0.4,
             location = 'right',
             ax=ax)"""
ax.view_init(elev=5,azim=30,roll=0)
plt.axis('off')
plt.savefig('Dis-BCur.png')
plt.show()
plt.close()

#####################################################################
fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                      figsize=(12,10))

x_axis_x = np.linspace(-1.8,1.8,Nx)
x_axis_y = np.linspace(0,0,Nx)
x_axis_z = np.linspace(0,0,Nx)
ax.plot(x_axis_x,x_axis_y,x_axis_z,color='black',linewidth=2)
ax.text(1.5,0.25,-0.32,r'$k_x$',fontsize=24)

y_axis_x = np.linspace(0,0,Ny)
y_axis_y = np.linspace(-1.8,1.8,Ny)
y_axis_z = np.linspace(0,0,Ny)
ax.plot(y_axis_x,y_axis_y,y_axis_z,color='black',linewidth=2)
ax.text(0,1.5,-0.4,r'$k_y$',fontsize=24)

z_axis_x = np.linspace(0,0,Nx) 
z_axis_y = np.linspace(0,0,Nx)
z_axis_z = np.linspace(-1.8,1.8,Nx)
ax.plot(z_axis_x,z_axis_y,z_axis_z,color='black',linewidth=2)
ax.text(0.2,0.2,1.5,r'$m$',fontsize=24)

ax.view_init(elev=5,azim=30,roll=0)
plt.axis('off')
plt.show()