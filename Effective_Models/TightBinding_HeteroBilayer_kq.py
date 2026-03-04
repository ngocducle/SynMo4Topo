import numpy as np 
import scipy 
import cmath 
import scipy.linalg as sla 
import matplotlib.pyplot as plt
from matplotlib import cm,colors 

# The geometrical parameters 
a = 0.2 # lattice parameter (nm)
d = a/2.2 # Size of basis 
alpha = -0.1 # C4 symmetry breaking parameter 
D = 0.25 # interlayer distance (nm) 
xAtoms = [0,d,d,0] # arrays of x-coordinates of atoms 
yAtoms = [0,0,d,d] # arrays of y-coordiantes of atoms 

# Tight-binding Minh model for monolayer 
t1 = -np.exp(-d**2/a**2)       # intracell hopping integrals along edge 
t2 = -np.exp(-(a-d)**2/a**2)   # intercell hopping integrals 
t3 = -0.5*np.exp(-2*d**2/a**2) # intracell hopping integrals along diagonals 
gap = 0.03                     # gap parameter 

print('t1 = '+str(t1))
print('t2 = '+str(t2)) 
print('t3 = '+str(t3))
print('Gap parameter = '+str(gap))

delta_x = 0 # shift (nm)
delta_y = 0 # shift (nm) 

# The list of momenta 
Nk = 21 
k_array = np.linspace(0,0.2,Nk)

Nq = 51 
q_array = np.linspace(0,0.5,Nq)

##### The bilayer Hamiltonian 
def TightBinding_Bilayer_Hamiltonian(t1,t2,t3,gap,alpha,a,d,D,kx,ky,delta_x,delta_y):
    Hamiltonian = np.zeros((8,8),dtype=complex)

    xAtoms = [0,d,d,0]
    yAtoms = [0,0,d,d] 

    t1up = t1*(1+gap)
    t2up = t2*(1+gap)
    t3up = t3*(1+gap)
    t1down = t1*(1-gap)
    t2down = t2*(1-gap)
    t3down = t3*(1-gap)

    # Block (0,0) 
    Hamiltonian[0,0] = 0.0 
    Hamiltonian[0,1] = t1up*cmath.exp(1j*kx*d) + t2up*cmath.exp(1j*kx*(d-a))
    Hamiltonian[0,2] = t3up*(1-alpha)*cmath.exp(1j*(kx*d+ky*d))
    Hamiltonian[0,3] = t1up*cmath.exp(1j*ky*d) + t2up*cmath.exp(1j*ky*(d-a))

    Hamiltonian[1,0] = np.conjugate(Hamiltonian[0,1])
    Hamiltonian[1,1] = 0.0
    Hamiltonian[1,2] = t1up*cmath.exp(1j*ky*d) + t2up*cmath.exp(1j*ky*(d-a))
    Hamiltonian[1,3] = t3up*(1+alpha)*cmath.exp(1j*(-kx*d+ky*d))

    Hamiltonian[2,0] = np.conjugate(Hamiltonian[0,2])
    Hamiltonian[2,1] = np.conjugate(Hamiltonian[1,2])
    Hamiltonian[2,2] = 0.0 
    Hamiltonian[2,3] = t1up*cmath.exp(-1j*kx*d) + t2up*cmath.exp(-1j*kx*(d-a))

    Hamiltonian[3,0] = np.conjugate(Hamiltonian[0,3])
    Hamiltonian[3,1] = np.conjugate(Hamiltonian[1,3])
    Hamiltonian[3,2] = np.conjugate(Hamiltonian[2,3])
    Hamiltonian[3,3] = 0.0 

    # Block (0,1) 
    for j in range(0,4):
        for k in range(0,4):
            Hamiltonian[j,4+k] = 0 

            # The sum over R 
            for nx in range(-3,4):
                for ny in range(-3,4):
                    #disp_x = nx*a + a/2 + xAtoms[j] + delta_x - xAtoms[k] 
                    #disp_y = ny*a + a/2 + yAtoms[j] + delta_y - yAtoms[k] 

                    disp_x = xAtoms[k] - xAtoms[j] - a/2 - delta_x + nx*a 
                    disp_y = yAtoms[k] - yAtoms[j] - a/2 - delta_y + ny*a 

                    c2 = D**2 / (disp_x**2 + disp_y**2 + D**2) 

                    if np.sqrt(disp_x**2 + disp_y**2) < 4*a:
                        t = np.exp(-(disp_x**2 + disp_y**2 + D**2)/a**2) * (2.5*c2-2) 
                    else:
                        t = 0

                    Hamiltonian[j,4+k] = Hamiltonian[j,4+k] + t*cmath.exp(1j*(kx*disp_x+ky*disp_y))

    # Block (1,0)
    for j in range(0,4):
        for k in range(0,4):
            Hamiltonian[4+j,k] = np.conjugate(Hamiltonian[k,4+j]) 

    # Block (1,1)
    Hamiltonian[4,4] = 0.0
    Hamiltonian[4,5] = t1down*cmath.exp(1j*kx*d) + t2down*cmath.exp(1j*kx*(d-a))
    Hamiltonian[4,6] = t3down*(1-alpha)*cmath.exp(1j*(kx*d+ky*d))
    Hamiltonian[4,7] = t1down*cmath.exp(1j*ky*d) + t2down*cmath.exp(1j*ky*(d-a))

    Hamiltonian[5,4] = np.conjugate(Hamiltonian[4,5])
    Hamiltonian[5,5] = 0.0 
    Hamiltonian[5,6] = t1down*cmath.exp(1j*ky*d) + t2down*cmath.exp(1j*ky*(d-a))
    Hamiltonian[5,7] = t3down*(1+alpha)*cmath.exp(1j*(-kx*d+ky*d))

    Hamiltonian[6,4] = np.conjugate(Hamiltonian[4,6])
    Hamiltonian[6,5] = np.conjugate(Hamiltonian[5,6])
    Hamiltonian[6,6] = 0.0 
    Hamiltonian[6,7] = t1down*cmath.exp(-1j*kx*d) + t2down*cmath.exp(-1j*kx*(d-a))

    Hamiltonian[7,4] = np.conjugate(Hamiltonian[4,7])
    Hamiltonian[7,5] = np.conjugate(Hamiltonian[5,7]) 
    Hamiltonian[7,6] = np.conjugate(Hamiltonian[6,7])
    Hamiltonian[7,7] = 0.0 


    return Hamiltonian 

### ======================================================================########################
# The list of momenta 
Gamma_x = 0 
Gamma_y = 0 
X_x = np.pi/a 
X_y = 0 
M_x = np.pi/a 
M_y = np.pi/a 

Nk = 1001 

kx_array = np.concatenate( (np.linspace(Gamma_x,X_x,Nk), np.linspace(X_x,M_x,Nk)[1:Nk],np.linspace(M_x,Gamma_x,Nk)[1:Nk]) )
ky_array = np.concatenate( (np.linspace(Gamma_y,X_y,Nk), np.linspace(X_y,M_y,Nk)[1:Nk],np.linspace(M_y,Gamma_y,Nk)[1:Nk]) )

# Arrays of energy 
E_array = np.zeros((3*Nk-2,8)) 

# We scan the momentum arrays 
for i in range(0,3*Nk-2):
    # The momenta 
    kx = kx_array[i]
    ky = ky_array[i]

    # The Hamiltonian 
    Hamiltonian = TightBinding_Bilayer_Hamiltonian(t1,t2,t3,gap,alpha,a,d,D,kx,ky,delta_x,delta_y)

    # Eigenvalues and eigenvectors 
    Evals,Evecs = sla.eigh(Hamiltonian)

    E_array[i,:] = Evals 

# Plot the band structure 
fig,ax = plt.subplots()
ax.plot(E_array[:,0:8],color='purple') 
ax.vlines(0,-8,8,color='black')
ax.vlines(Nk-1,-8,8,color='black')
ax.vlines(2*(Nk-1),-8,8,color='black')
ax.vlines(3*(Nk-1),-8,8,color='black')
tick_locs = [i*(Nk-1) for i in range(4)]
tick_labs = [r'$\Gamma$','X','M',r'$\Gamma$']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs,size=16)
ax.set_xlim(0,3*(Nk-1))
ax.set_ylim(-4,4)
ax.set_ylabel('Energy (eV)')
ax.set_title(r'$\alpha = $'+str(alpha)) 
#plt.show()

### ==============================================================================================
### Plot the 2D band structure 
Nk = 21 
k_array = np.linspace(-0.5,0.5,Nk)

Nq = 31
q_array = np.linspace(-0.5,0.5,Nq) 

E_array = np.zeros((Nk,Nq,8))
E_states = np.zeros((Nk,Nq,8,8)) 
F_array = np.zeros((Nk,Nq,8)) 

# The linked variables 
U1 = np.zeros((Nk,Nq,8),dtype=complex)
U2 = np.zeros((Nk,Nq,8),dtype=complex)

for ik in range(0,Nk):
    for iq in range(0,Nq):
        k = k_array[ik] 
        q = q_array[iq] 

        kx = k*2*np.pi/a
        ky = k*2*np.pi/a

        delta_x = q*a
        delta_y = q*a

        # The Hamiltonian 
        Hamiltonian = TightBinding_Bilayer_Hamiltonian(t1,t2,t3,gap,alpha,a,d,D,kx,ky,delta_x,delta_y)

        # Eigenvalues and eigenvectors 
        Evals,Evecs = sla.eigh(Hamiltonian)
        E_array[ik,iq,:] = Evals 
        E_states[ik,iq,:,:] = Evecs

# Calculate the Berry curvature 
for ik in range(0,Nk):
    for iq in range(0,Nq):
        # The index of the ik+1 element
        ir = (ik+1)%Nk 

        # The index of the iq+1 element 
        ju = (iq+1)%Nq 

        # Calculate U1: we sum over the bands 
        for band in range(0,8):
            U1[ik,iq,band] = np.matmul(E_states[ik,iq,:,band].conjugate(),E_states[ir,iq,:,band])
            U1[ik,iq,band] = U1[ik,iq,band]/abs(U1[ik,iq,band])

        # Calculate U2: we sum over the bands
        for band in range(0,8):
            U2[ik,iq,band] = np.matmul(E_states[ik,iq,:,band].conjugate(),E_states[ik,ju,:,band])
            U2[ik,iq,band] = U2[ik,iq,band]/abs(U2[ik,iq,band])

for ik in range(0,Nk):
    for iq in range(0,Nq):
        ir = (ik+1)%Nk 
        ju = (iq+1)%Nq 

        for band in range(0,8):
            Product = U1[ik,iq,band]*U2[ir,iq,band]*np.conjugate(U1[ik,ju,band])*np.conjugate(U2[ik,iq,band])
            F_array[ik,iq,band] = np.atan2(np.imag(Product),np.real(Product))

# Plot the band structure 
X,Y = np.meshgrid(k_array,q_array)

fig,ax = plt.subplots(subplot_kw={'projection':'3d'},figsize=(12,10))
ax.plot_surface(X,Y,E_array[:,:,0].T)
ax.plot_surface(X,Y,E_array[:,:,1].T)
"""ax.plot_surface(X,Y,E_array[:,:,2].T)
ax.plot_surface(X,Y,E_array[:,:,3].T)
ax.plot_surface(X,Y,E_array[:,:,4].T)
ax.plot_surface(X,Y,E_array[:,:,5].T)
ax.plot_surface(X,Y,E_array[:,:,6].T)
ax.plot_surface(X,Y,E_array[:,:,7].T)"""
ax.set_xlabel('ka/(2pi)')
ax.set_ylabel('q')
ax.set_title(r'$\alpha = $'+str(alpha))
#plt.show()

# Plot the Berry curvature maps 
cmap = 'seismic'
maxabs = abs(F_array[:,:,0]).max()
print(maxabs)

#rcmap = plt.colormaps.get_cmap('Reds')
#bcmap = plt.colormaps.get_cmap('Blues')
#bcmap = bcmap.reversed()

fig,ax = plt.subplots(2,1,sharex=True,figsize=(12,5))
vmin,vmax = -maxabs,maxabs 
norm = colors.Normalize(vmin=vmin,vmax=vmax)

ax[0].pcolormesh(X,Y,F_array[:,:,0].T,shading='gouraud',cmap=cmap)
ax[0].set_title('Band 1',fontsize=14)
ax[0].set_xlabel('k',fontsize=16)

fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),
             orientation='vertical',
             shrink=1.0,
             ax=ax[0])

ax[1].pcolormesh(X,Y,F_array[:,:,1].T,shading='gouraud',cmap=cmap)
ax[1].set_title('Band 2',fontsize=14)
ax[1].set_xlabel('k',fontsize=16)

ax[0].tick_params(axis='y',labelsize=14)
ax[1].tick_params(axis='y',labelsize=14)
ax[1].tick_params(axis='x',labelsize=14)

fig.colorbar(cm.ScalarMappable(norm=norm,cmap=cmap),
             orientation='vertical',
             shrink=1.0,
             ax=ax[1])

### ===============================================================================================
### Plot the band structure along a path 
Nq = 1001 
Nk = 201 
k_array = np.concatenate( (np.linspace(0,0,Nq),np.linspace(0,-0.1,Nk)[1:Nk] ) )
q_array = np.concatenate( (np.linspace(-0.5,0,Nq), np.linspace(0,0,Nk)[1:Nk] ))

Nmomentum = len(k_array)

E_array = np.zeros((Nmomentum,8))

for i in range(Nmomentum):
    k = k_array[i]
    q = q_array[i] 

    kx = (k+0.5)*2*np.pi/a  
    ky = (k+0.5)*2*np.pi/a

    delta_x = q*a 
    delta_y = q*a 

    # The Hamiltonian 
    Hamiltonian = TightBinding_Bilayer_Hamiltonian(t1,t2,t3,gap,alpha,a,d,D,kx,ky,delta_x,delta_y)

    # Eigenvalues and eigenvectors 
    Evals,Evecs = sla.eigh(Hamiltonian)
    E_array[i,:] = Evals 

# Plot the band structure 
fig,ax = plt.subplots()
ax.plot(E_array[:,0:4])
ax.vlines(0,-8,8,color='black')
ax.vlines(Nq-1,-8,8,color='black')
ax.vlines(Nq+Nk-1,-8,8,color='black')
ax.set_xlim(0,Nq+Nk-1)
ax.set_ylim(-1.5,0)
tick_locs = [0,Nq-1,Nq+Nk-1]
tick_labs = ['Q','O','K']
ax.set_xticks(tick_locs)
ax.set_xticklabels(tick_labs,size=16)
ax.set_title(r'$\alpha = $'+str(alpha)) 
plt.show()

