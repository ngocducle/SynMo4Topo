import numpy as np
from numpy.linalg import eigvalsh,det
import scipy.linalg as sla
import matplotlib.pyplot as plt
from cmath import exp
from time import time

### ===================================================================================
### Start counting time
start_time = time()

##### =================================================================================
##### FUNCTION: Bulk Hamiltonian
##### The Hamiltonian at (k,q1,q2)
def Hamiltonian(k,q1,q2,domega,v,dv,U,delta,V1,V2):
    U1 = U
    U3 = U
    U2 = U + delta

    return np.array(
        [
            [v*k,U1*exp(1j*q1),V1,0,0,0],
            [U1*exp(-1j*q1),-v*k,0,V1,0,0],
            [V1,0,domega+(v+dv)*k,U2,V2,0],
            [0,V1,U2,domega-(v+dv)*k,0,V2],
            [0,0,V2,0,v*k,U3*exp(-1j*q2)],
            [0,0,0,V2,U3*exp(1j*q2),-v*k]
        ]
    )

##### FUNCTION: Matrix A of the generalized eigenvalue problem
def Mat(E,q1,q2,domega,U,delta,V1,V2):
    U1 = U
    U3 = U
    U2 = U + delta

    return np.array(
        [
            [-E,  U1*exp(1j*q1), V1, 0,          0, 0],
            [U1*exp(-1j*q1), -E, 0, V1,          0, 0],
            [V1, 0,              domega-E, U2,   V2, 0],
            [0, V1,              U2, domega-E,   0, V2],
            [0, 0,               V2, 0,          -E, U3*exp(-1j*q2)],
            [0, 0,               0, V2,          U3*exp(1j*q2), -E]
        ]
    )

##### FUNCTION: Find momentum eigenvalues
##### H(k)*x = E*x  <=>  A*x = k*B*x  with  A = H(0)-E  and  B = -diag(velocities)
##### The eigenvectors are normalized column by column
def FindMomentumEig(E,q1,q2,domega,v,dv,U,delta,V1,V2):
    A = Mat(E,q1,q2,domega,U,delta,V1,V2)
    B = np.diag([-v,v,-(v+dv),v+dv,-v,v])

    eigvals,eigvecs = sla.eig(A,B)
    eigvecs = eigvecs/np.linalg.norm(eigvecs,axis=0)

    return eigvals,eigvecs


##### ===================================================================================
##### MAIN program
### Parameters
v = 0.3283
U = 0.0207
omega0 = 0.2413+0.0053
domega = -omega0*0.55*(0.764-0.8)
dv = -v*0.37*(0.764-0.8)

Delta = 0.2*U
V1 = U
V2 = U

gap = 3 # 3rd gap

### The array of genuine momenta
Nk = 101
k_array = np.linspace(-0.25,0.25,Nk)

### The array of synthetic momenta
Nq = 101
q_array = np.linspace(-0.03*2*np.pi,0.03*2*np.pi,Nq)

### Criterion to check if the determinant is zero
epsilon = 1e-3

### Increment in energy from the band edge
epsilonE = 1e-4

### Number of E-value to scan
NE = 501

### Arrays of obstructed bands, corresponding to transmission spectrum
ObsMax_array = np.zeros(Nq)
ObsMin_array = np.zeros(Nq)

### Arrays of smaller and greater of band 3 between L and R
BulkMax_array = np.zeros(Nq)
BulkMin_array = np.zeros(Nq)

### Lists of edge states: synthetic momentum and energy
Qedge = []
EdgeStates = []

##### We scan the synthetic momentum
for iq in range(Nq):
    ### The value of the synthetic momentum
    q1 = q_array[iq]
    q2 = -q1

    ### Calculate the bulk band structure
    E_bulk_L = np.zeros((Nk,6))
    E_bulk_R = np.zeros((Nk,6))

    for ik in range(Nk):
        ### Genuine momentum
        k = k_array[ik]

        ### The bulk band structure (H is Hermitian: use eigvalsh)
        E_bulk_L[ik,:] = eigvalsh(Hamiltonian(k,q1,q2,domega,v,dv,U,Delta,V1,V2))
        E_bulk_R[ik,:] = eigvalsh(Hamiltonian(k,q1,q2,domega,v,dv,U,-Delta,V1,V2))

    ### The obstructed bands, corresponding to transmission spectrum
    ObsMax_array[iq] = max(np.amin(E_bulk_L[:,gap]),np.amin(E_bulk_R[:,gap]))
    ObsMin_array[iq] = min(np.amax(E_bulk_L[:,gap-1]),np.amax(E_bulk_R[:,gap-1]))

    ### The closest bands
    BulkMax_array[iq] = min(np.amin(E_bulk_L[:,gap]),np.amin(E_bulk_R[:,gap]))
    BulkMin_array[iq] = max(np.amax(E_bulk_L[:,gap-1]),np.amax(E_bulk_R[:,gap-1]))

    ### The arrays of energy to scan, inside the common gap
    E_array = np.linspace(BulkMin_array[iq]+epsilonE,BulkMax_array[iq]-epsilonE,NE)

    ### Array of determinants
    S = np.zeros(NE)

    ### We scan the array E
    for ie in range(NE):
        ### Value of energy
        E = E_array[ie]

        ### Left-hand side (x<0): the modes decaying at x -> -infinity, Im(k)<0
        eigvalsL, eigvecsL = FindMomentumEig(E,q1,q2,domega,v,dv,U,Delta,V1,V2)

        ### Right-hand side (x>0): the modes decaying at x -> +infinity, Im(k)>0
        eigvalsR, eigvecsR = FindMomentumEig(E,q1,q2,domega,v,dv,U,-Delta,V1,V2)

        ### Matching at x = 0: W has 3 left and 3 right evanescent modes
        W = np.concatenate((eigvecsL[:,np.imag(eigvalsL)<0],-eigvecsR[:,np.imag(eigvalsR)>0]),axis=1)

        ### Inside the gap there are 3 modes on each side; otherwise skip E
        if (np.shape(W) != (6,6)):
            S[ie] = np.inf
            continue

        S[ie] = np.abs(det(W))

    ### Edge states: the local minima of S below the criterion epsilon
    for ie in range(1,NE-1):
        if ((S[ie]<epsilon) and (S[ie]<=S[ie-1]) and (S[ie]<=S[ie+1])):
            Qedge.append(q1)
            EdgeStates.append(E_array[ie])

print('Number of edge states found: '+str(len(EdgeStates)))
print('Elapsed time: '+str(time()-start_time)+' s')

##### Plot the figure
fig,ax = plt.subplots(figsize=(9,12))
ax.plot(q_array,ObsMax_array+omega0,color='green')
ax.plot(q_array,ObsMin_array+omega0,color='green')
ax.plot(q_array,BulkMax_array+omega0,color='darkgrey',linewidth=4)
ax.plot(q_array,BulkMin_array+omega0,color='darkgrey',linewidth=4)
ax.plot(Qedge,np.array(EdgeStates)+omega0,'o',markerfacecolor='red',markeredgecolor='red')
ax.set_xlabel(r'$q_1$',fontsize=32)
ax.set_ylabel(r'$\omega (2\pi c /\lambda)$',fontsize=32)
plt.show()
