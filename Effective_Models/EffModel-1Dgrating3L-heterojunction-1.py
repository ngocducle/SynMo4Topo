import numpy as np
from numpy.linalg import eigvalsh,det  
import scipy
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
def FindMomentumEig(E,q1,q2,domega,v,dv,U,delta,V1,V2):
    A = Mat(E,q1,q2,domega,U,delta,V1,V2)
    B = np.diag([-v,v,-(v+dv),v+dv,-v,v])

    eigvals,eigvecs = sla.eig(A,B)

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
Nk = 51
k_array = np.linspace(-0.25,0.25,Nk)

### The array of synthetic momenta 
Nq = 51 
q_array = np.linspace(-0.03*2*np.pi,0.03*2*np.pi,Nq)

### Criterion to check if the determinant is zero 
epsilon = 1e-3

### Increment in energy from the band edge 
epsilonE = 1e-4

### Number of E-value to scan
NE = 101

### Arrays of obstructed bands, corresponding to transmission spectrum  
ObsMax_array = np.zeros(Nq)
ObsMin_array = np.zeros(Nq)

### Arrays of smaller and greater of band 3 between L and R 
BulkMax_array = np.zeros(Nq)
BulkMin_array = np.zeros(Nq)

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

        ### The bulk band structure 
        ### ATTENTION! sla.eig gives wrong eigenvalues, leading to constant bands
        eigvals = eigvalsh(Hamiltonian(k,q1,q2,domega,v,dv,U,Delta,V1,V2))
        E_bulk_L[ik,:] = eigvals 

        eigvals = eigvalsh(Hamiltonian(k,q1,q2,domega,v,dv,U,-Delta,V1,V2))
        E_bulk_R[ik,:] = eigvals 

    ### The obstructed bands, corresponding to transmission spectrum 
    ObsMax_array[iq] = max(np.amin(E_bulk_L[:,gap]),np.amin(E_bulk_R[:,gap]))
    ObsMin_array[iq] = min(np.amax(E_bulk_L[:,gap-1]),np.amax(E_bulk_R[:,gap-1]))

    ### The closest bands 
    BulkMax_array[iq] = min(np.amin(E_bulk_L[:,gap]),np.amin(E_bulk_R[:,gap]))
    BulkMin_array[iq] = max(np.amax(E_bulk_L[:,gap-1]),np.amax(E_bulk_R[:,gap-1]))

    ### The arrays of energy to scan
    E_array = np.linspace(BulkMin_array[iq]+epsilonE,BulkMin_array[iq]-epsilonE,NE)

    ### We scan the array E 
    for ie in range(NE):
        ### Value of energy 
        E = E_array[ie]

        ### Left-hand side
        eigvalsL, eigvecsL = FindMomentumEig(E,q1,q2,domega,v,dv,U,Delta,V1,V2)
        #print('Eigenvalues: ')
        #print(eigvalsL)

        ### Right-hand side
        eigvalsR, eigvecsR = FindMomentumEig(E,q1,q2,domega,v,dv,U,-Delta,V1,V2)
        print('Eigenvalues: ')
        print(eigvalsR)






##### Plot the figure
fig,ax = plt.subplots(figsize=(9,12))
ax.plot(q_array,ObsMax_array+omega0)
ax.plot(q_array,ObsMin_array+omega0)
ax.plot(q_array,BulkMax_array+omega0)
ax.plot(q_array,BulkMin_array+omega0)
plt.show()



