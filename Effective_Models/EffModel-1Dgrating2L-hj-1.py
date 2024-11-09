import numpy as np 
import scipy 
import scipy.linalg as sla 
import matplotlib.pyplot as plt 
from cmath import exp 

##### ======================================================================
##### FUNCTION: Hamiltonian 
def Hamiltonian(k,q,omega1,v1,U1,omega2,v2,U2,V):
    H = np.array([[omega1+v1*k,U1,V*exp(-1j*np.pi*q),0],
                     [U1,omega1-v1*k,0,V*exp(1j*np.pi*q)],
                     [V*exp(1j*np.pi*q),0,omega2+v2*k,U2],
                     [0,V*exp(-1j*np.pi*q),U2,omega2-v2*k]])

    return H 

##### ========================================================================
##### FUNCTION: keig 
##### If this functions is put in list comprehension, it can either calculate 
##### kvals OR kvecs, but not both 
def keig(E,q,omega1,v1,U1,omega2,v2,U2,V):
    Vi = exp(-1j*np.pi*q)
    Vc = np.conjugate(Vi)

    ### 0th order matrix H0 
    H0 = np.array(
        [
            [omega1-E,U1,Vi,0],
            [U1,omega1-E,0,Vc],
            [Vc,0,omega2-E,U2],
            [0,Vi,U2,omega2-E]
        ]
    )

    ### 1st order matrix H1 
    H1 = np.diag([-v1,v1,-v2,v2])

    ### Find the k-eigenvalues 
    kvals,kvecs = sla.eig(H0,H1)

    ### Move the arguments of kvals to the range [0,2*pi[ 
    args = np.angle(kvals*exp(1e-9*1j))
    args[args<0] += 2*np.pi

    ### Sort kvals in increasing order of arguments
    kvals = kvals[np.argsort(args)]

    ### Sort kvecs in increasing order of arguments of kvals 
    kvecs = kvecs[:,np.argsort(args)]

    ### Normalize kvecs 
    kvecs = kvecs/np.linalg.norm(kvecs,axis=0) 
 
    return kvecs

##### ========================================================================
##### The parameters 
omega0 = 0.27815 
U0 = 0.02232 
v0 = 0.37602 

omega1 = 0.27999 
U1 = 0.02621 
v1 = 0.40590 
V1 = 0.04511 

omega2 = 0.28010 
U2 = 0.02142
v2 = 0.40454 
V2 = 0.04728

### The gap where we calculate the edge states 
gap = 0 

### Criterion to zero
epsilon = 1e-3 

### Increment in band edge 
epsilonE = 1e-4 

### Array of synthetic momenta 
Nq = 101 
q_array = np.linspace(0.45,0.55,Nq)

### Array of genuine momenta 
Nk = 251 
k_array = np.linspace(-0.12,0.12,Nk)

### Number of energies 
NE = 101

### Arrays of closest bands 
bulk1 = np.zeros(Nq)
bulk2 = np.zeros(Nq)

### List of edge states
EdgeStates = []
Qedge = []

##### ===========================================================================
for iq in range(Nq):
    ### The synthetic momentum 
    q = q_array[iq]

    ### HL: 1D array of size (Nk) whose elements are the (4,4) Hamiltonians 
    ### Shape of HL = (Nk,4,4)
    ### => SIMD style
    HL = np.array([Hamiltonian(k,q,omega0,v0,U0,omega1,v1,U1,V1) for k in k_array])

    #print(np.shape(HL_array)) 

    # 1D array of size (Nk) of eigvals of HL 
    # Shape of EL = (Nk,4)  
    EL = np.linalg.eigvalsh(HL)
    #print(np.shape(EL))

    ### HR: 1D array of size (Nk) whose elements are the (4,4) Hamiltonians 
    ### Shape of HR = (Nk,4,4)
    ### => SIMD style 
    HR = np.array([Hamiltonian(k,q,omega0,v0,U0,omega2,v2,U2,V2) for k in k_array])

    # 1D array of size (Nk) of eigvals of HR 
    # Shape of ER = (Nk,4)
    ER = np.linalg.eigvalsh(HR)

    ### Calculate the closest bands 
    bulk1[iq] = min(np.amin(EL[:,gap+1]),np.amin(ER[:,gap+1]))
    bulk2[iq] = max(np.amax(EL[:,gap]),np.amax(ER[:,gap]))

    ### Array of energy 
    E_array = np.linspace(bulk2[iq]-epsilonE,bulk1[iq]+epsilonE,NE) 

    #for iE in range(NE):
        # The value of the energy 
    #    E = E_array[iE]

        # The left kvals and kvecs 
    #    kL,WL = keig(E,q,omega0,v0,U0,omega1,v1,U1,V1)

        # The right kvals and kvecs 
    #    kR,WR = keig(E,q,omega0,v0,U0,omega2,v2,U2,V2)

    ### The left kvecs 
    ### This is an 1D array of size (NE) whose elements are (4,4) arrays of kvecs 
    WL = np.array([keig(E,q,omega0,v0,U0,omega1,v1,U1,V1) for E in E_array])

    ### The right kvecs
    WR = np.array([keig(E,q,omega0,v0,U0,omega2,v2,U2,V2) for E in E_array])

    ### Concatenate WL and WR to WW
    ### WW: 1D array of size (NE) whose elements are (4,4) arrays of kvecs 
    WW = np.concatenate((WL[:,:,0:2],-WR[:,:,2:4]),axis=2)

    ### S: 1D array of size (NE) whose elements are determinants of the elements 
    ### of WW 
    S = np.abs(np.linalg.det(WW))

    print(S)
    
    ### Edgestate 
    for iE in range(NE):
        if (S[iE] < epsilon):
            Qedge.append(q)
            EdgeStates.append(E_array[iE])

print(EdgeStates)

##### =================================================================================
##### Plot the figure 
fig,ax = plt.subplots(figsize=(9,12))
ax.plot(q_array,bulk1,color='darkgrey',linewidth=4)
ax.plot(q_array,bulk2,color='darkgrey',linewidth=4)
ax.plot(Qedge,EdgeStates,color='orange',linewidth=4)
ax.set_xlabel(r'q',fontsize=32)
ax.set_ylabel(r'$\omega (2\pi c /\lambda)$',fontsize=32)
plt.show()