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

##### =======================================================================
##### FUNCTION: PlaneWaves 
def keig(q,E,omega1,v1,U1,omega2,v2,U2,V):
    QH = np.array(
        [
            [omega1-E,U1,V*exp(-1j*np.pi*q),0],
            [U1,omega1-E,0,V*exp(1j*np.pi*q)],
            [V*exp(1j*np.pi*q),0,omega2-E,U2],
            [0,V*exp(-1j*np.pi*q),U2,omega2-E]
        ]
    )

    QS = np.diag([-v1,v1,-v2,v2])

    kvals,kvecs = sla.eig(QH,QS)

    args = np.angle(kvals*exp(1e-9*1j))
    args[args<0] += 2*np.pi 

    GH = kvecs[:,np.argsort(args)]

    return GH/np.linalg.norm(GH,axis=0)

##### Parameters 
omega0 = 0.27815 
v0 = 0.37602 
U0 = 0.02232 

omega1 = 0.27999
v1 = 0.40590 
U1 = 0.02621 
V1 = 0.04511 

omega2 = 0.28010 
v2 = 0.40454 
U2 = 0.02142 
V2 = 0.04728 

### Gap where we calculate the edge states 
gap = 0 

### Criterion for 0 
epsilon = 1e-3 

### Increment of band edge 
epsilonE = 1e-4

### Array of synthetic momenta 
Nq = 501 
q_array = np.linspace(0.45,0.55,Nq)

### Array of genuine momenta 
Nk = 101 
k_array = np.linspace(-0.12,0.12,Nk)

### Number of energy values 
NE = 101 

bulk1 = np.zeros(Nq)
bulk2 = np.zeros(Nq)

EdgeStates = []
Qedge = []

##### =================================================================================
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

    ### Array of energies 
    E_array = np.linspace(bulk2[iq]-epsilonE,bulk1[iq]+epsilonE,NE) 

    ### Arrays of k-eigenstates 
    Wa = np.array([keig(q,E,omega0,v0,U0,omega1,v1,U1,V1) for E in E_array])
    Wb = np.array([keig(q,E,omega0,v0,U0,omega2,v2,U2,V2) for E in E_array])

    WW = np.concatenate((Wb[:,:,0:2],-Wa[:,:,2:4]),axis=2)
    S = np.abs(np.linalg.det(WW))

    ### Scan E_array again 
    for iE in range(NE):
        if (S[iE]<epsilon):
            EdgeStates.append(E_array[iE])
            Qedge.append(q)

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