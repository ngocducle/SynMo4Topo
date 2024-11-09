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
    H0 = np.array(
        [
            [omega1-E,U1,V*exp(-1j*np.pi*q),0],
            [U1,omega1-E,0,V*exp(1j*np.pi*q)],
            [V*exp(1j*np.pi*q),0,omega2-E,U2],
            [0,V*exp(-1j*np.pi*q),U2,omega2-E]
        ]
    )

    H1 = np.diag([-v1,v1,-v2,v2])

    #print('H0 = ')
    #print(H0)

    #print('H1 = ')
    #print(H1)

    kvals,kvecs = sla.eig(H0,H1)

    #print('kvals = ')
    #print(kvals)
    #print('kvecs = ')
    #print(kvecs)

    args = np.angle(kvals*exp(1e-9*1j))
    args[args<0] += 2*np.pi 

    GH = kvecs[:,np.argsort(args)]

    #print('argsort = ')
    #print(np.argsort(args))

    #print('arg_list = ')
    #print(args[np.argsort(args)])

    #print('kvecs = ')
    #print(GH/np.linalg.norm(GH,axis=0))

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
epsilonE = 1e-5 

### Array of synthetic momenta 
Nq = 1001 
q_array = np.linspace(0.0,1.0,Nq)

### Array of genuine momenta 
Nk = 101 
k_array = np.linspace(-0.12,0.12,Nk)

### Number of energy values 
NE = 251 

EL = np.zeros((Nk,4))
ER = np.zeros((Nk,4))

bulk1 = np.zeros(Nq)
bulk2 = np.zeros(Nq)

EdgeStates = []
Qedge = []

count = 1 

##### =================================================================================
for iq in range(Nq):

    ### The synthetic momentum 
    q = q_array[iq]

    ### Scan the k_array 
    for ik in range(Nk):
        ### The genuine momentum 
        k = k_array[ik]

        ### The left-hand Hamiltonian 
        HL = Hamiltonian(k,q,omega0,v0,U0,omega1,v1,U1,V1) 

        ### Diagonalize the left-hand Hamiltonian 
        EL[ik,:] = np.linalg.eigvalsh(HL)

        ### The right-hand Hamiltonian 
        HR = Hamiltonian(k,q,omega0,v0,U0,omega2,v2,U2,V2)

        ### Diagonalize the right-hand Hamiltonian 
        ER[ik,:] = np.linalg.eigvalsh(HR)

    ### Calculate the closest bands 
    bulk1[iq] = min(np.amin(EL[:,gap+1]),np.amin(ER[:,gap+1]))
    bulk2[iq] = max(np.amax(EL[:,gap]),np.amax(ER[:,gap]))

    ### Array of energies 
    E_array = np.linspace(bulk2[iq]-epsilonE,bulk1[iq]+epsilonE,NE) 

    ### Array of determinants 
    S = np.zeros(NE)

    for iE in range(NE):
        #print('count = '+str(count))

        # The value of the energy 
        E = E_array[iE]

        # The left-hand side 
        Wa = keig(q,E,omega0,v0,U0,omega1,v1,U1,V1)

        # The right-hand side 
        Wb = keig(q,E,omega0,v0,U0,omega2,v2,U2,V2)

        # Concatentate to the determinant 
        WW = np.concatenate((Wb[:,0:2],-Wa[:,2:4]),axis=1)

        #print('WW = ')
        #print(WW)

        # The determinant 
        S[iE] = np.abs(np.linalg.det(WW))
        #print('q = '+str(q))
        #print('k = '+str(k))
        #print('E = '+str(E))
        #print('S = '+str(S[iE]))
        print(S[iE])

        count = count + 1

    #print(S)

    ### Scan E_array again 
    for iE in range(NE):
        if (S[iE]<epsilon):
            EdgeStates.append(E_array[iE])
            Qedge.append(q)

#print(str(i)+','+str(S[i]) for i in range(NE))    

##### =================================================================================
##### Plot the figure 
fig,ax = plt.subplots(figsize=(9,12))
ax.plot(q_array,bulk1,color='darkgrey',linewidth=4)
ax.plot(q_array,bulk2,color='darkgrey',linewidth=4)
#ax.plot(Qedge,EdgeStates,color='orange',linewidth=4)
ax.plot(Qedge,EdgeStates,'o',markerfacecolor='red',markeredgecolor='red')
ax.set_xlabel(r'q',fontsize=32)
ax.set_ylabel(r'$\omega (2\pi c /\lambda)$',fontsize=32)
plt.show()