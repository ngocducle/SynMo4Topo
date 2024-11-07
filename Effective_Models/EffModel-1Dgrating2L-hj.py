import numpy as np 
import scipy 
import scipy.linalg as sla
import matplotlib.pyplot as plt 
from cmath import exp 

##### =========================================================================
##### FUNCTION: Hamiltonian 
def Hamiltonian(k,q,omega0,v,U,Delta,V):
    H = np.array([[omega0+v*k,U+Delta,V*exp(-1j*np.pi*q),0],
                     [U+Delta,omega0-v*k,0,V*exp(1j*np.pi*q)],
                     [V*exp(1j*np.pi*q),0,omega0+v*k,U-Delta],
                     [0,V*exp(-1j*np.pi*q),U-Delta,omega0-v*k]])

    return H 

##### =========================================================================-
##### FUNCTION: keig 
def keig(E,q,omega0,v,U,Delta,V):
    H0 = np.zeros((4,4),dtype=complex)
    H0[0,0] = omega0 - E 
    H0[0,1] = U+Delta 
    H0[0,2] = V*exp(-1j*np.pi*q)
    H0[0,3] = 0 
    H0[1,1] = omega0 - E 
    H0[1,2] = 0 
    H0[1,3] = V*exp(1j*np.pi*q)
    H0[2,2] = omega0 - E 
    H0[2,3] = U-Delta
    H0[3,3] = omega0 - E 

    for i in range(1,4):
        for j in range(i):
            H0[i,j] = np.conj(H0[j,i])
    
    H1 = np.diag([-v,v,-v,v])

    kvals,kvecs = sla.eig(H0,H1)

    #print(kvals)
    #print(kvecs)

    ### Move the arguments to the range 0 <= argument <= 2*pi 
    arg_array = np.angle(kvals)+1e-9*1j
    arg_array[arg_array<0] += 2*np.pi 

    #print(arg_array)

    #print(np.argsort(arg_array))

    kvals = kvals[np.argsort(arg_array)]
    #print(kvals)

    kvecs = kvecs[:,np.argsort(arg_array)]
    kvecs = kvecs/np.linalg.norm(kvecs,axis=0)
    #print(kvecs)

    return kvals,kvecs 

##### Parameters 
omega0 = 0.25 
v = 0.30 
U = 0.02 
V = 0.05 
Delta = 0.2*U 

### Array of k 
Kmax = 0.05 
Nk = 101 
k_array = np.linspace(-Kmax,Kmax,Nk)

### Array of q 
Qmax = 1.0 
Nq = 101 
q_array = np.linspace(-Qmax,Qmax,Nq)

### The gap we calculate the edge states 
gap = 0 

### Criterion for 0 
epsilon = 1e-3 

### Number of E values to scan 
NE = 501 

### Small increment in band edge 
epsilonE = 1e-4 

### Initialize the edge states to be empty 
edge_state = []

### The arrays of obstructed and bulk bands 
allEmax = np.zeros(Nq)
allEmin = np.zeros(Nq)
bulk1 = np.zeros(Nq)
bulk2 = np.zeros(Nq)

S_array = np.zeros((NE,Nq))

##### ==================================================================================
##### Scan over the q_array 
for iq in range(Nq): 
    ### The synthetic momentum 
    q = q_array[iq]

    ### Calculate the bulk band structure 
    EL = np.zeros((Nk,4))
    ER = np.zeros((Nk,4))

    for ik in range(Nk):
        ### Genuine momentum 
        k = k_array[ik]

        ### The bulk band structure 
        ### The left-hand Hamiltonian 
        eigvals = sla.eigvalsh(Hamiltonian(k,q,omega0,v,U,Delta,V))
        EL[ik,:] = eigvals 

        ### The right-hand Hamiltonian 
        eigvals = sla.eigvalsh(Hamiltonian(k,q,omega0,v,U,-Delta,V))
        ER[ik,:] = eigvals 

    #HL = [Hamiltonian(k,q,omega0,v,U,Delta,V) for k in k_array]
    #print(np.shape(HL))

    # Left-hand Hamiltonian eigenvalues: shape(EL) = (Nk,4)
    #EL = sla.eigvalsh([Hamiltonian(k,q,omega0,v,U,Delta,V) for k in k_array])
    #EL = eigvals 

    # Right-hand Hamiltonian eigenvalues: shape(ER) = (Nk,4)
    #ER = sla.eigvalsh([Hamiltonian(k,q,omega0,v,U,-Delta,V) for k in k_array]) 
    #ER = eigvals  

    ### Calculate the obstructed bands 
    #allEmax[iq] = min(min(EL[:,gap+1]),min(ER[:,gap+1]))
    #allEmin[iq] = max(max(EL[:,gap]),max(ER[:,gap]))

    ### Calculate the bulk bands 
    bulk1[iq] = min(np.amin(EL[:,gap+1]),np.amin(ER[:,gap+1]))
    bulk2[iq] = max(np.amax(EL[:,gap]),np.amax(ER[:,gap]))

    ##### ==============================================================================
    ##### Calculate the edge states at the synthetic momentum q 
    ### Array of energies 
    E_array = np.linspace(bulk2[iq]-epsilonE,bulk1[iq]+epsilonE,NE)
    #print(E_array) 

    #S_array = np.zeros(NE)

    ### Scan the energy array E_array 
    for iE in range(NE):
        # The energy 
        E = E_array[iE]
        print(E)

        # The matrix of eigenstates 
        W = np.zeros((4,4))

        # The left-hand Hamiltonian 
        print('Left-hand Hamiltonian')
        kL,WL = keig(E,q,omega0,v,U,Delta,V)

        #print(WL)

        # The right-hand Hamiltonian 
        print('Right-hand Hamiltonian')
        kR,WR = keig(E,q,omega0,v,U,-Delta,V)

        ### Combine WL and WR to the matrix of eigenstates 
        W[:,0] = WL[:,2]
        W[:,1] = WL[:,3]
        W[:,2] = -WR[:,0]
        W[:,3] = -WR[:,1]

        ### The determinant of W 
        S = np.abs(np.linalg.det(W))

        S_array[iE,iq] = S 

        print('S = '+str(S))

        print(E)

    ### If S == 0 then add to the list of edge states 
    #for iE in range(1,NE-1):
    #    if ((S_array[iE] < epsilon) and (S_array[iE-1]>=S_array[iE]) and (S_array[iE+1]>=S_array[iE])):
    #        edge_state.append([q,E])

#edge_state = np.array(edge_state)
#print(edge_state)

##### Plot the figure 
fig,ax = plt.subplots(figsize=(8,10))
#ax.plot(q_array,allEmax)
#ax.plot(q_array,allEmin)
ax.plot(q_array,bulk1)
ax.plot(q_array,bulk2)
#ax.plot(edge_state[:,0],edge_state[:,1],'o',markerfacecolor='red',markeredgecolor='red')
ax.set_xlabel('q',fontsize=16)
ax.set_ylabel(r'$\omega$',fontsize=16)
plt.show()

##### Plot the map of S 
fig,ax = plt.subplots(figsize=(8,10))
ax.pcolormesh(S_array)
plt.show()
