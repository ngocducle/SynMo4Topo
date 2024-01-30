import numpy as np
import scipy 
import scipy.linalg as sla 
import matplotlib.pyplot as plt 

### Function to diagonalize the Hamiltonian 
def HamiltonianMpoint(qx, qy, omegaM, v, U, W):
    Hamiltonian = np.array( [ [omegaM + v*(qx+qy)/np.sqrt(2), W, W , U], \
                          [W, omegaM + v*(qx-qy)/np.sqrt(2), U, W],  \
                          [W, U, omegaM + v*(-qx+qy)/np.sqrt(2), W], \
                          [U, W, W, omegaM + v*(-qx-qy)/np.sqrt(2)] ] )

    evalues, evectors = sla.eigh(Hamiltonian) 

    return evalues, evectors 

# Energy at point M 
omegaM = 0.42  

# The group velocity 
v = 1  
print('Group velocity v = '+str(v)) 

# The counter-propagating coupling strength 
U = 0.05   
print('U = '+str(U)) 

# The orthogonally-propagating coupling strength 
W = 0.1   
print('W = '+str(W))

# The array of momenta 
# Gamma' = (-0.2, -0.2)
# X' = (-0.2, 0)
# M  = (0, 0) 

Nk = 50 # number of momenta 
Kmax = 0.05 # maximum value of q from M-point 

qxarray = np.concatenate( (np.linspace(-Kmax, 0, Nk), # Gamma' -> X'
                           np.linspace(0, 0, Nk),    # X' -> M 
                           np.linspace(0, -Kmax, Nk)  # M -> Gamma'  
                          ),
                          axis = 0 ) 

qyarray = np.concatenate( (np.linspace(-Kmax,-Kmax,Nk), # Gamma' -> X'
                           np.linspace(-Kmax, 0, Nk),  # X' -> M  
                           np.linspace(0, -Kmax, Nk)   # M -> Gamma'  
                          ), 
                          axis = 0 )  

# Diagonalize the Hamiltonian 
E1 = np.empty(3*Nk, dtype = complex) 
E2 = np.empty(3*Nk, dtype = complex)
E3 = np.empty(3*Nk, dtype = complex)
E4 = np.empty(3*Nk, dtype = complex)  

index = np.arange(3*Nk) 

for i in range(3*Nk): 
    qx = qxarray[i]  
    qy = qyarray[i]

    evalues, evectors = HamiltonianMpoint(qx, qy, omegaM, v, U, W) 

    E1[i] = evalues[0]
    E2[i] = evalues[1]
    E3[i] = evalues[2]
    E4[i] = evalues[3] 

# Plot the band structure 
fig, ax = plt.subplots() 
plt.plot(index, E1)  
plt.plot(index, E2) 
plt.plot(index, E3)
plt.plot(index, E4) 
plt.vlines(Nk, 0, 1, linestyle = 'dashed', color = 'black') 
plt.vlines(2*Nk,0,1, linestyle = 'dashed', color = 'black') 
plt.ylim(0,0.8) 
tick_locs = [0, Nk, 2*Nk, 3*Nk] 
tick_labs = [ r'$\Gamma \prime$', r'$X\prime$', 'M', r'$\Gamma \prime$' ]
ax.set_xticks(tick_locs) 
ax.set_xticklabels(tick_labs, size = 16) 
#ax.set_xlabel(r'$k a / (2 \pi)$', fontsize = 14) 
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.show() 
