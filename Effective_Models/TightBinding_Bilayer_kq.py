import numpy as np 
import scipy 
import cmath 
import scipy.linalg as sla 
import matplotlib.pyplot as plt

# The geometrical parameters 
a = 0.2 # lattice parameter (nm)
d = a/2.2 # Size of basis 
alpha = 0 # C4 symmetry breaking parameter 
D = 0.25 # interlayer distance (nm) 
xAtoms = [0,d,d,0] # arrays of x-coordinates of atoms 
yAtoms = [0,0,d,d] # arrays of y-coordiantes of atoms 

# Tight-binding Minh model for monolayer 
t1 = -np.exp(-d**2/a**2)       # intracell hopping integrals along edge 
t2 = -np.exp(-(a-d)**2/a**2)   # intercell hopping integrals 
t3 = -0.5*np.exp(-2*d**2/a**2) # intracell hopping integrals along diagonals 

print('t1 = '+str(t1))
print('t2 = '+str(t2)) 
print('t3 = '+str(t3))

delta_x = 0 # shift (nm)
delta_y = 0 # shift (nm) 

# The list of momenta 
Nk = 21 
k_array = np.linspace(0,0.2,Nk)

Nq = 51 
q_array = np.linspace(0,0.5,Nq)

##### The bilayer Hamiltonian 
def TightBinding_Bilayer_Hamiltonian(t1,t2,t3,alpha,a,d,D,kx,ky,delta_x,delta_y):
    Hamiltonian = np.zeros((8,8),dtype=complex)

    xAtoms = [0,d,d,0]
    yAtoms = [0,0,d,d] 

    # Block (0,0) 
    Hamiltonian[0,0] = 0.0 
    Hamiltonian[0,1] = t1*cmath.exp(1j*kx*d) + t2*cmath.exp(1j*kx*(d-a))
    Hamiltonian[0,2] = t3*(1-alpha)*cmath.exp(1j*(kx*d+ky*d))
    Hamiltonian[0,3] = t1*cmath.exp(1j*ky*d) + t2*cmath.exp(1j*ky*(d-a))

    Hamiltonian[1,0] = np.conjugate(Hamiltonian[0,1])
    Hamiltonian[1,1] = 0.0
    Hamiltonian[1,2] = t1*cmath.exp(1j*ky*d) + t2*cmath.exp(1j*ky*(d-a))
    Hamiltonian[1,3] = t3*(1+alpha)*cmath.exp(1j*(-kx*d+ky*d))

    Hamiltonian[2,0] = np.conjugate(Hamiltonian[0,2])
    Hamiltonian[2,1] = np.conjugate(Hamiltonian[1,2])
    Hamiltonian[2,2] = 0.0 
    Hamiltonian[2,3] = t1*cmath.exp(-1j*kx*d) + t2*cmath.exp(-1j*kx*(d-a))

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
                    disp_x = nx*a + a/2 + xAtoms[j] + delta_x - xAtoms[k] 
                    disp_y = ny*a + a/2 + yAtoms[j] + delta_y - yAtoms[k] 

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
    Hamiltonian[4,5] = t1*cmath.exp(1j*kx*d) + t2*cmath.exp(1j*kx*(d-a))
    Hamiltonian[4,6] = t3*(1-alpha)*cmath.exp(1j*(kx*d+ky*d))
    Hamiltonian[4,7] = t1*cmath.exp(1j*ky*d) + t2*cmath.exp(1j*ky*(d-a))

    Hamiltonian[5,4] = np.conjugate(Hamiltonian[4,5])
    Hamiltonian[5,5] = 0.0 
    Hamiltonian[5,6] = t1*cmath.exp(1j*ky*d) + t2*cmath.exp(1j*ky*(d-a))
    Hamiltonian[5,7] = t3*(1+alpha)*cmath.exp(1j*(-kx*d+ky*d))

    Hamiltonian[6,4] = np.conjugate(Hamiltonian[4,6])
    Hamiltonian[6,5] = np.conjugate(Hamiltonian[5,6])
    Hamiltonian[6,6] = 0.0 
    Hamiltonian[6,7] = t1*cmath.exp(-1j*kx*d) + t2*cmath.exp(-1j*kx*(d-a))

    Hamiltonian[7,4] = np.conjugate(Hamiltonian[4,7])
    Hamiltonian[7,5] = np.conjugate(Hamiltonian[5,7]) 
    Hamiltonian[7,6] = np.conjugate(Hamiltonian[6,7])
    Hamiltonian[7,7] = 0.0 


    return Hamiltonian 
