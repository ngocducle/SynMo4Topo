import numpy as np
import scipy
import scipy.linalg as sla
import cmath 

########################################################################################
#                                                                                      # 
#       This module contains the Hamiltonians for different effective models           #
#                                                                                      #
########################################################################################


### FUNCTION: Hamiltonian of 2D photonic crystal slab monolayer at M-point
def EffModel_2DSlab1L_M(qx,qy,omegaM,v,U,W):
    Hamiltonian = np.array(
        [
            [omegaM+v*(qx+qy)/np.sqrt(2),W,W,U],
            [W,omegaM+v*(qx-qy)/np.sqrt(2),U,W],
            [W,U,omegaM+v*(-qx+qy)/np.sqrt(2),W],
            [U,W,W,omegaM+v*(-qx-qy)/np.sqrt(2)]
        ]
    )

    evalues, evectors = sla.eigh(Hamiltonian)

    return evalues, evectors

### FUNCTION: Hamiltonian of 2D photonic crystal slab monolayer at M-point
def EffModel_2DSlab2L_M(qx,qy,deltax,deltay,omegaM,v,U,W,V):
    omegapp = omegaM + v*(qx+qy)/np.sqrt(2) 
    omegapm = omegaM + v*(qx-qy)/np.sqrt(2)
    omegamp = omegaM + v*(-qx+qy)/np.sqrt(2)
    omegamm = omegaM + v*(-qx-qy)/np.sqrt(2)

    K = 2.0*np.pi

    Vpp = V * cmath.exp(- 0.5 * 1j * (deltax+deltay) * K)
    Vpm = V * cmath.exp(- 0.5 * 1j * (deltax-deltay) * K)
    Vmp = V * cmath.exp(- 0.5 * 1j * (-deltax+deltay) * K)
    Vmm = V * cmath.exp(- 0.5 * 1j * (-deltax-deltay) * K)

    Hamiltonian = np.array(
        [
        [omegapp, W, W, U, Vpp, 0, 0, 0],
        [W, omegapm, U, W, 0, Vpm, 0, 0],
        [W, U, omegamp, W, 0, 0, Vmp, 0],
        [U, W, W, omegamm, 0, 0, 0, Vmm],
        [Vmm, 0, 0, 0, omegapp, W, W, U],
        [0, Vmp, 0, 0, W, omegapm, U, W],
        [0, 0, Vpm, 0, W, U, omegamp, W],
        [0, 0, 0, Vpp, U, W, W, omegamm]
        ]
    )

    evalues, evectors = sla.eigh(Hamiltonian)

    return evalues, evectors 

### FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer with rhombus hole
def EffHam_2DSlab2L_RHole_M(kx,ky,qx,qy,
                            omega1,v1,U1,W1,alpha1,
                            omega2,v2,U2,W2,alpha2,V):
    Hamiltonian = np.zeros((8,8),dtype=complex)

    K = np.pi 

    ### Block (1,1)
    Hamiltonian[0,0] = omega1 + v1*(kx+ky)/np.sqrt(2.0)
    Hamiltonian[0,1] = W1 
    Hamiltonian[0,2] = W1 
    Hamiltonian[0,3] = U1*(1+alpha1) 

    Hamiltonian[1,0] = W1 
    Hamiltonian[1,1] = omega1 + v1*(kx-ky)/np.sqrt(2.0)
    Hamiltonian[1,2] = U1*(1-alpha1) 
    Hamiltonian[1,3] = W1 

    Hamiltonian[2,0] = W1 
    Hamiltonian[2,1] = U1*(1-alpha1) 
    Hamiltonian[2,2] = omega1 + v1*(-kx+ky)/np.sqrt(2.0)
    Hamiltonian[2,3] = W1 

    Hamiltonian[3,0] = U1*(1+alpha1) 
    Hamiltonian[3,1] = W1 
    Hamiltonian[3,2] = W1 
    Hamiltonian[3,3] = omega1 + v1*(-kx-ky)/np.sqrt(2.0) 

    ### Block (1,2)
    Hamiltonian[0,4] = V*cmath.exp(-1j*(qx+qy)*K)
    Hamiltonian[1,5] = V*cmath.exp(-1j*(qx-qy)*K)
    Hamiltonian[2,6] = V*cmath.exp(-1j*(-qx+qy)*K)
    Hamiltonian[3,7] = V*cmath.exp(-1j*(-qx-qy)*K) 

    ### Block (2,1)
    Hamiltonian[4,0] = V*cmath.exp(1j*(qx+qy)*K)
    Hamiltonian[5,1] = V*cmath.exp(1j*(qx-qy)*K)
    Hamiltonian[6,2] = V*cmath.exp(1j*(-qx+qy)*K)
    Hamiltonian[7,3] = V*cmath.exp(1j*(-qx-qy)*K)

    ### Block (2,2)
    Hamiltonian[4,4] = omega2 + v2*(kx+ky)/np.sqrt(2.0)
    Hamiltonian[4,5] = W2 
    Hamiltonian[4,6] = W2 
    Hamiltonian[4,7] = U2*(1+alpha2) 

    Hamiltonian[5,4] = W2 
    Hamiltonian[5,5] = omega2 + v2*(kx-ky)/np.sqrt(2.0)
    Hamiltonian[5,6] = U2*(1-alpha2) 
    Hamiltonian[5,7] = W2 

    Hamiltonian[6,4] = W2 
    Hamiltonian[6,5] = U2*(1-alpha2) 
    Hamiltonian[6,6] = omega2 + v2*(-kx+ky)/np.sqrt(2.0)
    Hamiltonian[6,7] = W2 

    Hamiltonian[7,4] = U2*(1+alpha2) 
    Hamiltonian[7,5] = W2 
    Hamiltonian[7,6] = W2 
    Hamiltonian[7,7] = omega2 + v2*(-kx-ky)/np.sqrt(2.0)

    return Hamiltonian 