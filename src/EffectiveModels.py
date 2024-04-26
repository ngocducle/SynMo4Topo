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