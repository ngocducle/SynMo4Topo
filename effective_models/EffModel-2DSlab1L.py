import numpy as np 
import scipy 
import scipy.linalg as sla
import matplotlib.pyplot as plt

### ================================================================================ ###
#                                                                                      #
#           Fit the effective model for 2D photonic crystal slab monolayer             #
#                                                                                      #
### ================================================================================ ###

### FUNCTION: Hamiltonian of 2D photonic crystal slab monolayer at M-point
def EffModel_2DSlab1L-M(qx,qy,omegaM,v,U,W):
    Hamiltonian = np.array(
        [
            [omegaM+v*(qx+qy)/np.sqrt(2),W,W,U],
            [W,omegaM+v(qx-qy)/np.sqrt(2),U,W],
            [W,U,omegaM+v*(-qx+qy)/np.sqrt(2),W],
            [U,W,W,omegaM+v*(-qx-qy)/np.sqrt(2)]
        ]
    )

    evalues, evectors = sla.eigh(Hamiltonian)

    return evalues, evectors

