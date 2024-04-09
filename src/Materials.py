import meep as mp

### ============================================================================== ###
###                                                                                ###
###          This module contains all the materials used in the project            ###
###                                                                                ###
### ============================================================================== ###

Si = mp.Medium(index = 3.54)
SiO2 = mp.Medium(index = 1.46)
PMMA = mp.Medium(index = 1.46)
Dielectric = mp.Medium(epsilon = 12.0)
Air = mp.Medium(epsilon = 1.0)