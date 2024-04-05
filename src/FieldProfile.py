import numpy as np
import scipy 
import meep as mp
from meep import mpb

### ============================================================================== ###
###                                                                                ###
### This module contains all the functions to calculate the field profiles         ###
###                                                                                ###
### ============================================================================== ###

##### FUNCTION: Get the E-fields from a mode solver
def get_fields(ms,band,Bloch_Phase):
    efields.append(ms.get_efield(band,bloch_phase=Bloch_Phase))

##### FUNCTION: Calculate the E-fields
def EField_Profile(ModeSolver,k_field,Bloch_Phase):
    ### Initiate the array of E-fields
    efields = []


