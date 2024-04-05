import numpy as np
import scipy 
import meep as mp
from meep import mpb

### ============================================================================== ###
###                                                                                ###
### This module contains all the functions to calculate the field profiles         ###
###                                                                                ###
### ============================================================================== ###

##### FUNCTION: Calculate the E-fields
def EField_Profile(ModeSolver,k_field,polarization,Bloch_Phase):
    ### Initiate the array of E-fields
    efields = []

    ##### FUNCTION: Get the E-fields from a mode solver
    if Bloch_Phase == 'True':
        def get_e_fields(ModeSolver,band):
            efields.append(ModeSolver.get_efield(band,bloch_phase=True))
    else:
        def get_e_fields(ModeSolver,band):
            efields.append(ModeSolver.get_efield(band,bloch_phase=False))

    ### Rerun the simulation with respect to the polarization
    if polarization == 'all':
        ModeSolver.run(
            mpb.output_at_kpoint(
                         k_field,
                         mpb.fix_efield_phase,
                         get_e_fields
                        ) 
        )
    elif polarization == 'zeven':
        ModeSolver.run_zeven(
            mpb.output_at_kpoint(
                         k_field,
                         mpb.fix_efield_phase,
                         get_e_fields
                        ) 
        )
    elif polarization == 'zodd':
        ModeSolver.run_zodd(
            mpb.output_at_kpoint(
                         k_field,
                         mpb.fix_efield_phase,
                         get_e_fields
                        ) 
        )
    else:
        print('ERROR! The polarization does not belong to the allowed list')
        exit()

    


