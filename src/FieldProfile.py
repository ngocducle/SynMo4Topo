import numpy as np
import scipy 
import meep as mp
from meep import mpb

import sys
sys.path.insert(0,'.')
from DielectricProfile import DielectricProfileZvalue

### ============================================================================== ###
###                                                                                ###
### This module contains all the functions to calculate the field profiles         ###
###                                                                                ###
### ============================================================================== ###

##### FUNCTION: Calculate the E-fields
def EField_Profile(ModeSolver,k_field,Lz,zvalue,polarization,Bloch_Phase):
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

    ### The dielectric profile in the plane z = zvalue
    resolution_eps = 81 # Number of pixels per a 
    num_periods = 3 # Number of periods along each direction 

    eps_Oxy = DielectricProfileZvalue(ModeSolver,zvalue,Lz,resolution_eps,num_periods)

    # The meshgrid of (x,y)
    Nx = resolution_eps*num_periods 
    Ny = resolution_eps*num_periods 
    Xlim = 0.5*num_periods   
    Ylim = 0.5*num_periods  
    X, Y = np.meshgrid( np.linspace(-Xlim,Xlim,Nx), 
                        np.linspace(-Ylim,Ylim,Ny) ) 
    
    ### Get the fields as a MPBArray
    resolution_field = 81 
    ElectricField = mpb.MPBData(rectify = True, 
                                resolution = resolution_field, 
                                periods = num_periods) 
    
    return ElectricField
