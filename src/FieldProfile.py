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

    ### Get the epsilon profile
    resolution_eps = 81 # Number of pixels per a 
    num_periods = 3 # Number of periods along each direction 
    md = mpb.MPBData(rectify = True, 
                     periods = num_periods,  
                     resolution = resolution_eps)
    eps = ModeSolver.get_epsilon()
    converted_eps = md.convert(eps)

    print('The shape of converted_eps: '+str(np.shape(converted_eps)))

    ### ATTENTION! Be careful that the structure is also copied num_periods
    # time along the z-axis 

    ### The epsilon profile in the Oxy plane (parallel to the plane Oxy: z = 0)
    # Array of zvalues to plot eps
    # The shape of converted_eps is 
    # (num_periods*resolution_eps*length_unit_cell_x,
    # num_periods*resolution_eps*length_unit_cell_y,
    # num_periods*resolution_eps*length_unit_cell_z)
    # Here:
    # length_unit_cell_x = length_unit_cell_y = 1
    # length_unit_cell_z = Lz 
    shape_eps = np.shape(converted_eps)
    len_zarray_eps = int(shape_eps[2] / num_periods)
    print('len_zarray_eps = '+str(len_zarray_eps)) 
    zarray_eps = np.linspace(-0.5*Lz,0.5*Lz,len_zarray_eps) 

    # We look for the index of the element of zarray_eps which is closest
    # to zvalue 
    for i in range(len_zarray_eps-1): 
        if ((zarray_eps[i] <= zvalue) and (zvalue < zarray_eps[i+1])):
            zindex_eps = i 
            break 
        else:
            zindex_eps = len_zarray_eps - 1

    # The index of the slice in the z-axis
    print('zindex_eps = '+str(zindex_eps))