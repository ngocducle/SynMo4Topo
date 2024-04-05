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
def EField_Profile(ModeSolver,k_field,Lz,zvalue,polarization,
                   resolution_eps,resolution_field,
                   num_periods,Bloch_Phase):
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
    eps_Oxy = DielectricProfileZvalue(ModeSolver,zvalue,Lz,resolution_eps,num_periods)

    # The meshgrid of (x,y)
    Nx = resolution_eps*num_periods 
    Ny = resolution_eps*num_periods 
    Xlim = 0.5*num_periods   
    Ylim = 0.5*num_periods  
    X, Y = np.meshgrid( np.linspace(-Xlim,Xlim,Nx), 
                        np.linspace(-Ylim,Ylim,Ny) ) 
    
    ### Get the fields as a MPBArray
    EField = mpb.MPBData(rectify = True, 
                                resolution = resolution_field, 
                                periods = num_periods) 
    
    return efields,EField

##### FUNCTION: Plot the field profile
def ExtractEField_Profile(efields,EField,Lz,zvalue,
                      resolution,resolution_eps,resolution_field,num_periods):
    # Define the arrays for X and Y to plot the fields
    Xlim = 0.5*num_periods   
    Ylim = 0.5*num_periods  
    Nx = resolution_field*num_periods
    Ny = resolution_field*num_periods

    Xfield, Yfield = np.meshgrid( 
    np.linspace(-Xlim, Xlim, Nx),
    np.linspace(-Ylim, Ylim, Ny) 
    )

    # Index for z: the z-coordinate of the slice 
    # Array of the z-coordinate for field pattern (different from array of 
    # dielectric constant)
    zarray_field_len = int(Lz*resolution.z) 
    zarray_field = np.linspace(-0.5*Lz,0.5*Lz,zarray_field_len)  

    for i in range(zarray_field_len-1):
        if (zarray_field[i]<=zvalue) and (zvalue<zarray_field[i+1]):
            zindex_field = i 
            break
        else:
            zindex_field = zarray_field_len - 1 
  
    print('zindex_field = '+str(zindex_field))  

    # efields and hfields contain num_bands array (datasheets?), 
    # each corresponds to one band. 
    # Each array has 4 indices (x,y,z,E-mp.Vector3)
    # First, swap the axes y and z, i.e. swapping indices 1 and 2 
    # We need to swap to move the indices of the y-coordinate and 
    # of the field to the two last positions, so that the other indices
    # are at the leftmost side to be able to be replaced by the ellipsis.
    # This is because we can put only 1 ellipsis into the formula. 
    Efieldx = []
    Efieldy = []
    Efieldz = [] 

    for f in efields:
        # Get the x component of the E-fields 
        print('The shape of f: '+str(np.shape(f))) 

        # Take the slice (band,x,y) for Ex at z
        Ex = f[:,:,zindex_field,0]  

        # Save the data Ex to Efieldx
        Efieldx.append(EField.convert(Ex))  

        # Get the y component of the E-fields 
        print('The shape of f: '+str(np.shape(f))) 

        # Take the slice (band,x,y) for Ey at z
        Ey = f[:,:,zindex_field,1]   

        # Save the data Ey to converted
        Efieldy.append(EField.convert(Ey))  

        # Get the z component of the E-fields 
        print('The shape of f: '+str(np.shape(f))) 

        # Take the slice (band,x,y) for Ez at z
        Ez = f[:,:,zindex_field,2] 

        # Save the data Ez to converted
        Efieldz.append(EField.convert(Ez))  

    return Efieldx,Efieldy,Efieldz