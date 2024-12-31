import numpy as np
import scipy
import tidy3d as td 
from tidy3d.constants import C_0

##### ============================================================================= #####
###                                                                                   ###
###         This module contains the functions related to the structures              ###
###                  of the FDTD simulations using Tidy3D                             ###
###                                                                                   ###
##### ============================================================================= #####

##### ===================================================================================
##### FUNCTION: define the geometry of the 2D photonic crystal slab monolayer with 
##### square unit cell and rhombus hole. The structure is rotated so that 
##### the source emits the electromagnetic wave along the diagonal direction 
##### that we choose to be the x-direction 
###
###     d: diagonal length of the primitive cell
###     b: edge length of undeformed square hole
###     h: thickness of the slab
###     e: deformation parameter 
###     mat_envir: the material forming the environment
###     mat_slab: the material forming the slab
###     Ncell: the number of unit cells forming the slab

def td_2DSlab1L_RHoleP(d,b,h,e,mat_envir,mat_slab,Ncell):

    ### The environment
    envir = td.Structure(
        geometry = td.Box(
            center = (0,0,0),
            size = (td.inf,td.inf,td.inf),
            ),
        medium = mat_envir,
        name = 'envir',
        )

    ### Define the slab
    slab = td.Structure(
        geometry = td.Box(
            center = (0,0,0),
            size = (td.inf,td.inf,h),
            ),
        medium = mat_slab,
        name = 'slab',
        )

    sim_structures = [envir,slab]

    ##### We divide 2 cases: Ncell is odd and Ncell is even
    if Ncell%2 == 1:
        Nhalf = int((Ncell-1)/2)
        print('Nhalf = '+str(Nhalf))

        # Central line
        for j in np.arange(-Nhalf,Nhalf+1):
            holename = 'holec_'+str(j)

            vertices = np.array(
                [
                    (j*d + b*(1+e)/(1-e)/np.sqrt(2),0),
                    (j*d, b*(1-e)/(1+e)/np.sqrt(2)),
                    (j*d - b*(1+e)/(1-e)/np.sqrt(2),0),
                    (j*d, -b*(1-e)/(1+e)/np.sqrt(2))
                ]
            )

            hole = td.Structure(
                geometry = td.PolySlab(
                    axis = 2,
                    reference_plane = 'middle',
                    slab_bounds = [-0.5*h,0.5*h],
                    vertices = vertices,
                    ),
                medium = mat_envir,
                name = holename,
                )

            sim_structures.append(hole)

        # Upper line
        for j in np.arange(-Nhalf,Nhalf):
            holename = 'holeu_'+str(j)

            vertices = np.array(
                [
                    ( (j+0.5)*d + b*(1+e)/(1-e)/np.sqrt(2),0.5*d),
                    ( (j+0.5)*d, b*(1-e)/(1+e)/np.sqrt(2) +0.5*d),
                    ( (j+0.5)*d - b*(1+e)/(1-e)/np.sqrt(2),0.5*d),
                    ( (j+0.5)*d, -b*(1-e)/(1+e)/np.sqrt(2)+0.5*d)
                ]
            )

            hole = td.Structure(
                geometry = td.PolySlab(
                    axis = 2,
                    reference_plane = 'middle',
                    slab_bounds = [-0.5*h,0.5*h],
                    vertices = vertices,
                    ),
                medium = mat_envir,
                name = holename,
                )

            sim_structures.append(hole)

        # Lower line 
        for j in np.arange(-Nhalf,Nhalf):
            holename = 'holel_'+str(j)

            vertices = np.array(
                [
                    ( (j+0.5)*d + b*(1+e)/(1-e)/np.sqrt(2),-0.5*d),
                    ( (j+0.5)*d, b*(1-e)/(1+e)/np.sqrt(2) -0.5*d),
                    ( (j+0.5)*d - b*(1+e)/(1-e)/np.sqrt(2),-0.5*d),
                    ( (j+0.5)*d, -b*(1-e)/(1+e)/np.sqrt(2)-0.5*d)
                ]
            )

            hole = td.Structure(
                geometry = td.PolySlab(
                    axis = 2,
                    reference_plane = 'middle',
                    slab_bounds = [-0.5*h,0.5*h],
                    vertices = vertices,
                    ),
                medium = mat_envir,
                name = holename,
            )

            sim_structures.append(hole)

    else:
        Nhalf = int(Ncell/2)
        print('Nhalf = '+str(Nhalf))

        # Central line
        for j in np.arange(-Nhalf,Nhalf):
            holename = 'holec_'+str(j)

            vertices = np.array(
                [
                    ((j+0.5)*d + b*(1+e)/(1-e)/np.sqrt(2),0),
                    ((j+0.5)*d, b*(1-e)/(1+e)/np.sqrt(2)),
                    ((j+0.5)*d - b*(1+e)/(1-e)/np.sqrt(2),0),
                    ((j+0.5)*d, -b*(1-e)/(1+e)/np.sqrt(2))
                ]
            )

            hole = td.Structure(
                geometry = td.PolySlab(
                    axis = 2,
                    reference_plane = 'middle',
                    slab_bounds = [-0.5*h,0.5*h],
                    vertices = vertices,
                    ),
                medium = mat_envir,
                name = holename,
                )

            sim_structures.append(hole)

        # Upper line
        for j in np.arange(-Nhalf,Nhalf):
            holename = 'holeu_'+str(j)

            vertices = np.array(
                [
                    ( j*d + b*(1+e)/(1-e)/np.sqrt(2),0.5*d),
                    ( j*d, b*(1-e)/(1+e)/np.sqrt(2) +0.5*d),
                    ( j*d - b*(1+e)/(1-e)/np.sqrt(2),0.5*d),
                    ( j*d, -b*(1-e)/(1+e)/np.sqrt(2)+0.5*d)
                ]
            )

            hole = td.Structure(
                geometry = td.PolySlab(
                    axis = 2,
                    reference_plane = 'middle',
                    slab_bounds = [-0.5*h,0.5*h],
                    vertices = vertices,
                    ),
                medium = mat_envir,
                name = holename,
                )

            sim_structures.append(hole)

        # Lower line
        for j in np.arange(-Nhalf,Nhalf):
            holename = 'holel_'+str(j)

            vertices = np.array(
                [
                    ( j*d + b*(1+e)/(1-e)/np.sqrt(2),-0.5*d),
                    ( j*d, b*(1-e)/(1+e)/np.sqrt(2) -0.5*d),
                    ( j*d - b*(1+e)/(1-e)/np.sqrt(2),-0.5*d),
                    ( j*d, -b*(1-e)/(1+e)/np.sqrt(2)-0.5*d)
                ]
            )

            hole = td.Structure(
                geometry = td.PolySlab(
                    axis = 2,
                    reference_plane = 'middle',
                    slab_bounds = [-0.5*h,0.5*h],
                    vertices = vertices,
                    ),
                medium = mat_envir,
                name = holename,
                )

            sim_structures.append(hole)

    return sim_structures



##### ===================================================================================
##### FUNCTION: define the geometry of the heterojunction between two
##### 2D photonic crystal slab monolayers with 
##### square unit cell and rhombus hole. The structure is rotated so that 
##### the source emits the electromagnetic wave along the diagonal direction 
##### that we choose to be the x-direction 
##### ATTENTION! The two photonic crystal slabs are made of the SAME material

def td_2DSlab1L_RHoleP_hj(d,h,b1,e1,b2,e2,dist,mat_envir,mat_slab,Ncell):

    ### The environment
    envir = td.Structure(
        geometry = td.Box(
            center = (0,0,0),
            size = (td.inf,td.inf,td.inf),
            ),
        medium = mat_envir,
        name = 'envir',
        )
    
    ### Define the slab
    slab_up = td.Structure(
        geometry = td.Box(
            center = (0,0,0.5*(h+dist)),
            size = (td.inf,td.inf,h),
            ),
        medium = mat_slab,
        name = 'slab_up',
        )

    sim_structures = [envir,slab_up]

    return sim_structures