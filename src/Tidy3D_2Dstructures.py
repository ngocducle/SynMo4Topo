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
##### 2D photonic crystal slab bilayers with 
##### square unit cell and rhombus hole. The structure is rotated so that 
##### the source emits the electromagnetic wave along the diagonal direction 
##### that we choose to be the x-direction 
##### ATTENTION! The two photonic crystal slabs are made of the SAME material
###
###     d: diagonal of unit cell 
###     h: thickness of each slab
###     b1: edge length of undeformed square hole type 1 
###     e1: deformation parameter of undeformed square hole type 1
###     b2: edge length of undeformed square hole type 2
###     e2: deformation parameter of undeformed square hole type 2
###     Lx: size of the structure along x-axis
###     pad: length of the pad at the two sides
###     dist: distance separating the two slabs
###     delta: relative displacement between the two slabs
###     mat_envir: materials of the environment
###     mat_slab: materials of the slab
###     Ncell: number of unit cells at each side 

def td_2DSlab2L_RHoleP_hj(d,h,b1,e1,b2,e2,Lx,pad,dist,delta,mat_envir,mat_slab,Ncell):

    ### The environment
    envir = td.Structure(
        geometry = td.Box(
            center = (0,0,0),
            size = (td.inf,td.inf,td.inf),
            ),
        medium = mat_envir,
        name = 'envir',
        )
    
    ### Slab 
    slab = td.Structure(
        geometry = td.Box(
            center = (0,0,0),
            size = (td.inf,td.inf,2*h+dist),
        ),
        medium = mat_slab,
        name = 'slab',
    )
    
    ### Air space
    air_space = td.Structure(
        geometry = td.Box(
            center = (0,0,0),
            size = (Lx-2*pad,td.inf,dist)
        ),
        medium = mat_envir,
        name = 'air_space',
    )
    
    """### Define the pads
    pad_L = td.Structure(
        geometry = td.Box(
            center = (-0.5*Lx,0,0),
            size = (2*pad,td.inf,2*h+dist),
        ),
        medium = mat_slab,
        name = 'pad_L',
    )

    pad_R = td.Structure(
        geometry = td.Box(
            center = (0.5*Lx,0,0),
            size = (2*pad,td.inf,2*h+dist),
        ),
        medium = mat_slab,
        name = 'pad_R',
    )
    
    ### Define the upper slab
    slab_up = td.Structure(
        geometry = td.Box(
            center = (0,0,0.5*(h+dist)),
            size = (td.inf,td.inf,h),
            ),
        medium = mat_slab,
        name = 'slab_up',
        )
    
    ### Define the lower slab
    slab_down = td.Structure(
        geometry = td.Box(
            center = (0,0,-0.5*(h+dist)),
            size = (td.inf,td.inf,h),
            ),
        medium = mat_slab,
        name = 'slab_down',
        )"""

    #sim_structures = [envir,pad_L,pad_R,slab_up,slab_down]
    sim_structures = [envir,slab,air_space]

    ##### =============================================================================
    ##### The top slab
    ### The right-hand side slab
    # Central line 
    for j in np.arange(0,Ncell):
        holename = 'holetrc_'+str(j)

        vertices = np.array(
            [
                (0.5*delta*d+0.25*d+j*d + b2*(1+e2)/(1-e2)/np.sqrt(2),0),
                (0.5*delta*d+0.25*d+j*d,  b2*(1-e2)/(1+e2)/np.sqrt(2)),
                (0.5*delta*d+0.25*d+j*d - b2*(1+e2)/(1-e2)/np.sqrt(2),0),
                (0.5*delta*d+0.25*d+j*d,- b2*(1-e2)/(1+e2)/np.sqrt(2)),
            ]
        )

        hole = td.Structure(
            geometry = td.PolySlab(
                axis = 2,
                reference_plane = 'middle',
                slab_bounds = [0.5*dist,0.5*dist+h],
                vertices = vertices,
            ),
            medium = mat_envir,
            name = holename,
        )

        sim_structures.append(hole)

    # Upper line 
    for j in np.arange(0,Ncell):
        holename = 'holetru_'+str(j)

        vertices = np.array(
            [
                (0.5*delta*d+0.75*d+j*d + b2*(1+e2)/(1-e2)/np.sqrt(2),0.5*d),
                (0.5*delta*d+0.75*d+j*d,  b2*(1-e2)/(1+e2)/np.sqrt(2)+0.5*d),
                (0.5*delta*d+0.75*d+j*d - b2*(1+e2)/(1-e2)/np.sqrt(2),0.5*d),
                (0.5*delta*d+0.75*d+j*d,- b2*(1-e2)/(1+e2)/np.sqrt(2)+0.5*d),
            ]
        )

        hole = td.Structure(
            geometry = td.PolySlab(
                axis = 2,
                reference_plane = 'middle',
                slab_bounds = [0.5*dist,0.5*dist+h],
                vertices = vertices,
            ),
            medium = mat_envir,
            name = holename,
        )

        sim_structures.append(hole)

    # Lower line 
    for j in np.arange(0,Ncell):
        holename = 'holetrl_'+str(j)

        vertices = np.array(
            [
                (0.5*delta*d+0.75*d+j*d + b2*(1+e2)/(1-e2)/np.sqrt(2),-0.5*d),
                (0.5*delta*d+0.75*d+j*d,  b2*(1-e2)/(1+e2)/np.sqrt(2)-0.5*d),
                (0.5*delta*d+0.75*d+j*d - b2*(1+e2)/(1-e2)/np.sqrt(2),-0.5*d),
                (0.5*delta*d+0.75*d+j*d,- b2*(1-e2)/(1+e2)/np.sqrt(2)-0.5*d),
            ]
        )

        hole = td.Structure(
            geometry = td.PolySlab(
                axis = 2,
                reference_plane = 'middle',
                slab_bounds = [0.5*dist,0.5*dist+h],
                vertices = vertices,
            ),
            medium = mat_envir,
            name = holename,
        )

        sim_structures.append(hole)

    ### The left-hand side slab
    # Central line 
    for j in np.arange(0,Ncell):
        holename = 'holetlc_'+str(j)

        vertices = np.array(
            [
                (0.5*delta*d-0.75*d-j*d + b1*(1+e1)/(1-e1)/np.sqrt(2),0),
                (0.5*delta*d-0.75*d-j*d,  b1*(1-e1)/(1+e1)/np.sqrt(2)),
                (0.5*delta*d-0.75*d-j*d - b1*(1+e1)/(1-e1)/np.sqrt(2),0),
                (0.5*delta*d-0.75*d-j*d,- b1*(1-e1)/(1+e1)/np.sqrt(2)),
            ]
        )

        hole = td.Structure(
            geometry = td.PolySlab(
                axis = 2,
                reference_plane = 'middle',
                slab_bounds = [0.5*dist,0.5*dist+h],
                vertices = vertices,
            ),
            medium = mat_envir,
            name = holename,
        )

        sim_structures.append(hole)

    # Upper line 
    for j in np.arange(0,Ncell):
        holename = 'holetlu_'+str(j)

        vertices = np.array(
            [
                (0.5*delta*d-0.25*d-j*d + b1*(1+e1)/(1-e1)/np.sqrt(2),0.5*d),
                (0.5*delta*d-0.25*d-j*d,  b1*(1-e1)/(1+e1)/np.sqrt(2)+0.5*d),
                (0.5*delta*d-0.25*d-j*d - b1*(1+e1)/(1-e1)/np.sqrt(2),0.5*d),
                (0.5*delta*d-0.25*d-j*d,- b1*(1-e1)/(1+e1)/np.sqrt(2)+0.5*d),
            ]
        )

        hole = td.Structure(
            geometry = td.PolySlab(
                axis = 2,
                reference_plane = 'middle',
                slab_bounds = [0.5*dist,0.5*dist+h],
                vertices = vertices,
            ),
            medium = mat_envir,
            name = holename,
        )

        sim_structures.append(hole)

    # Lower line 
    for j in np.arange(0,Ncell):
        holename = 'holetll_'+str(j)

        vertices = np.array(
            [
                (0.5*delta*d-0.25*d-j*d + b1*(1+e1)/(1-e1)/np.sqrt(2),-0.5*d),
                (0.5*delta*d-0.25*d-j*d,  b1*(1-e1)/(1+e1)/np.sqrt(2)-0.5*d),
                (0.5*delta*d-0.25*d-j*d - b1*(1+e1)/(1-e1)/np.sqrt(2),-0.5*d),
                (0.5*delta*d-0.25*d-j*d,- b1*(1-e1)/(1+e1)/np.sqrt(2)-0.5*d),
            ]
        )

        hole = td.Structure(
            geometry = td.PolySlab(
                axis = 2,
                reference_plane = 'middle',
                slab_bounds = [0.5*dist,0.5*dist+h],
                vertices = vertices,
            ),
            medium = mat_envir,
            name = holename,
        )

        sim_structures.append(hole)

    ##### ==============================================================================
    ##### The bottom slab
    ### The right-hand side slab
    # Central line 
    for j in np.arange(0,Ncell):
        holename = 'holebrc_'+str(j)

        vertices = np.array(
            [
                (-0.5*delta*d+0.25*d+j*d + b1*(1+e1)/(1-e1)/np.sqrt(2),0),
                (-0.5*delta*d+0.25*d+j*d,  b1*(1-e1)/(1+e1)/np.sqrt(2)),
                (-0.5*delta*d+0.25*d+j*d - b1*(1+e1)/(1-e1)/np.sqrt(2),0),
                (-0.5*delta*d+0.25*d+j*d,- b1*(1-e1)/(1+e1)/np.sqrt(2)),
            ]
        )

        hole = td.Structure(
            geometry = td.PolySlab(
                axis = 2,
                reference_plane = 'middle',
                slab_bounds = [-0.5*dist-h,-0.5*dist],
                vertices = vertices,
            ),
            medium = mat_envir,
            name = holename,
        )

        sim_structures.append(hole)

    # Upper line 
    for j in np.arange(0,Ncell):
        holename = 'holebru_'+str(j)

        vertices = np.array(
            [
                (-0.5*delta*d+0.75*d+j*d + b1*(1+e1)/(1-e1)/np.sqrt(2),0.5*d),
                (-0.5*delta*d+0.75*d+j*d,  b1*(1-e1)/(1+e1)/np.sqrt(2)+0.5*d),
                (-0.5*delta*d+0.75*d+j*d - b1*(1+e1)/(1-e1)/np.sqrt(2),0.5*d),
                (-0.5*delta*d+0.75*d+j*d,- b1*(1-e1)/(1+e1)/np.sqrt(2)+0.5*d),
            ]
        )

        hole = td.Structure(
            geometry = td.PolySlab(
                axis = 2,
                reference_plane = 'middle',
                slab_bounds = [-0.5*dist-h,-0.5*dist],
                vertices = vertices,
            ),
            medium = mat_envir,
            name = holename,
        )

        sim_structures.append(hole)

    # Lower line 
    for j in np.arange(0,Ncell):
        holename = 'holebrl_'+str(j)

        vertices = np.array(
            [
                (-0.5*delta*d+0.75*d+j*d + b1*(1+e1)/(1-e1)/np.sqrt(2),-0.5*d),
                (-0.5*delta*d+0.75*d+j*d,  b1*(1-e1)/(1+e1)/np.sqrt(2)-0.5*d),
                (-0.5*delta*d+0.75*d+j*d - b1*(1+e1)/(1-e1)/np.sqrt(2),-0.5*d),
                (-0.5*delta*d+0.75*d+j*d,- b1*(1-e1)/(1+e1)/np.sqrt(2)-0.5*d),
            ]
        )

        hole = td.Structure(
            geometry = td.PolySlab(
                axis = 2,
                reference_plane = 'middle',
                slab_bounds = [-0.5*dist-h,-0.5*dist],
                vertices = vertices,
            ),
            medium = mat_envir,
            name = holename,
        )

        sim_structures.append(hole)

    ### The left-hand side slab
    # Central line 
    for j in np.arange(0,Ncell):
        holename = 'holeblc_'+str(j)

        vertices = np.array(
            [
                (-0.5*delta*d-0.75*d-j*d + b2*(1+e2)/(1-e2)/np.sqrt(2),0),
                (-0.5*delta*d-0.75*d-j*d,  b2*(1-e2)/(1+e2)/np.sqrt(2)),
                (-0.5*delta*d-0.75*d-j*d - b2*(1+e2)/(1-e2)/np.sqrt(2),0),
                (-0.5*delta*d-0.75*d-j*d,- b2*(1-e2)/(1+e2)/np.sqrt(2)),
            ]
        )

        hole = td.Structure(
            geometry = td.PolySlab(
                axis = 2,
                reference_plane = 'middle',
                slab_bounds = [-0.5*dist-h,-0.5*dist],
                vertices = vertices,
            ),
            medium = mat_envir,
            name = holename,
        )

        sim_structures.append(hole)

    # Upper line 
    for j in np.arange(0,Ncell):
        holename = 'holeblu_'+str(j)

        vertices = np.array(
            [
                (-0.5*delta*d-0.25*d-j*d + b2*(1+e2)/(1-e2)/np.sqrt(2),0.5*d),
                (-0.5*delta*d-0.25*d-j*d,  b2*(1-e2)/(1+e2)/np.sqrt(2)+0.5*d),
                (-0.5*delta*d-0.25*d-j*d - b2*(1+e2)/(1-e2)/np.sqrt(2),0.5*d),
                (-0.5*delta*d-0.25*d-j*d,- b2*(1-e2)/(1+e2)/np.sqrt(2)+0.5*d),
            ]
        )

        hole = td.Structure(
            geometry = td.PolySlab(
                axis = 2,
                reference_plane = 'middle',
                slab_bounds = [-0.5*dist-h,-0.5*dist],
                vertices = vertices,
            ),
            medium = mat_envir,
            name = holename,
        )

        sim_structures.append(hole)

    # Lower line 
    for j in np.arange(0,Ncell):
        holename = 'holebll_'+str(j)

        vertices = np.array(
            [
                (-0.5*delta*d-0.25*d-j*d + b2*(1+e2)/(1-e2)/np.sqrt(2),-0.5*d),
                (-0.5*delta*d-0.25*d-j*d,  b2*(1-e2)/(1+e2)/np.sqrt(2)-0.5*d),
                (-0.5*delta*d-0.25*d-j*d - b2*(1+e2)/(1-e2)/np.sqrt(2),-0.5*d),
                (-0.5*delta*d-0.25*d-j*d,- b2*(1-e2)/(1+e2)/np.sqrt(2)-0.5*d),
            ]
        )

        hole = td.Structure(
            geometry = td.PolySlab(
                axis = 2,
                reference_plane = 'middle',
                slab_bounds = [-0.5*dist-h,-0.5*dist],
                vertices = vertices,
            ),
            medium = mat_envir,
            name = holename,
        )

        sim_structures.append(hole)

    return sim_structures