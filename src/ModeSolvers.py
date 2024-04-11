import meep as mp
from meep import mpb 

### ================================================================================ ###
###                                                                                  ###
###     This module contains the functions to define the mode solver in              ### 
###     different structures, with specified materials                               ###
###                                                                                  ###
### ================================================================================ ###

##### FUNCTION: ModeSolver 2D photonic crystal slab monolayer with circular holes
def _2DSlab1LCircularHole(h,Lz,radius,num_bands,Nk,resolution,kSpace,Mater,Envir):
    ### Define the k-points
    if kSpace == 'BZ':  # The whole Brillouin zone
        k_points = [
            mp.Vector3(0.0,0.0,0.0),    # Gamma
            mp.Vector3(0.5,0.0,0.0),    # X
            mp.Vector3(0.5,0.5,0.0),    # M 
            mp.Vector3(0.0,0.0,0.0)     # Gamma 
        ]

        ToInterpolate = 'Yes'

    elif kSpace == 'M-vicinity': # In the vicinity of the M-point
        k_points = [
            mp.Vector3(0.45,0.45,0.0),  # Gamma'
            mp.Vector3(0.5,0.45,0.0),   # X'
            mp.Vector3(0.5,0.5,0.0),    # M 
            mp.Vector3(0.45,0.45,0.0)   # Gamma 
        ]

        ToInterpolate = 'Yes'

    elif kSpace == 'M': # Consider the M-point only 
        k_points = [
            mp.Vector3(0.5,0.5,0.0)     # M
        ]

        ToInterpolate = 'No'

    else:
        print('ERROR! The k-point has not been in the allowed list yet')
        exit()

    # Interpolate to get the points on the lines
    if ToInterpolate == 'Yes':
        k_points = mp.interpolate(Nk,k_points)

    ### Define the lattice
    geometry_lattice = mp.Lattice(
        size = mp.Vector3(1.0,1.0,Lz),
        basis1 = mp.Vector3(1.0,0.0),
        basis2 = mp.Vector3(0.0,1.0)
    )

    ### Define the geometry
    geometry = [
        mp.Block(
            center = mp.Vector3(0.0,0.0,0.0),
            size = mp.Vector3(mp.inf,mp.inf,mp.inf),
            material = Envir 
        ),

        mp.Block(
            center = mp.Vector3(0.0,0.0,0.0),
            size = mp.Vector3(1.0,1.0,h),
            material = Mater 
        ),

        mp.Cylinder(
            center = mp.Vector3(0.0,0.0,0.0),
            radius = radius,
            height = h,
            axis = mp.Vector3(0,0,1),
            material = Envir 
        )
    ]

    ### The ModeSolver
    ModeSolver = mpb.ModeSolver(
        geometry = geometry,
        geometry_lattice = geometry_lattice,
        k_points = k_points,
        resolution = resolution,
        num_bands = num_bands
    )

    ### Return the mode solver
    return ModeSolver

##### FUNCTION: ModeSolver 2D photonic crystal slab bilayer with circular holes
def _2DSlab2LCircularHole(h,Lz,radius,dist,delta1,delta2,
                          num_bands,Nk,resolution,kSpace,Mater,Envir):
    ### Define the k-points
    if kSpace == 'BZ':  # The whole Brillouin zone
        k_points = [
            mp.Vector3(0.0,0.0,0.0),    # Gamma
            mp.Vector3(0.5,0.0,0.0),    # X
            mp.Vector3(0.5,0.5,0.0),    # M 
            mp.Vector3(0.0,0.0,0.0)     # Gamma 
        ]

        ToInterpolate = 'Yes'

    elif kSpace == 'M-vicinity': # In the vicinity of the M-point
        k_points = [
            mp.Vector3(0.3,0.3,0.0),    # Gamma'
            mp.Vector3(0.5,0.3,0.0),    # X
            mp.Vector3(0.5,0.5,0.0),    # M 
            mp.Vector3(0.3,0.3,0.0)     # Gamma 
        ]

        ToInterpolate = 'Yes'

    elif kSpace == 'M': # Consider the M-point only 
        k_points = [
            mp.Vector3(0.5,0.5,0.0)     # M
        ]

        ToInterpolate = 'No'

    else:
        print('ERROR! The k-point has not been in the allowed list yet')
        exit()

    # Interpolate to get the points on the lines
    if ToInterpolate == 'Yes':
        k_points = mp.interpolate(Nk,k_points)

    ### Define the lattice
    geometry_lattice = mp.Lattice(
        size = mp.Vector3(1.0,1.0,Lz),
        basis1 = mp.Vector3(1.0,0.0),
        basis2 = mp.Vector3(0.0,1.0)
    )

    # Define the geometry 
    geometry = [ 
        mp.Block(
            center = mp.Vector3(0,0,0),  
            size = mp.Vector3(mp.inf, mp.inf, mp.inf), 
            material = Envir
            ),

        mp.Block(
            center = mp.Vector3(0,0,0.5*(h+dist)), 
            size = mp.Vector3(1,1,h), 
            material = Mater 
            ), 

        mp.Block(
            center = mp.Vector3(0,0,-0.5*(h+dist)),
            size = mp.Vector3(1,1,h),
            material = Mater 
            ), 

        mp.Cylinder(
            center = (0.5*delta1,0.5*delta2,0.5*(dist+h)), 
            radius = radius,  
            height = h, 
            axis = mp.Vector3(0,0,1), 
            material = Envir
            ),

        mp.Cylinder(
            center = (-0.5*delta1,-0.5*delta2,-0.5*(dist+h)), 
            radius = radius,  
            height = h, 
            axis = mp.Vector3(0,0,1), 
            material = Envir
            )    
    ] 

    ### The ModeSolver
    ModeSolver = mpb.ModeSolver(
        geometry = geometry,
        geometry_lattice = geometry_lattice,
        k_points = k_points,
        resolution = resolution,
        num_bands = num_bands
    )

    ### Return the mode solver
    return ModeSolver