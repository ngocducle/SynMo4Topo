import meep as mp
from meep import mpb 

### ================================================================================ ###
###                                                                                  ###
###     This module contains the functions to define the mode solver in              ### 
###     different structures, with specified materials                               ###
###                                                                                  ###
### ================================================================================ ###

##### FUNCTION: ModeSolver 2D photonic crystal slab monolayer with circular holes
def _2DSlab1LCircularHole(h,Lz,radius,num_bands,Nk,resolution,kSpace):
    ### Define the k-points
    if kSpace == 'BZ':  # The whole Brillouin zone
        k_points = [
            mp.Vector3(0.0,0.0,0.0),    # Gamma
            mp.Vector3(0.5,0.0,0.0),    # X
            mp.Vector3(0.5,0.5,0.0),    # M 
            mp.Vector3(0.0,0.0,0.0)     # Gamma 
        ]
    elif kSpace == 'M': # In the vicinity of the M-point
        k_points = [
            mp.Vector3(0.3,0.3,0.0),    # Gamma'
            mp.Vector3(0.5,0.3,0.0),    # X
            mp.Vector3(0.5,0.5,0.0),    # M 
            mp.Vector3(0.3,0.3,0.0)     # Gamma 
        ]
    else:
        print('ERROR! The k-point has not been in the allowed list yet')
        exit()

    # Interpolate to get the points on the lines
    k_points = mp.interpolate(Nk,k_points)

    ### Define the materials
    Si = mp.Medium(index = 3.54)
    SiO2 = mp.Medium(index = 1.46)
    PMMA = mp.Medium(index = 1.46)
    Dielectric = mp.Medium(epsilon = 12.0)
    Air = mp.Medium(epsilon = 1.0)

    Environment = PMMA 

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
            material = Environment
        ),

        mp.Block(
            center = mp.Vector3(0.0,0.0,0.0),
            size = mp.Vector3(1.0,1.0,h),
            material = Si
        ),

        mp.Cylinder(
            center = mp.Vector3(0.0,0.0,0.0),
            radius = radius,
            height = h,
            axis = mp.Vector3(0,0,1),
            material = Environment
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
                          num_bands,Nk,resolution,kSpace):
    ### Define the k-points
    if kSpace == 'BZ':  # The whole Brillouin zone
        k_points = [
            mp.Vector3(0.0,0.0,0.0),    # Gamma
            mp.Vector3(0.5,0.0,0.0),    # X
            mp.Vector3(0.5,0.5,0.0),    # M 
            mp.Vector3(0.0,0.0,0.0)     # Gamma 
        ]
    elif kSpace == 'M': # In the vicinity of the M-point
        k_points = [
            mp.Vector3(0.3,0.3,0.0),    # Gamma'
            mp.Vector3(0.5,0.3,0.0),    # X
            mp.Vector3(0.5,0.5,0.0),    # M 
            mp.Vector3(0.3,0.3,0.0)     # Gamma 
        ]
    else:
        print('ERROR! The k-point has not been in the allowed list yet')
        exit()

    # Interpolate to get the points on the lines
    k_points = mp.interpolate(Nk,k_points)

    ### Define the materials
    Si = mp.Medium(index = 3.54)
    SiO2 = mp.Medium(index = 1.46)
    Dielectric = mp.Medium(epsilon = 12.0)
    Air = mp.Medium(epsilon = 1.0)

    Environment = Air 

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
            material = Environment
            ),

        mp.Block(
            center = mp.Vector3(0,0,0.5*(h+dist)), 
            size = mp.Vector3(1,1,h), 
            material = Si
            ), 

        mp.Block(
            center = mp.Vector3(0,0,-0.5*(h+dist)),
            size = mp.Vector3(1,1,h),
            material = Si
            ), 

        mp.Cylinder(
            center = (0.5*delta1,0.5*delta2,0.5*(dist+h)), 
            radius = radius,  
            height = h, 
            axis = mp.Vector3(0,0,1), 
            material = Environment
            ),

        mp.Cylinder(
            center = (-0.5*delta1,-0.5*delta2,-0.5*(dist+h)), 
            radius = radius,  
            height = h, 
            axis = mp.Vector3(0,0,1), 
            material = Environment
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