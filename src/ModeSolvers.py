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
            mp.Vector3(0.45,0.45,0.0)   # Gamma'
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
            mp.Vector3(0.45,0.45,0.0),    # Gamma'
            mp.Vector3(0.5,0.45,0.0),     # X
            mp.Vector3(0.5,0.5,0.0),      # M 
            mp.Vector3(0.45,0.45,0.0)     # Gamma 
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

##### FUNCTION: ModeSolver 2D photonic crystal slab monolayer with square holes
def _2DSlab1LSquareHole(h,Lz,edge,num_bands,Nk,resolution,kSpace,Mater,Envir):
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
            mp.Vector3(0.45,0.45,0.0)   # Gamma'
        ]

        ToInterpolate = 'Yes'

    elif kSpace == 'M': # Consider the M-point only 
        k_points = [
            mp.Vector3(0.5,0.5,0.0)     # M
        ]

        ToInterpolate = 'No'

    else:
        #print('ERROR! The k-point has not been in the allowed list yet')
        #exit()
        k_points = kSpace 
        ToInterpolate = 'No'

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

        mp.Block(
            center = mp.Vector3(0.0,0.0,0.0),
            size = mp.Vector3(edge,edge,h),
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
        #tolerance = 1.e-7,
        #mesh_size = 10
    )

    ### Return the mode solver
    return ModeSolver

##### FUNCTION: ModeSolver 2D photonic crystal slab bilayer with square holes
def _2DSlab2LSquareHole(h,Lz,edge,dist,delta1,delta2,
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
            mp.Vector3(0.45,0.45,0.0),    # Gamma'
            mp.Vector3(0.5,0.45,0.0),     # X
            mp.Vector3(0.5,0.5,0.0),      # M 
            mp.Vector3(0.45,0.45,0.0)     # Gamma 
        ]

        ToInterpolate = 'Yes'

    elif kSpace == 'Gamma-M':   # Run along the line Gamma-M: kx = ky 
        k_points = [
            mp.Vector3(0.0,0.0,0.0),      # Gamma 
            mp.Vector3(0.5,0.5,0.0)       # M 
        ]

        ToInterpolate = 'Yes'

    elif kSpace == 'M': # Consider the M-point only 
        k_points = [
            mp.Vector3(0.5,0.5,0.0)     # M
        ]

        ToInterpolate = 'No'

    elif kSpace == 'Gamma-M-Mvicinity': # Run along the line Gamma-M: kx=ky, close to M
        k_points = [
            mp.Vector3(0.4,0.4,0.0),    # A point in the line Gamma-M
            mp.Vector3(0.5,0.5,0.0)     # M
        ]

        ToInterpolate = 'Yes'

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

        mp.Block(
            center = (0.5*delta1,0.5*delta2,0.5*(dist+h)),    
            size = mp.Vector3(edge,edge,h), 
            material = Envir
            ),

        mp.Block(
            center = (-0.5*delta1,-0.5*delta2,-0.5*(dist+h)),  
            size = mp.Vector3(edge,edge,h), 
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

##### FUNCTION: ModeSolver 2D photonic crystal slab monolayer with rhombus holes
def _2DSlab1LRhombusHole(h,Lz,b,e,vertices,
                         num_bands,k_points,resolution,Mater,Envir):
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
            size = mp.Vector3(1.0,1.0,mp.inf),
            material = Envir
        ),

        mp.Block(
            center = mp.Vector3(0.0,0.0,0.0),
            size = mp.Vector3(1.0,1.0,h),
            material = Mater 
        ),

        mp.Prism(
            vertices = vertices,
            height = h, 
            axis = mp.Vector3(0.0,0.0,1.0),
            center = mp.Vector3(0.0,0.0,0.0),
            material = Envir 
        )
    ]

    ### The mode solver 
    ModeSolver = mpb.ModeSolver(
        geometry = geometry,
        geometry_lattice = geometry_lattice,
        k_points = k_points,
        resolution = resolution,
        num_bands = num_bands 
    )

    ### Return the mode solver 
    return ModeSolver

##### FUNCTION: ModeSolver 2D photonic crystal slab bilayer with rhombus holes
def _2DSlab2LRhombusHole(h1,b1,e1,vertices1,
                         h2,b2,e2,vertices2,
                         Lz,dist,deltax,deltay,
                         num_bands,k_points,resolution,Mater,Envir):
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
            size = mp.Vector3(1.0,1.0,mp.inf),
            material = Envir 
        ),

        mp.Block(
            center = mp.Vector3(0.0,0.0,0.5*(h1+dist)),
            size = mp.Vector3(1.0,1.0,h1),
            material = Mater
        ),

        mp.Block(
            center = mp.Vector3(0.0,0.0,-0.5*(h2+dist)),
            size = mp.Vector3(1.0,1.0,h2),
            material = Mater 
        ),

        mp.Prism(
            vertices = vertices1,
            height = h1,
            axis = mp.Vector3(0.0,0.0,1.0),
            center = (deltax,deltay,0.5*(dist+h1)),
            material = Envir 
        ),

        mp.Prism(
            vertices = vertices2,
            height = h2,
            axis = mp.Vector3(0.0,0.0,1.0),
            center = (-deltax,-deltay,-0.5*(dist+h2)),
            material = Envir 
        )
    ]

    ### The mode solver 
    ModeSolver = mpb.ModeSolver(
        geometry = geometry,
        geometry_lattice = geometry_lattice,
        k_points = k_points,
        resolution = resolution,
        num_bands = num_bands
    )

    ### Return the mode solver 
    return ModeSolver 