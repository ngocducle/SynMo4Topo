import numpy as np
import scipy 
import meep as mp 

##### ============================================================================= #####
###                                                                                   ###
###         This module contains the functions related to the structures              ###
###                             of the FDTD simulations                               ###
###                                                                                   ###
##### ============================================================================= #####

##### ===================================================================================
##### FUNCTION: define the geometry of the 2D photonic crystal slab monolayer with 
##### square unit cell and rhombus hole. The structure is rotated so that 
##### the source emits the electromagnetic wave along the diagonal direction 
##### that we choose to be the x-direction 
###   
###   d: diagonal length of the unit cell 
###   h: thickness of the 2D photonic crystal slab
###   b: edge length of the undeformed square hole
###   e: deformation parameter (e = 0 for square hole, e !=0 for rhombus hole)
###   vertice: the set of vertices of the hole, assuming that its center is place at (0,0,0)
###   Mater: material of the slab
###   Envir: the environment
###   Ncell: the number of cells (Ncell = 0 if there is no structure)
###   sx: size of the simulation region along the x-direction
###   structurey: size of the structure along the y-direction

def geo_2DSlab1L_RHoleP(d,h,b,e,vertice,Mater,Envir,Ncell,sx,structurey):
    ### Initialize the geometry with environment 
    geometry = []

    geometry.append(mp.Block(
        center = mp.Vector3(0,0,0),
        size = mp.Vector3(mp.inf,mp.inf,mp.inf),
        material = Envir 
        )
    )

    ### Add the 2D slab of thickness h 
    geometry.append(mp.Block(
        center = mp.Vector3(0,0,0),
        size = mp.Vector3(1.5*sx,1.5*structurey,h),
        material = Mater 
        )
    )

    ### We divide 2 cases: 
    # If Ncell is odd: the central hole is located at the point (0,0,0)
    # If Ncell is even: there is no central block, and the holes start to 
    # be placed at points (0.5*d,0,0) and (-0.5*d,0,0)
    if Ncell%2 == 1:
        # The number of unit cells in each half
        Nhalf = int((Ncell-1)/2)
        print('Nhalf = '+str(Nhalf))

        # The cells on the central line 
        for j in np.arange(-Nhalf,Nhalf+1):
            geometry.append(
                mp.Prism(
                    vertices = vertice,
                    height = h,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3(j*d,0,0),
                    material = Envir 
                )
            )

        # The cells on the upper line and the lower line
        for j in np.arange(-Nhalf,Nhalf):
            # The cells on the upper line 
            geometry.append(
                mp.Prism(
                    vertices = vertice,
                    height = h,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3((j+0.5)*d,0.5*d,0),
                    material = Envir 
                )
            )

            # The cells on the lower line
            geometry.append(
                mp.Prism(
                    vertices = vertice,
                    height = h,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3((j+0.5)*d,-0.5*d,0),
                    material = Envir 
                )
            )


    else:
        # The number of unit cells in each half 
        Nhalf = int(Ncell/2)
        print('Nhalf = '+str(Nhalf))

        # The cells on the central line 
        for j in np.arange(-Nhalf,Nhalf):
            geometry.append(
                mp.Prism(
                    vertices = vertice,
                    height = h,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3((j+0.5)*d,0,0),
                    material = Envir 
                )
            )

        # The cells on the upper line and the lower line
        for j in np.arange(-Nhalf+1,Nhalf):
            # The upper line 
            geometry.append(
                mp.Prism(
                    vertices = vertice,
                    height = h,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3(j*d,0.5*d,0),
                    material = Envir 
                )
            )

            # The lower line 
            geometry.append(
                mp.Prism(
                    vertices = vertice,
                    height = h,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3(j*d,-0.5*d,0),
                    material = Envir 
                )
            )
        
    return geometry

##### ===================================================================================
##### FUNCTION: define the geometry of the heterojunction between two
##### 2D photonic crystal slab monolayers with 
##### square unit cell and rhombus hole. The structure is rotated so that 
##### the source emits the electromagnetic wave along the diagonal direction 
##### that we choose to be the x-direction 
##### ATTENTION! The two photonic crystal slabs are made of the SAME material
###   
###   d: diagonal length of the unit cell 
###   h: thickness of the 2D photonic crystal slab
###   b: edge length of the undeformed square hole
###   e: deformation parameter (e = 0 for square hole, e !=0 for rhombus hole)
###   vertice: the set of vertices of the hole, assuming that its center is place at (0,0,0)
###   Mater: material of the slab
###   Envir: the environment
###   Ncell: the number of cells (Ncell = 0 if there is no structure)
###   sx: size of the simulation region along the x-direction
###   structurey: size of the structure along the y-direction

def geo_2DSlab1L_RholeP_hj_sameMater(d,h,
                                     vertice1,vertice2,
                                     Mater,Envir,Ncell,sx,sy):
    ### Initialize the geometry with environment 
    geometry = []

    geometry.append(
        mp.Block(
            center = mp.Vector3(0,0,0),
            size = mp.Vector3(mp.inf,mp.inf,mp.inf),
            material = Envir 
        )
    )

    ### Add the 2D slab of thickness h 
    geometry.append(
        mp.Block(
            center = mp.Vector3(0,0,0),
            size = mp.Vector3(1.5*sx,1.5*sy,h),
            material = Mater
        )
    )

    ### The right-hand side slab
    for j in np.arange(0,Ncell):
        # Central line 
        geometry.append(
            mp.Prism(
                vertices = vertice2,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.25*d+j*d,0,0),
                material = Envir
            )
        )

        # Upper line 
        geometry.append(
            mp.Prism(
                vertices = vertice2,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.75*d+j*d,0.5*d,0),
                material = Envir
            )
        )

        # Lower line 
        geometry.append(
            mp.Prism(
                vertices = vertice2,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.75*d+j*d,-0.5*d,0),
                material = Envir
            )
        )

    ### The left-hand side slab
    for j in np.arange(0,Ncell):
        # The central line 
        geometry.append(
            mp.Prism(
                vertices = vertice1,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.75*d-j*d,0,0),
                material = Envir
            )
        )

        # The upper line 
        geometry.append(
            mp.Prism(
                vertices = vertice1,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.25*d-j*d,0.5*d,0),
                material = Envir
            )
        )

        # The lower line 
        geometry.append(
            mp.Prism(
                vertices = vertice1,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.25*d-j*d,-0.5*d,0),
                material = Envir
            )
        )

    return geometry

##### =======================================================================
##### FUNCTION: define the geometry of 2D photonic crystal slab bilayer 
##### square unit cell and rhombus hole. The structure is rotated so that
##### the source emits the electromagnetic wave along the diagonal direction
##### that we choose to be the x-direction
###      
###     d: diagonal length of the unit cell
###     h1: thickness of layer 1
###     h2: thickness of layer 2 
###     hbilayer = h1 + dist + h2: total thickness of bilayer
###     delta: shift 
###     vertice1: set of vertices of layer 1 
###     vertice2: set of vertices of layer 2 
###     Mater: material of the bilayer 
###     Envir: environment 
###     Ncell: number of unit cells along the x-direction (diagonal)
###     sx: cell structure along the x-direction
###     sy: cell structure along the y-direction 

def geo_2DSlab2L_RHoleP(d,h1,h2,hbilayer,delta,
                        vertice1,vertice2,
                        Mater,Envir,Ncell,sx,sy):
    ### Initialize the geometry with environment
    geometry = []

    geometry.append(mp.Block(
        center = mp.Vector3(0,0,0),
        size = mp.Vector3(mp.inf,mp.inf,mp.inf),
        material = Envir
        )
    )

    ### Add the 2D slabs of thickness h1 and h2 
    geometry.append(
        mp.Block(
            center = mp.Vector3(0,0,0.5*(hbilayer-h1)),
            size = mp.Vector3(1.5*sx,1.5*sy,h1),
            material = Mater 
        )
    )

    geometry.append(
        mp.Block(
            center = mp.Vector3(0,0,-0.5*(hbilayer-h2)),
            size = mp.Vector3(1.5*sx,1.5*sy,h2),
            material = Mater 
        )
    )

    ##### LAYER 1:
    ### We divide 2 cases: 
    # If Ncell is odd: the central hole is located at the point (0.5*delta,0,0.5*(hbilayer-h1))
    # If Ncell is even: there is no central block, and the holes start to 
    # be placed at points (0.5*d+0.5*delta,0,0.5*(hbilayer-h1)) and (-0.5*d+0.5*delta,0,0.5*(hbilayer-h1))
    if Ncell%2 == 1:
        # The number of unit cells in each half 
        Nhalf = int((Ncell-1)/2)
        print('Nhalf = '+str(Nhalf))

        # The cells on the central line 
        for j in np.arange(-Nhalf,Nhalf+1):
            geometry.append(
                mp.Prism(
                    vertices = vertice1,
                    height = h1,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3(j*d-0.5*delta,0,0.5*(hbilayer-h1)),
                    material = Envir 
                )
            )

        # The cells on the upper line and the lower line
        for j in np.arange(-Nhalf,Nhalf):
            # The cells on the upper line 
            geometry.append(
                mp.Prism(
                    vertices = vertice1,
                    height = h1,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3((j+0.5)*d-0.5*delta,0.5*d,0.5*(hbilayer-h1)),
                    material = Envir 
                )
            )

            # The cells on the lower line
            geometry.append(
                mp.Prism(
                    vertices = vertice1,
                    height = h1,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3((j+0.5)*d-0.5*delta,-0.5*d,0.5*(hbilayer-h1)),
                    material = Envir 
                )
            )

    else:
        # The number of unit cells in each half
        Nhalf = int(Ncell/2)
        print('Nhalf = '+str(Nhalf))

        # The cells on the central line 
        for j in np.arange(-Nhalf,Nhalf):
            geometry.append(
                mp.Prism(
                    vertices = vertice1,
                    height = h1,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3((j+0.5)*d-0.5*delta,0,0.5*(hbilayer-h1)),
                    material = Envir 
                )
            )

        # The cells on the upper line and the lower line
        for j in np.arange(-Nhalf+1,Nhalf):
            # The upper line 
            geometry.append(
                mp.Prism(
                    vertices = vertice1,
                    height = h1,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3(j*d-0.5*delta,0.5*d,0.5*(hbilayer-h1)),
                    material = Envir 
                )
            )

            # The lower line 
            geometry.append(
                mp.Prism(
                    vertices = vertice1,
                    height = h1,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3(j*d-0.5*delta,-0.5*d,0.5*(hbilayer-h1)),
                    material = Envir 
                )
            )

    ##### LAYER 2:
    ### We divide 2 cases: 
    # If Ncell is odd: the central hole is located at the point (-0.5*delta,0,-0.5*(hbilayer-h2))
    # If Ncell is even: there is no central block, and the holes start to 
    # be placed at points (0.5*d-0.5*delta,0,-0.5*(hbilayer-h2)) and (-0.5*d-0.5*delta,0,-0.5*(hbilayer-h2))
    if Ncell%2 == 1:
        # The number of unit cells in each half
        Nhalf = int((Ncell-1)/2)
        print('Nhalf = '+str(Nhalf))

        # The cells on the central line 
        for j in np.arange(-Nhalf,Nhalf+1):
            geometry.append(
                mp.Prism(
                    vertices = vertice2,
                    height = h2,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3(j*d+0.5*delta,0,-0.5*(hbilayer-h2)),
                    material = Envir 
                )
            )

        # The cells on the upper line and the lower line
        for j in np.arange(-Nhalf,Nhalf):
            # The cells on the upper line 
            geometry.append(
                mp.Prism(
                    vertices = vertice2,
                    height = h2,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3((j+0.5)*d+0.5*delta,0.5*d,-0.5*(hbilayer-h2)),
                    material = Envir 
                )
            )

            # The cells on the lower line
            geometry.append(
                mp.Prism(
                    vertices = vertice2,
                    height = h2,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3((j+0.5)*d+0.5*delta,-0.5*d,-0.5*(hbilayer-h2)),
                    material = Envir 
                )
            )

    else:
        # The number of unit cells in each half
        Nhalf = int(Ncell/2)
        print('Nhalf = '+str(Nhalf))

        # The cells on the central line 
        for j in np.arange(-Nhalf,Nhalf):
            geometry.append(
                mp.Prism(
                    vertices = vertice2,
                    height = h2,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3((j+0.5)*d+0.5*delta,0,-0.5*(hbilayer-h2)),
                    material = Envir 
                )
            )

        # The cells on the upper line and the lower line
        for j in np.arange(-Nhalf+1,Nhalf):
            # The upper line 
            geometry.append(
                mp.Prism(
                    vertices = vertice2,
                    height = h2,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3(j*d+0.5*delta,0.5*d,-0.5*(hbilayer-h2)),
                    material = Envir 
                )
            )

            # The lower line 
            geometry.append(
                mp.Prism(
                    vertices = vertice2,
                    height = h2,
                    axis = mp.Vector3(0,0,1),
                    center = mp.Vector3(j*d+0.5*delta,-0.5*d,-0.5*(hbilayer-h2)),
                    material = Envir 
                )
            )

    return geometry

##### =======================================================================
##### FUNCTION: define the geometry of 2D photonic crystal slab bilayer 
##### square unit cell and rhombus hole. The structure is rotated so that
##### the source emits the electromagnetic wave along the diagonal direction
##### that we choose to be the x-direction
###      
###     d: diagonal length of the unit cell
###     h1: thickness of layer 1
###     h2: thickness of layer 2 
###     hbilayer = h1 + dist + h2: total thickness of bilayer
###     delta: shift 
###     vertice1: set of vertices of layer 1 
###     vertice2: set of vertices of layer 2 
###     Mater: material of the bilayer 
###     Envir: environment 
###     Ncell: number of unit cells along the x-direction (diagonal)
###     sx: cell structure along the x-direction
###     sy: cell structure along the y-direction
###
### The lower layer is shifted by an amount delta with respect to the upper layer

def geo_2DSlab2L_RHoleP_hj_sameMater(d,h,hbilayer,delta,
                                     vertice1,vertice2,
                                     Mater,Envir,Ncell,sx,sy):
    ### Initialize the geometry with environment 
    geometry = []

    geometry.append(
        mp.Block(
            center = mp.Vector3(0,0,0),
            size = mp.Vector3(mp.inf,mp.inf,mp.inf),
            material = Envir 
        )
    )

    ### Add the 2D slabs of thickness 2*h+dist 
    """geometry.append(
        mp.Block(
            center = mp.Vector3(0,0,0.5*(hbilayer-h)),
            size = mp.Vector3(1.5*sx,1.5*sy,h),
            material = Mater 
        )
    )

    geometry.append(
        mp.Block(
            center = mp.Vector3(0,0,-0.5*(hbilayer-h)),
            size = mp.Vector3(1.5*sx,1.5*sy,h),
            material = Mater 
        )
    )"""

    geometry.append(
        mp.Block(
            center = mp.Vector3(0,0,0),
            size = mp.Vector3(mp.inf,mp.inf,hbilayer),
            material = Mater
        )
    )

    dist = hbilayer - 2*h

    geometry.append(
        mp.Block(
            center = mp.Vector3(0,0,0),
            size = mp.Vector3((2*Ncell+1)*d,mp.inf,dist),
            material = Envir 
        )
    )

    ##### LAYER 1:
    ### The right-hand side slab
    for j in np.arange(0,Ncell):
        # Central line 
        geometry.append(
            mp.Prism(
                vertices = vertice2,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.25*d+j*d-0.5*delta,0,0.5*(hbilayer-h)),
                material = Envir
            )
        )

        # Upper line 
        geometry.append(
            mp.Prism(
                vertices = vertice2,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.75*d+j*d-0.5*delta,0.5*d,0.5*(hbilayer-h)),
                material = Envir
            )
        )

        # Lower line 
        geometry.append(
            mp.Prism(
                vertices = vertice2,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.75*d+j*d-0.5*delta,-0.5*d,0.5*(hbilayer-h)),
                material = Envir
            )
        )

    ### The left-hand side slab
    for j in np.arange(0,Ncell):
        # The central line 
        geometry.append(
            mp.Prism(
                vertices = vertice1,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.75*d-j*d-0.5*delta,0,0.5*(hbilayer-h)),
                material = Envir
            )
        )

        # The upper line 
        geometry.append(
            mp.Prism(
                vertices = vertice1,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.25*d-j*d-0.5*delta,0.5*d,0.5*(hbilayer-h)),
                material = Envir
            )
        )

        # The lower line 
        geometry.append(
            mp.Prism(
                vertices = vertice1,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.25*d-j*d-0.5*delta,-0.5*d,0.5*(hbilayer-h)),
                material = Envir
            )
        )

    ##### LAYER 2:
    ### The right-hand side slab
    for j in np.arange(0,Ncell):
        # Central line 
        geometry.append(
            mp.Prism(
                vertices = vertice1,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.25*d+j*d+0.5*delta,0,-0.5*(hbilayer-h)),
                material = Envir
            )
        )

         # Upper line 
        geometry.append(
            mp.Prism(
                vertices = vertice1,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.75*d+j*d+0.5*delta,0.5*d,-0.5*(hbilayer-h)),
                material = Envir
            )
        )

        # Lower line 
        geometry.append(
            mp.Prism(
                vertices = vertice1,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.75*d+j*d+0.5*delta,-0.5*d,-0.5*(hbilayer-h)),
                material = Envir
            )
        )

    ### The left-hand side slab
    for j in np.arange(0,Ncell):
        # The central line 
        geometry.append(
            mp.Prism(
                vertices = vertice2,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.75*d-j*d+0.5*delta,0,-0.5*(hbilayer-h)),
                material = Envir
            )
        )

        # The upper line 
        geometry.append(
            mp.Prism(
                vertices = vertice2,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.25*d-j*d+0.5*delta,0.5*d,-0.5*(hbilayer-h)),
                material = Envir
            )
        )

        # The lower line 
        geometry.append(
            mp.Prism(
                vertices = vertice2,
                height = h,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.25*d-j*d+0.5*delta,-0.5*d,-0.5*(hbilayer-h)),
                material = Envir
            )
        )

    return geometry 