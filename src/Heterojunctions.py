import numpy as np 
import scipy 
import meep as mp 

##### ======================================================================== ##### 
###                                                                              ###
###     This module contains the functions related to the FDTD simulations       ###
###                            of heterojunctions                                ###
###                                                                              ###
##### ======================================================================== #####

##### ==============================================================================
##### FUNCTION: define the geometry of heterojunction made of two 2DSlab2L-RHole
##### with the same lattice parameter 
### 
###   The flux propagates along the x direction 
###   d = diagonal of square unit cell 
###   h1,b1,e1: parameters of the left-hand side (the side of the incoming wave)
###   h2,b2,e2: parameters of the right-hand side (the side of the outgoing wave)
###   dist: interlayer distance (assumed to be the same for both sides) 
###   vertice_cell: the vertices of the unit cell (square rotated by 90 degrees)
###   vertice1: the vertices of the rhombus hole of the layer 1 
###   vertice2: the vertices of the rhombus hole of the layer 2 
###   Mater1: material for the left-hand side structure (both layers are made of the 
###           same material)
###   Mater2: material for the right-hand side structure (both layers are made of the 
###           same material)
###   Ncellx: The number of unit cells along the horizontal direction for each half
###   Ncelly: The number of unit cells along the vertical direction
###   sx,sy: The size of the simulation cell along the x and y directions, respectively 
###   structurex,structurey: The size of the structure along the x and y directions, 
###            respectively 
###   htotal: total height of the structure 
###   shift: the relative displacment between the two layers

def geo_2DSlab2L_RHole_hj(d,h1,b1,e1,h2,b2,e2,dist,
                          vertice_cell,vertice1,vertice2,
                          Mater1,Mater2,Envir,
                          Ncellx,Ncelly,sx,sy,
                          structurex,structurey,htotal,delta):
    ### Initialize the geometry with environment 
    geometry = [] 

    geometry.append(mp.Block(
        center = mp.Vector3(0,0,0),
        size = mp.Vector3(mp.inf,mp.inf,mp.inf),
        material = Envir
    )) 

    ### Add the structure (already includes the padding blocks)
    geometry.append(
        mp.Block(
            center = mp.Vector3(-0.25*sx,0,0.5*(htotal-h1)),
            size = mp.Vector3(0.5*sx,sy,h1),
            material = Mater1 
        )
    )

    geometry.append(
        mp.Block(
            center = mp.Vector3(-0.25*sx,0,0.5*(-htotal+h2)),
            size = mp.Vector3(0.5*sx,sy,h2),
            material = Mater2 
        )
    )

    geometry.append(
        mp.Block(
            center = mp.Vector3(0.25*sx,0,0.5*(htotal-h2)),
            size = mp.Vector3(0.5*sx,sy,h2),
            material = Mater2 
        )
    )

    geometry.append(
        mp.Block(
            center = mp.Vector3(0.25*sx,0,0.5*(-htotal+h1)),
            size = mp.Vector3(0.5*sx,sy,h1),
            material = Mater1 
        )
    )

    ##### Add the unit cells at the heterojunction boundary 
    for j in range(Ncelly):
        geometry.append(
            mp.Prism(
                vertices = vertice_cell,
                height = h1,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.25*d,0.5*structurey-0.5*d-j*d,0.5*(htotal-h1)),
                material = Mater1 
            )
        )

        geometry.append(
            mp.Prism(
                vertices = vertice_cell,
                height = h2,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.25*d,0.5*structurey-0.5*d-j*d,0.5*(-htotal+h2)),
                material = Mater2 
            )
        )

    for j in range(Ncelly+1):
        geometry.append(
            mp.Prism(
                vertices = vertice_cell,
                height = h2,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.25*d,0.5*structurey-j*d,0.5*(htotal-h2)),
                material = Mater2
            )
        )

        geometry.append(
            mp.Prism(
                vertices = vertice_cell,
                height = h1,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.25*d,0.5*structurey-j*d,0.5*(-htotal+h1)),
                material = Mater1
            )
        )

    ##### Add the holes 
    for j in range(Ncelly):
        for i in range(Ncellx):
            geometry.append(mp.Prism(
                vertices = vertice1,
                height=h1,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.5*structurex+0.5*d+i*d+0.5*delta,
                          0.5*structurey-0.5*d-j*d,
                          0.5*(htotal-h1)),
                material = Envir 
            ))

            geometry.append(mp.Prism(
                vertices = vertice2,
                height=h2,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.5*structurex+0.5*d+i*d-0.5*delta,
                          0.5*structurey-0.5*d-j*d,
                        0.5*(-htotal+h2)),
                material = Envir 
            ))

        
        for i in range(Ncellx):
            geometry.append(mp.Prism(
                vertices = vertice2,
                height=h2,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.75*d+i*d+0.5*delta,
                          0.5*structurey-0.5*d-j*d,
                        0.5*(htotal-h2)),
                material = Envir 
            ))

            geometry.append(mp.Prism(
                vertices = vertice1,
                height=h1,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.75*d+i*d-0.5*delta,
                          0.5*structurey-0.5*d-j*d,
                          0.5*(-htotal+h1)),
                material = Envir 
            ))

    for j in range(Ncelly+1):
        for i in range(Ncellx):
            geometry.append(mp.Prism(
                vertices = vertice1,
                height=h1,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.5*structurex+i*d+0.5*delta,
                        0.5*structurey-j*d,
                        0.5*(htotal-h1)),
                material = Envir 
            ))

            geometry.append(mp.Prism(
                vertices = vertice2,
                height=h2,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(-0.5*structurex+i*d-0.5*delta,
                        0.5*structurey-j*d,
                        0.5*(-htotal+h2)),
                material = Envir 
            ))

        for i in range(Ncellx+1):
            geometry.append(mp.Prism(
                vertices = vertice2,
                height=h2,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.25*d+i*d+0.5*delta,
                        0.5*structurey-j*d,
                        0.5*(htotal-h2)),
                material = Envir 
            ))

            geometry.append(mp.Prism(
                vertices = vertice1,
                height=h1,
                axis = mp.Vector3(0,0,1),
                center = mp.Vector3(0.25*d+i*d-0.5*delta,
                          0.5*structurey-j*d,
                          0.5*(-htotal+h1)),
                material = Envir 
            )) 

    return geometry 



##### ==============================================================================
##### FUNCTION: define the geometry of heterojunction made of two 2DSlab2L-RHole
##### with the same lattice parameter 
### 
###   The flux propagates along the x direction 
###   d = diagonal of square unit cell 
###   h1,b1,e1: parameters of the left-hand side (the side of the incoming wave)
###   h2,b2,e2: parameters of the right-hand side (the side of the outgoing wave)
###   dist: interlayer distance (assumed to be the same for both sides) 
###   vertice_cell: the vertices of the unit cell (square rotated by 90 degrees)
###   vertice1: the vertices of the rhombus hole of the layer 1 
###   vertice2: the vertices of the rhombus hole of the layer 2 
###   Mater1: material for the left-hand side structure (both layers are made of the 
###           same material)
###   Mater2: material for the right-hand side structure (both layers are made of the 
###           same material)
###   Ncell:  The number of unit cells along the horizontal direction for each half
###   sx,sy: The size of the simulation cell along the x and y directions, respectively 
###   structurex,structurey: The size of the structure along the x and y directions, 
###            respectively 
###   hbilayer: total height of the structure 
###   shift: the relative displacment between the two layers

def geo_2DSlab2L_RHole_hj_PBCy(d,h1,b1,e1,h2,b2,e2,dist,
                          vertice_cell,vertice1,vertice2,
                          Mater1,Mater2,Envir,
                          Ncell,sx,sy,
                          structurex,structurey,hbilayer,delta):
    
    ### -----------------------------------------------------------------------
    ### Initialize the geometry with environment     
    geometry = []

    geometry.append(mp.Block(
        center = mp.Vector3(0,0,0),
        size = mp.Vector3(mp.inf,mp.inf,mp.inf),
        material = Envir 
    ))

    ### -----------------------------------------------------------------------
    ### Add the structure (already includes the padding blocks)
    geometry.append(
        mp.Block(
            center = mp.Vector3(-0.25*sx,0,0.5*(hbilayer-h1)),
            size = mp.Vector3(0.5*sx,sy,h1),
            material = Mater1 
        )
    )

    geometry.append(
        mp.Block(
            center = mp.Vector3(-0.25*sx,0,0.5*(-hbilayer+h2)),
            size = mp.Vector3(0.5*sx,sy,h2),
            material = Mater2 
        )
    )

    geometry.append(
        mp.Block(
            center = mp.Vector3(0.25*sx,0,0.5*(hbilayer-h2)),
            size = mp.Vector3(0.5*sx,sy,h2),
            material = Mater2 
        )
    )

    geometry.append(
        mp.Block(
            center = mp.Vector3(0.25*sx,0,0.5*(-hbilayer+h1)),
            size = mp.Vector3(0.5*sx,sy,h1),
            material = Mater1 
        )
    )

    ##### --------------------------------------------------------------------------
    ##### Add the unit cells at the heterojunction boundary 
    geometry.append(
        mp.Prism(
            vertices = vertice_cell,
            height = h1,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(-0.25*d,0,0.5*(hbilayer-h1)),
            material = Mater1 
        )
    )

    geometry.append(
        mp.Prism(
            vertices = vertice_cell,
            height = h2,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(-0.25*d,0,0.5*(-hbilayer+h2)),
            material = Mater2
        )
    )

    geometry.append(
        mp.Prism(
            vertices = vertice_cell,
            height = h2,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(0.25*d,0,0.5*(hbilayer-h2)),
            material = Mater2 
        )
    )

    geometry.append(
        mp.Prism(
            vertices = vertice_cell,
            height = h1,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(0.25*d,0,0.5*(-hbilayer+h1)),
            material = Mater1
        )
    )

    ##### -------------------------------------------------------------------------------
    ##### Add the holes 
    for i in range(Ncell):
        geometry.append(mp.Prism(
            vertices = vertice1,
            height = h1,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(-0.5*structurex+0.5*d+i*d+0.5*delta,
                                0,0.5*(hbilayer-h1)),
            material = Envir
        ))


    for i in range(Ncell):
        geometry.append(mp.Prism(
            vertices = vertice2,
            height = h2,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(-0.5*structurex+0.5*d+i*d-0.5*delta,
                                0,0.5*(-hbilayer+h2)),
            material = Envir
        )) 

    for i in range(Ncell):
        geometry.append(mp.Prism(
            vertices = vertice2,
            height = h2,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(0.25*d+i*d+0.5*delta,
                                0,0.5*(hbilayer-h2)),
            material = Envir 
        ))

    for i in range(Ncell):
        geometry.append(mp.Prism(
            vertices = vertice1,
            height = h1,
            axis = mp.Vector3(0,0,1),
            center = mp.Vector3(0.25*d+i*d-0.5*delta,
                                0,0.5*(-hbilayer+h1)),
            material = Envir 
        ))

    return geometry 