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
###   

def geo_2DSlab1L_RHole(d,h,b,e,vertice_cell,
                       Mater,Envir,
                       Ncell,sx,sy,Lz,
                       structurex,structurey):
    
    ### Initialize the geometry with environment 
    geometry = []


    return geometry 