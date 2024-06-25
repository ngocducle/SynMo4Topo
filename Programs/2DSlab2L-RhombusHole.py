import numpy as np 
import meep as mp 
from meep import mpb 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab2LRhombusHole
from DielectricProfile import * 
from BandStructure import * 
from ExportData import *  
from Materials import * 
from LightCone import LightCone 

#################################################################################
#                                                                               #
#   This code studies a 2DSlab2L rhombus holes with:                            #
#   Slab thickness h, the parent square holes has edge length b                 #
#   and the deformation parameter of the rhombus is e                           #
#                                                                               #
#   The two slabs are displaced by distances delta1 and delta2                  #
#   along the x and y directions, respectively                                  #
#                                                                               #
#################################################################################

##### The main program goes here
def main():
    ###################################################################################
    #                                                                                 #
    #           Please define the parameters of the simulation here                   #
    #                                                                                 #
    ###################################################################################

    ### The light polarization 
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'all' # For 2D multilayer cases, always choose 'all'
    print('# Polarization = '+polarization)

    ### Resolution 
    resolution = mp.Vector3(32,32,32)  # pixels/a 
    print('# The resolution:'+str(resolution))

    ### Number of bands
    num_bands = 20 
    print('# The number of bands to simulate: '+str(num_bands))

    ### Geometrical parameters 
    # The upper layer 
    h1 = 0.36   # Thickness of the upper layer
    b1 = 0.3    # The edge length of the undeformed square hole 
    #d1 = b1/np.sqrt(2) # The projection of half the mean diagonal on the x and y axes
    e1 = 0.05   # The anisotropy between the two diagonals of the upper layer
    
    print('\n# Upper slab:')
    print('# Thickness h1 = '+str(h1))
    print('# The edge length of the undeformed square hole b1 = '+str(b1))
    print('# The anisotropy between the two diagonals e1 = '+str(e1))

    # The lower layer 
    h2 = 0.35   # Thickness of the lower layer
    b2 = 0.40   # The edge length of the undeformed square hole    
    #d2 = b2/np.sqrt(2) # The projection of the mean diagonal on the x and y axes
    e2 = 0.05   # The anisotropy between the two diagonals of the lower layer

    print('\n# Lower slab:')
    print('# Thickness h2 = '+str(h2))
    print('# The edge length of the undeformed square hole b2 = '+str(b2))
    print('# The anisotropy between the two diagonals e2 = '+str(e2))

    # The interlayer parameters 
    dist = 0.05 # Distance between the two layers 
    deltax = 0.0 # Relative displacement along the x-direction 
    deltay = 0.0 # Relative displacement along the y-direction 
    Lz = 10.0     # The height of the unit cell along the z-direction 

    print('\n# The distance between the two layers d = '+str(dist))
    print('# Relative displacement along x-direction deltax = '+str(deltax))
    print('# Relative displacement along y-direction deltay = '+str(deltay))
    print('# The height of the simulation cell Lz = '+str(Lz))

    ### The vertices of the rhombus holes
    vertices1 = [
        mp.Vector3( 0.5*(1+e1)*b1+0.5*deltax,  0.5*(1+e1)*b1+0.5*deltay, 0.0 ),
        mp.Vector3( 0.5*(1-e1)*b1+0.5*deltax, -0.5*(1-e1)*b1+0.5*deltay, 0.0 ),
        mp.Vector3(-0.5*(1+e1)*b1+0.5*deltax, -0.5*(1+e1)*b1+0.5*deltay, 0.0 ),
        mp.Vector3(-0.5*(1-e1)*b1+0.5*deltax,  0.5*(1-e1)*b1+0.5*deltay, 0.0 ) 
    ]

    vertices2 = [
        mp.Vector3( 0.5*(1+e2)*b2-0.5*deltax,  0.5*(1+e2)*b2-0.5*deltay, 0.0 ),
        mp.Vector3( 0.5*(1-e2)*b2-0.5*deltax, -0.5*(1-e2)*b2-0.5*deltay, 0.0 ),
        mp.Vector3(-0.5*(1+e2)*b2-0.5*deltax, -0.5*(1+e2)*b2-0.5*deltay, 0.0 ),
        mp.Vector3(-0.5*(1-e2)*b2-0.5*deltax,  0.5*(1-e2)*b2-0.5*deltay, 0.0 ) 
    ]

    ### The k-point at which we plot the field profile
    k_field = mp.Vector3(0.5,0.5,0.0)   # M-point
    print('# The k-point at which we plot the field profile:'+str(k_field))

    ### Show figure (Yes/No)
    show_fig = 'Yes'
    print('\n# Show the figure: '+str(show_fig))

    ### Materials
    # The material forming the structure 
    Mater = Si 

    # The environment
    Envir = PMMA 

    ### The part of the momentum space to plot the band structure 
    # Choose between: 'BZ', 'Gamma-M', 'M', 'M-vicinity'
    kSpace = 'BZ'
    print('# The part of the momentum space to simulate: '+kSpace) 

    ### Number of k-points to interpolate between the 2 high-symmetry points
    Nk = 2
    print('# The number of points to interpolate the high-symmetry line Nk = '+str(Nk))

    # The set of k-points 
    # Due to the mirror symmetry with respect to the planes 
    # x = y and x = -y, the lines Gamma-X and Gamma-Y have 
    # the same dispersion, and the lines MX and MY have the same 
    # dispersion
    k_points = [
        mp.Vector3(0.0,0.0,0.0),    # Gamma 
        mp.Vector3(0.5,0.0,0.0),    # X 
        mp.Vector3(0.5,0.5,0.0),    # M+ 
        mp.Vector3(0.0,0.0,0.0),    # Gamma  
        mp.Vector3(-0.5,0.5,0.0)    # M- 
    ]

    k_points = mp.interpolate(Nk,k_points)

    ### Show the figure
    show_fig = 'Yes'

    ### The title and the name of the files
    namesave = '2DSlab2L-RhombusHole-h1_'+str(h1)+'-b1_'+str(b1)+'-e1_'+str(e1) \
        + '-h2_'+str(h2)+'-b2_'+str(b2)+'-e2_'+str(e2) \
        + '-d_'+str(dist)+'-deltax_'+str(deltax)+'-deltay_'+str(deltay) \
        + '-'+polarization

    ### Calculate the dispersion for the light cone with Envir and k_points
    lightcone = LightCone(Envir,k_points,resolution)

    #################################################################################
    #                                                                               #
    #                   Here the simulation starts                                  #
    #                                                                               #
    #################################################################################

    ### Define the mode solver 
    ms = _2DSlab2LRhombusHole(h1,b1,e1,vertices1,
                              h2,b2,e2,vertices2,
                              Lz,dist,deltax,deltay,
                              num_bands,k_points,resolution,Mater,Envir)
    
    ### Run the simulation 
    if polarization == 'all':
        ms.run()
    elif polarization == 'zeven':
        ms.run_zeven()
    elif polarization == 'zodd':
        ms.run_zodd()
    else:
        print('ERROR! The polarization does not belong to the allowed list')
        exit()

    ### Extract the frequencies of the modes from the ModeSolver
    freqs = ms.all_freqs 

    ### The number of elemenets in k_points
    number = np.arange(len(ms.k_points))

    ### Output the dielectric profile along the z-direction with projected (X,Y)
    x = 0.50 
    y = 0.50 
    zmin = -0.5*Lz 
    zmax = 0.5*Lz 
    Nz = 50 

    # Calculate the dielectric profile along the z-direction at fixed (x,y)
    z_array,epsilon_z_array = DielectricProfileZ(ms,x,y,zmin,zmax,Nz)

    # Plot the dielectric profile, the name of the figure is:
    #           namesave+'-epsilon-z.png'
    PlotDielectricProfileZ(x,y,z_array,epsilon_z_array,namesave,show_fig)

    # Print the dielectric profile to the file:
    #           namesave+'-epsilon-z.txt'
    PrintDielectricProfileZ(x,y,z_array,epsilon_z_array,namesave)



    ### Output the dielectric profile with Ncellx x Ncelly unit cells
    ### for zmin <= z <= zmax, with Nz values of z 
    Ncellx = 5
    Ncelly = 5
    zmin = -0.5*dist-h1 
    zmax = 0.5*dist+h2 
    Nx = 300 
    Ny = 300 
    Nz = 13 

    # Calculate the dielectric profile in planes parallel to Oxy 
    x_plot,y_plot,z_array,epsilon_xy_array \
        = DielectricProfileXY(ms,Ncellx,Ncelly,zmin,zmax,Nx,Ny,Nz)
    
    # Plot the dielectric profile, the name of the figure is:
    #           namesave+'-z_'+str(k)+'.png'
    # where k is the number of the value of z in the array z_array
    PlotDielectricProfileXY(x_plot,y_plot,z_array,epsilon_xy_array,namesave,show_fig)
    
    # Print the dielectric profile to the file:
    #           namesave+'-epsilon-xy.txt'
    PrintDielectricProfileXY(x_plot,y_plot,z_array,epsilon_xy_array,namesave)

    


    ### Print the band structure to file
    PrintBandStructure(freqs,number,namesave)

    ### Plot the band structure
    if kSpace == 'BZ':
        PlotBand_BrillouinZone_Scell_Rhole(number,freqs,Nk,lightcone,namesave,show_fig)
    elif kSpace == 'M-vicinity':
        PlotBand_M(number,freqs,Nk,namesave,show_fig)
    else:
        print('ERROR! The k-point has not been in the allowed list yet')
        exit()

##### Run the MAIN program 
if __name__ == "__main__":
    main()