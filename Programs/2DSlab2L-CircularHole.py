import numpy as np
import meep as mp 
from meep import mpb

import sys
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab2LCircularHole
from DielectricProfile import *
from BandStructure import *
from ExportData import *
from EFieldProfile import *
from HFieldProfile import *


##### The MAIN program goes here 
def main():
    ##########################################################################
    #                                                                        #
    #        Please define the parameters of the simulation here             #
    #                                                                        #
    ##########################################################################


    ### The light polarization
    # Choose between: 'all', 'zeven', 'zodd' 
    polarization = 'zeven'
    print('# Polarization = '+polarization)

    ### The part of the momentum space to plot the band structure
    # Choose between: 'BZ', 'M' 
    kSpace = 'BZ'
    print('# The part of the momentum space to simulate:'+kSpace)

    ### The k-point at which we plot the field profile
    k_field = mp.Vector3(0.5,0.5,0.0)   # M-point
    print('# The k-point at which we plot the field profile:'+str(k_field))

    ### Resolution 
    resolution = mp.Vector3(32,32,32)   # pixels/a
    print('# The resolution:'+str(resolution))

    ### Geometrical parameters
    h = 0.3      # Thickness of the slab
    radius = 0.4 # Radius of the hole
    dist = 2.0   # Distance between the two layers 
    delta1 = 0.0 # Relative displacement along the x-direction
    delta2 = 0.0 # Relative displacement along the y-direction
    Lz = 10.0    # The height of the unit cell along the z-direction

    print('# Thickness of the slab h = '+str(h))
    print('# Radius of the circular hole r = '+str(radius))
    print('# The height of the simulation cell Lz = '+str(Lz))
    print('# Relative shift along x-direction delta1 = '+str(delta1))
    print('# Relative shift along y-direction delta2 = '+str(delta2))
    print('# Interlayer distance dist = '+str(dist))

    ### Number of bands
    num_bands = 20
    print('# The number of bands to simulate: '+str(num_bands))

    ### Number of k-points to interpolate between 2 high-symmetry points
    Nk = 19  
    print('# The number of points to interpolate the high-symmetry line Nk = '+str(Nk))

    ### Show figure (Yes/No)
    show_fig = 'No'
    print('# Show the figure: '+str(show_fig))

    ### The parameters to calculate the E- and H-field profiles
    ### in the planes parallel to Oxy
    Bloch_Phase = 'True'   # Enable the Bloch phase: True/False
    resolution_eps = 81    # Number of pixels per a 
    resolution_field = 81  # Number of pixels per a 
    num_periods = 3        # Number of periods along each direction  

    # The value of z where we take a slice to plot the fields 
    zvalue = 1.225

    ##########################################################################
    #                                                                        #
    #                   Here the simulation starts                           #
    #                                                                        #
    ##########################################################################

    ### Define the mode solver
    ms = _2DSlab2LCircularHole(h,Lz,radius,dist,delta1,delta2,
                          num_bands,Nk,resolution,kSpace)

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

    ### The number of elements in k_points
    number = np.arange(len(ms.k_points))

    ### The title and the name of the files 
    namesave = '2DSlab2L-CircularHole-h_'+str(h)+'-r_'+str(radius) \
        +'-'+polarization+'-'+kSpace 
    


    ### Output the dielectric profile along the z-direction with projected (X,Y)
    x = 0.27 
    y = 0.36 
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
    zmin = -0.5*Lz
    zmax = 0.5*Lz 
    Nx = 300
    Ny = 300
    Nz = 9

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
        PlotBand_BrillouinZone(number,freqs,Nk,namesave,show_fig)
    elif kSpace == 'M':
        PlotBand_M(number,freqs,Nk,namesave,show_fig)
    else:
        print('ERROR! The k-point has not been in the allowed list yet')
        exit()

    # Define the mode solver to calculate the E-field
    Efieldx,Efieldy,Efieldz,X,Y,Xfield,Yfield,eps_Oxy = EFields_2DSlab2LCircularHole(h,Lz,
                                radius,dist,delta1,delta2,num_bands,resolution,
                                k_field,zvalue,
                                polarization,resolution_eps,resolution_field,
                                num_periods,Bloch_Phase)
    
    # Plot the E-field
    Plot_ReEfield_Profile(Efieldx,Efieldy,Efieldz,zvalue,
                       X,Y,eps_Oxy,Xfield,Yfield,num_periods,show_fig)
    
    Plot_ImEfield_Profile(Efieldx,Efieldy,Efieldz,zvalue,
                       X,Y,eps_Oxy,Xfield,Yfield,num_periods,show_fig)
    
    # Define the mode solver to calculate the H-field
    Hfieldx,Hfieldy,Hfieldz,X,Y,Xfield,Yfield,eps_Oxy = HFields_2DSlab2LCircularHole(h,Lz,
                                radius,dist,delta1,delta2,num_bands,resolution,
                                k_field,zvalue,
                                polarization,resolution_eps,resolution_field,
                                num_periods,Bloch_Phase)

    # Plot the H-field
    Plot_ReHfield_Profile(Hfieldx,Hfieldy,Hfieldz,zvalue,
                       X,Y,eps_Oxy,Xfield,Yfield,num_periods,show_fig)
    
    Plot_ImHfield_Profile(Hfieldx,Hfieldy,Hfieldz,zvalue,
                       X,Y,eps_Oxy,Xfield,Yfield,num_periods,show_fig)

##### Run the MAIN program
if __name__ == "__main__":
    main()