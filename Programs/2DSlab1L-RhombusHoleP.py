import numpy as np 
import meep as mp 
from meep import mpb 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab1LRhombusHole 
from DielectricProfile import *
from BandStructure import *  
from ExportData import * 
from EFieldProfile import *
from HFieldProfile import *
from Materials import * 
from LightCone import LightCone 

######################################################################################
#                                                                                    #
#       Calculation of the band structure for 2DSlab1L with rhombus hole             #
#       The hole is deformed from the C4 square, breaking the C4 symmetry            #
#       The product of the diagonal remains unchanged, so we add the suffix          #
#       P (product). It allows to keep the hole area, filling factor, and            #
#       effective refractive index unchanged                                         # 
#       That means let d be the diagonal of the square hole, the diagonals           #
#       of the rhombi are:                                                           #
#           d1 =  d*(1+e)/(1-e)                                                      #
#           d2 =  d*(1-e)/(1+e)                                                      #
#                                                                                    #
######################################################################################

##### The main program goes here 
def main():
    ##################################################################################
    #                                                                                # 
    #             Please define the parameters of the simulation here                #
    #                                                                                #
    ##################################################################################

    ### The light polarization 
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'all'
    print('# Polarization = '+polarization)

    ### Resolution 
    resolution = mp.Vector3(32,32,32)  # pixels/a 
    print('# The resolution:'+str(resolution))

    ### Number of bands 
    num_bands = 15
    print('# The number of bands to simulate: '+str(num_bands))

    ### Geometrical parameters 
    h = 0.35       # Thickness of the slab 
    Lz = 5.0       # The height of the unit cell along the z-direction 
    b = 0.38       # The mean of the projections of the diagonals vertices (0<= b <0.5
    e = -0.1        # The anisotropy between the two diagonals (-1<= e <=1)

    print('# Thickness of the slab h = '+str(h))
    print('# The projection of mean half-axis of the rhombus on the edges b = '+str(b))
    print('# The anisotropy between the two diagonals e = '+str(e))
    print('# The height of the simulation cell Lz = '+str(Lz))

    ### The vertices of the rhombus hole 
    # If e > 0 then the diagonal x=y is longer than the diagonal x=-y 
    # If e = 0 then the hole is a square 
    # If e < 0 then the diagonal x=y is shorter than the diagonal x=-y 
    vertices = [
        mp.Vector3( 0.5*b*(1+e)/(1-e),  0.5*b*(1+e)/(1-e), 0.0 ),
        mp.Vector3( 0.5*b*(1-e)/(1+e), -0.5*b*(1-e)/(1+e), 0.0 ),
        mp.Vector3(-0.5*b*(1+e)/(1-e), -0.5*b*(1+e)/(1-e), 0.0 ),
        mp.Vector3(-0.5*b*(1-e)/(1+e),  0.5*b*(1-e)/(1+e), 0.0 )
    ]

    ### The k-point at which we plot the field profile 
    k_field = mp.Vector3(0.5,0.5,0.0)   # M-point 
    print('# The k-point at which we plot the field profile:'+str(k_field))

    ### Show figure (Yes/No)
    show_fig = 'Yes'
    print('# Show the figure: '+str(show_fig))

    ### Materials 
    # The material forming the structure 
    Mater = Si 

    # The envionment 
    Envir = PMMA 

    # The part of the momentum space to plot the band structure 
    # Choose between: 'BZ', 'Gamma-M', 'M', 'M-vicinity'
    kSpace = 'BZ'
    print('# The part of the momentum space to simulate:'+kSpace) 

    ### Number of k-points to interpolate between the 2 high-symmetry points 
    Nk = 19 
    print('# The number of points to interpolate the high-symmetry line Nk = '+str(Nk))

    ### The set of k-points 
    # kSpace = BZ 
    k_points = [
        mp.Vector3(0.0,0.0,0.0),    # Gamma 
        mp.Vector3(0.5,0.0,0.0),    # X 
        mp.Vector3(0.5,0.5,0.0),    # M+ 
        mp.Vector3(0.0,0.0,0.0),    # Gamma  
        mp.Vector3(-0.5,0.5,0.0)    # M-  
    ]

    # kSpace = M-vinity 
    #k_points = [
    #    mp.Vector3(0.5,0.5,0.0),    # M 
    #    mp.Vector3(0.45,0.45,0.0),  # Gamma'
    #    mp.Vector3(0.5,0.45,0.0),   # X'
    #    mp.Vector3(0.5,0.5,0.0),    # M 
    #    mp.Vector3(0.55,0.5,0.0)    # Gamma'
    #]

    k_points = mp.interpolate(Nk,k_points)

    # Show the figure
    show_fig = 'No'

    ### Calculate the dispersion for the light cone with Envir and k_points
    lightcone = LightCone(Envir,k_points,resolution)

    ### The k-point at which we plot the field profile
    k_field = mp.Vector3(0.5,0.5,0.0)   # M-point
    print('# The k-point at which we plot the field profile:'+str(k_field))

    ### The parameters to calculate the E- and H-field profiles
    ### in the planes parallel to Oxy
    Bloch_Phase = 'True'   # Enable the Bloch phase: True/False
    resolution_eps = 81    # Number of pixels per a 
    resolution_field = 81  # Number of pixels per a 
    num_periods = 3        # Number of periods along each direction  

    # The value of z where we take a slice to plot the field
    zvalue = 0.25*h

    ############################################################################
    #                                                                          #
    #                       Here the simulation starts                         #
    #                                                                          #
    ############################################################################

    ### Define the mode solver 
    ms = _2DSlab1LRhombusHole(h,Lz,b,e,vertices,
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

    ### The number of elements in k_points 
    number = np.arange(len(ms.k_points))

    ### The title and the name of the files
    namesave = '2DSlab1L-RhombusHoleP-h_'+str(h)+'-b_'+str(b) \
        +'-e_'+str(e)+'-'+polarization 
    
    ### Output the dielectric profile along the z-direction with projected (X,Y)
    x = 0.499 
    y = 0.499 
    zmin = -0.5*Lz 
    zmax = 0.5*Lz 
    Nz = 50 

    # Calculate the dielectric profile along the z-direction at fixed (x,y)
    z_array,epsilon_z_array = DielectricProfileZ(ms,x,y,zmin,zmax,Nz)

    # Plot the dielectric profile, the name of the figure is:
    #           namesave+'-epsilon-z.png'
    PlotDielectricProfileZ(x,y,z_array,epsilon_z_array,namesave,show_fig)

    ### Output the dielectric profile with Ncellx x Ncelly unit cells 
    ### for zmin <= z <= zmax, with Nz values of z 
    Ncellx = 5
    Ncelly = 5 
    zmin = -0.2*Lz 
    zmax =  0.2*Lz 
    Nx = 300 
    Ny = 300 
    Nz = 9 

    # Calculate the dielectric profile in planes parallel to Oxy 
    x_plot,y_plot,z_array,epsilon_xy_array \
        = DielectricProfileXY(ms,Ncellx,Ncelly,zmin,zmax,Nx,Ny,Nz)
    
    # Plot the dielectric profile, the name of the figure is:
    #       namesave+'-z_'+str(k)+'.png'
    # where k is the number of the value of z in the array z_array 
    PlotDielectricProfileXY(x_plot,y_plot,z_array,epsilon_xy_array,namesave,show_fig)

    # Print the dielectric profile to the file:
    #       namesave+'-epsilon-xy.txt'
    PrintDielectricProfileXY(x_plot,y_plot,z_array,epsilon_xy_array,namesave)



    ### Print the band structure to file 
    #print(np.shape(lightcone))
    #print(lightcone)

    PrintBandStructure(freqs,number,namesave)

    ### Plot the band structure 
    if kSpace == 'BZ':
        PlotBand_BrillouinZone_Scell_Rhole(number,freqs,Nk,lightcone,namesave,show_fig)
    elif kSpace == 'M-vicinity':
        PlotBand_M_Scell_Rhole(number,freqs,Nk,lightcone,namesave,show_fig)
    else:
        print('ERROR! The k-point has not been in the allowed list yet')
        exit()

    # Define the mode solver to calculate the E-field
    Efieldx,Efieldy,Efieldz,X,Y,Xfield,Yfield,eps_Oxy = EFields_2DSlab1LRhombusHoleP(h,Lz,b,e,vertices,
                                num_bands,resolution,Mater,Envir,
                                k_field,zvalue,
                                polarization,resolution_eps,resolution_field,
                                num_periods,Bloch_Phase)
    
    # Plot the E-field
    Plot_Efield_Profile(Efieldx,Efieldy,Efieldz,zvalue,
                       X,Y,eps_Oxy,Xfield,Yfield,num_periods,show_fig)
    
    # Define the mode solver to calculate the H-field
    Hfieldx,Hfieldy,Hfieldz,X,Y,Xfield,Yfield,eps_Oxy = HFields_2DSlab1LRhombusHoleP(h,Lz,b,e,vertices,
                                num_bands,resolution,Mater,Envir,
                                k_field,zvalue,
                                polarization,resolution_eps,resolution_field,
                                num_periods,Bloch_Phase)

    # Plot the H-field
    Plot_Hfield_Profile(Hfieldx,Hfieldy,Hfieldz,zvalue,
                       X,Y,eps_Oxy,Xfield,Yfield,num_periods,show_fig)

##### Run the MAIN program 
if __name__ == "__main__":
    main() 