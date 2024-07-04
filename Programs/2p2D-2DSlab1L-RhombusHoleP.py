import numpy as np 
import meep as mp 
from meep import mpb 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab1LRhombusHole 
from DielectricProfile import * 
from BandStructure import * 
from ExportData import * 
from Materials import * 
from LightCone import LightCone 

import os 

##### The MAIN program goes here 
def main():
    ##################################################################################
    #                                                                                #
    #           Please define the parameters of the simulation here                  #
    #                                                                                #
    ##################################################################################

    ### The light polarization
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'zeven'
    print('# Polarization = '+polarization)

    ### Resolution 
    resolution = mp.Vector3(8,8,8)  # pixels/a 
    print('# The resolution:'+str(resolution))

    ### Number of bands 
    num_bands = 10 
    print('# The number of bands to simulate: '+str(num_bands))
    
    ### Geometrical parameters 
    h = 0.35        # Thickness of the slab 
    Lz = 5.0        # The height of the unit cell along the z-direction 
    b = 0.38        # The mean of the projections of the diagonal vertices (0 <= b < 0.5)
    
    print('# Thickness of the slab h = '+str(h))
    print('# The edge of the undeformed square hole b = '+str(b))
    print('# The height of the simulation cell Lz = '+str(Lz))

    ### Number of k-points to interpolate between the 2 high-symmetry points 
    Nk = 9 
    print('# The number of points to interpolate the high-symmetry line Nk = '+str(Nk))

    ### The array of synthetic momenta
    N1 = 5  
    k1_array = np.linspace(0.0,0.5,N1)

    N2 = 5 
    k2_array = np.linspace(0.0,0.5,N2)

    # The part of the momentum space to plot the band structure 
    # Choose between: 'BZ', 'Gamma-M', 'M', 'M-vicinity'
    kSpace = 'BZ'
    print('# The part of the momentum space to simulate:'+kSpace)

    ### The set of k-points 
    # kSpace = BZ
    k_points = [
        mp.Vector3(0.0,0.0,0.0),    # Gamma 
        mp.Vector3(0.5,0.0,0.0),    # X 
        mp.Vector3(0.5,0.5,0.0),    # M+
        mp.Vector3(0.0,0.0,0.0),    # Gamma 
        mp.Vector3(-0.5,0.5,0.0)    # M- 
    ]

    # kSpace = M-vicinity 
    #k_points = [
    #    mp.Vector3(0.5,0.5,0.0),    # M 
    #    mp.Vector3(0.45,0.45,0.0),  # Gamma'
    #    mp.Vector3(0.5,0.45,0.0),   # X'
    #    mp.Vector3(0.5,0.5,0.0),    # M 
    #    mp.Vector3(0.55,0.5,0.0)    # Gamma' 
    #] 

    k_points = mp.interpolate(Nk,k_points)

    ### The number of elements in k_points 
    number = np.arange(len(k_points))

    ### Show the figures (Yes/No)
    show_fig = 'No'
    print('# Show the figures: '+str(show_fig)) 

    ### Materials 
    # The material forming the structure 
    Mater = Si 

    # The environment
    Envir = PMMA 

    ### Calculate the dispersion for the light cone with Envir and k_points 
    lightcone = LightCone(Envir,k_points,resolution)

    #################################################################################
    #                                                                               #
    #                       Here the simulation starts                              #
    #                                                                               #
    #################################################################################

    ##### We scan over the synthetic momenta. For each synthetic momentum, 
    ###   we do a simulation. And then we write the results to files
    for i1 in range(N1):
        for i2 in range(N2):
            ### The values of the synthetic momenta 
            k1 = k1_array[i1]
            k2 = k2_array[i2]

            ### The deformation of the square into the rhombus 
            e = 0.1*np.sin(2*np.pi*k1)*np.cos(2*np.pi*k2) 

            ### Define the vertices of the rhombus hole 
            vertices = [
                mp.Vector3( 0.5*b*(1+e)/(1-e),  0.5*b*(1+e)/(1-e), 0.0),
                mp.Vector3( 0.5*b*(1-e)/(1+e), -0.5*b*(1-e)/(1+e), 0.0),
                mp.Vector3(-0.5*b*(1+e)/(1-e), -0.5*b*(1+e)/(1-e), 0.0),
                mp.Vector3(-0.5*b*(1-e)/(1+e),  0.5*b*(1-e)/(1+e), 0.0)
            ] 

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
            
            ### The title and the name of the files 
            namesave = '2p2D-2DSlab1L-RhombusHoleP-h_'+str(h)+'-b_'+str(b) \
                + "-k1_{0:.4f}".format(k1)+"-k2_{0:.4f}".format(k2) \
                + "-e_{0:.4f}".format(e)+'-'+polarization 
            
            print(namesave) 

            os.mkdir('./'+str(namesave))
            os.chdir('./'+str(namesave))

            ### Output the dielectric profile along the z-direction with projected (X,Y)
            x = 0.499 
            y = 0.499 
            zmin = -0.5*Lz 
            zmax = 0.5*Lz
            Nz = 50 

            # Calculate the dielectric profile along the z-direction at fixed (x,y)
            z_array,epsilon_z_array = DielectricProfileZ(ms,x,y,zmin,zmax,Nz) 

            # Plot the dielectric profile, the name of the figure is:
            #               namesave + '-epsilon-z.png'
            PlotDielectricProfileZ(x,y,z_array,epsilon_z_array,namesave,show_fig)

            ### Output the dielectric profile with Ncellx x Ncelly unit cells 
            ### for zmin <= z <= zmax, with Nz values of z 
            Ncellx = 5
            Ncelly = 5
            zmin = -0.2*Lz 
            zmax = 0.2*Lz 
            Nx = 200
            Ny = 200 
            Nz = 5 

            # Calculate the dielectric profile in planes parallel to Oxy 
            x_plot,y_plot,z_array,epsilon_xy_array \
                = DielectricProfileXY(ms,Ncellx,Ncelly,zmin,zmax,Nx,Ny,Nz)
            
            # Plot the dielectric profile, the name of the figure is:
            #       namesave + '-z_'+str(k)+'.png'
            # where k is the number of the value of z in the array z_array 
            PlotDielectricProfileXY(x_plot,y_plot,z_array,epsilon_xy_array,namesave,show_fig)

            # Print the dielectric profile to the file:
            #       namesave + '-epsilon-xy.txt'
            PrintDielectricProfileXY(x_plot,y_plot,z_array,epsilon_xy_array,namesave)


            ### Print the band structure to file 
            PrintBandStructure(freqs,number,namesave) 

            ### Plot the band structure 
            if kSpace == 'BZ':
                PlotBand_BrillouinZone_Scell_Rhole(number,freqs,Nk,lightcone,namesave,show_fig)
            elif kSpace == 'M-vicinity':
                PlotBand_M_Scell_Rhole(number,freqs,Nk,lightcone,namesave,show_fig)
            else:
                print('ERROR! The k-point has not been in the allowed list yet')
                exit()

            ### Move to the initial directory 
            os.chdir('..')

##### Run the MAIN program 
if __name__ == '__main__':
    main()