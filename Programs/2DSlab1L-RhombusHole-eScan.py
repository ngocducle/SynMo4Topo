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
    #             Please define the parameters of the simulation here                #
    #                                                                                #
    ##################################################################################

    ### The light polarization 
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'zeven'
    print('# Polarization = '+polarization)

    ### Resolution 
    resolution = mp.Vector3(8,8,8)   # pixels/a 
    print('# The resolution:'+str(resolution))

    ### Number of bands 
    num_bands = 25 
    print('# The number of bands to simulate: '+str(num_bands))

    ### The thickness of the slab
    h = 0.35 

    ### The edge length of the square hole when the two diagonals are equal 
    b = 0.4  

    ### The height of the unit cell along the z-direction 
    Lz = 5.0

    ### The array of the anisotropy between the two diagonals (-1 <= e <= 1)
    Ne = 6
    e_array = np.linspace(-0.5,0.5,Ne)

    ### Show the figure (Yes/No)
    show_fig = 'No'
    print('# Show the figure: '+str(show_fig))

    ### Materials 
    # The material forming the structure 
    Mater = Si 

    # The environment
    Envir = PMMA 

    # The part of the momentum space to plot the band structure 
    # Choose between: 'BZ', 'Gamma-M', 'M', 'M-vicinity'
    kSpace = 'BZ'
    print('# The part of the momentum space to simulate:'+kSpace) 

    ### Number of k-points to interpolate between 2 high-symmetry points 
    Nk = 19 
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

    ### Calculate the dispersion for the light cone with Envir and k_points 
    lightcone = LightCone(Envir,k_points,resolution)

    ############################################################################
    #                                                                          #
    #                       Here the simulation starts                         #
    #                                                                          #
    ############################################################################
    for e in e_array:
        ### Geometrical parameters 
        print('# Thickness of the slab h = '+str(h))
        print('# The edge length of the parent square is b = '+str(b))
        print('# The anisotropy between the two diagonals e = '+str(e))
        print('# The height of the simulation cell Lz = '+str(Lz))

        ### The vertices of the rhombus hole
        # If e>0 then the diagonal x=y is longer than the diagonal x=-y
        # If e=0 then the hole is a square 
        # If e<0 then the diagonal x=y is shorter than the diagonal x=-y 
        vertices = [
            mp.Vector3(  0.5*(1+e)*b,  0.5*(1+e)*b, 0.0 ),
            mp.Vector3(  0.5*(1-e)*b, -0.5*(1-e)*b, 0.0 ),
            mp.Vector3( -0.5*(1+e)*b, -0.5*(1+e)*b, 0.0 ),
            mp.Vector3( -0.5*(1-e)*b,  0.5*(1-e)*b, 0.0 )
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

        ### The number of elements in k_points 
        number = np.arange(len(ms.k_points))

        ### The title and the name of the files 
        #namesave = '2DSlab1L-RhombusHole-h_'+str(h)+'-b_'+str(b) \
        #        +'-e_'+str(e)+'-'+polarization 
        namesave = "2DSlab1L-SquareHole-h_{0:.4f}".format(h)
        namesave = namesave + "-b_{0:.4f}".format(b)
        namesave = namesave + "-e_{0:.4f}".format(e)
        namesave = namesave+'-'+polarization 
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
        PrintBandStructure(freqs,number,namesave)

        ### Plot the band structure 
        if kSpace == 'BZ':
            PlotBand_BrillouinZone_Scell_Rhole(number,freqs,Nk,lightcone,namesave,show_fig)
        elif kSpace == 'M-vicinity':
            PlotBand_M(number,freqs,Nk,namesave,show_fig)
        else:
            print('ERROR! The k-point has not been in the allowed list yet')
            exit()   

        ### Move to the initial directory 
        os.chdir('..')

##### Run the MAIN program 
if __name__ == "__main__":
    main()