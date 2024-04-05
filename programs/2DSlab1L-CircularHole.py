import numpy as np
import meep as mp 
from meep import mpb

import sys
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab1LCircularHole
from DielectricProfile import *
from BandStructure import *
from ExportData import *


##### The MAIN program goes here 
def main():
    ### The light polarization
    # Choose between: 'all', 'zeven', 'zodd' 
    polarization = 'all'

    ### The part of the momentum space to plot the band structure
    # Choose between: 'BZ', 'M' 
    kSpace = 'BZ'

    ### Resolution 
    resolution = mp.Vector3(8,8,8)   # pixels/a

    ### Geometrical parameters
    h = 0.40     # Thickness of the slab
    radius = 0.4 # Radius of the hole
    Lz = 5.0     # The height of the unit cell along the z-direction

    ### Number of bands
    num_bands = 8 

    ### Number of k-points to interpolate between 2 high-symmetry points
    Nk = 9 

    ### Define the mode solver
    ms = _2DSlab1LCircularHole(h,Lz,radius,num_bands,Nk,resolution,kSpace)

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
    namesave = '2DSlab1L-CircularHole-h_'+str(h)+'-r_'+str(radius) \
        +'-'+polarization+'-'+kSpace 
    
    ### Output the dielectric profile along the z-direction with projected (X,Y)
    x = 0.27 
    y = 0.36 
    zmin = -0.5*Lz
    zmax = 0.5*Lz 
    Nz = 500 

    # Calculate the dielectric profile along the z-direction at fixed (x,y)
    z_array,epsilon_z_array = DielectricProfileZ(ms,x,y,zmin,zmax,Nz)

    # Plot the dielectric profile, the name of the figure is:
    #           namesave+'-epsilon-z.png'
    PlotDielectricProfileZ(x,y,z_array,epsilon_z_array,namesave)

    # Print the dielectric profile to the file:
    #           namesave+'-epsilon-z.txt'
    PrintDielectricProfileZ(x,y,z_array,epsilon_z_array,namesave)

    ### Output the dielectric profile with Ncellx x Ncelly unit cells 
    ### for zmin <= z <= zmax, with Nz values of z 
    Ncellx = 5
    Ncelly = 5
    zmin = -0.5*Lz
    zmax = 0.5*Lz 
    Nz = 9
    
    ### Print the band structure to file
    PrintBandStructure(freqs,number,namesave)

    ### Plot the band structure
    if kSpace == 'BZ':
        PlotBand_BrillouinZone(number,freqs,Nk,namesave)
    elif kSpace == 'M':
        PlotBand_M(number,freqs,Nk,namesave)
    else:
        print('ERROR! The k-point has not been in the allowed list yet')
        exit()

##### Run the MAIN program
if __name__ == "__main__":
    main()