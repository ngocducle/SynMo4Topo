import numpy as np 
import meep as mp 
from meep import mpb 
import matplotlib.pyplot as plt

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab2LCircularHole 
from Materials import * 

##### FUNCTION: Plot the figure
def PlotFigure2p1D(delta_Array,Bands,num_bands,namesave,show_fig):
    fig,ax = plt.subplots()

    for j in range(num_bands):
        plt.plot(delta_Array,Bands[:,j],'o',markerfacecolor='black',markersize=1)

    ax.set_xlabel('delta/a',fontsize=14)
    ax.set_ylabel(r'$\omega a / (2 \pi c)$',fontsize=14)
    plt.savefig(namesave+'.png')
    if show_fig == 'Yes':
        plt.show()


##### The MAIN program goes here
def main():
    #################################################################################
    #                                                                               #
    #             Please define the parameters of the simulation here               #
    #                                                                               #
    #################################################################################

    ### The light polarization
    #   Choose between: 'all', 'zeven', 'zodd'
    polarization = 'all'
    print('# Polarization = '+polarization)

    ### The part of the 2D momentum space to plot the band structure 
    #   Choose between: 'BZ', 'M', 'M-vicinity'
    kSpace = 'M'
    print('# The part of the momentum space to simulate:'+kSpace) 

    ### Resolution 
    resolution = mp.Vector3(16,16,16)   # pixels/a
    print('# The resolution:'+str(resolution))

    ### Number of bands 
    num_bands = 25
    print('# The number of bands to simulate: '+str(num_bands))

    ### Number of momenta/synthetic momenta to interpolate 
    ### between 2 high-symmetry points
    Nk = 29
    print('# The number of points to interpolate between 2 high-symmetry points Nk = '+str(Nk))

    ### Array of synthetic momenta
    delta_Array = np.linspace(0.0,0.5,Nk)
    print('# The array of synthetic momenta: ')
    print(delta_Array)

    ### Show the figure (Yes/No)
    show_fig = 'Yes'
    print('# Show the figure: '+str(show_fig))

    ### Initiate the array of frequencies
    Bands = np.zeros((Nk,num_bands))
    #print('# The array of frequencies: ')
    #print(freq_arrays)

    ##### We scan over the delta_array
    for i in np.arange(Nk):
        ### Extract the synthetic momentum from delta_array 
        delta = delta_Array[i]
        print('# j = '+str(i)+': delta = '+str(delta)) 

        ### Geometrical parameters
        h = 0.3             # Thickness of the slab
        radius = 0.4        # Radius of the hole 
        dist = 0.05         # Distance between the two layers 
        delta1 = delta      # Relative displacement along the x-direction
        delta2 = delta      # Relative displacement along the y-direction 
        Lz = 10.0           # The height of the unit cell along the z-direction 

        print('# Thickness of the slab h = '+str(h))
        print('# Radius of teh circular hole r = '+str(radius))
        print('# The height of the simulation cell Lz = '+str(Lz))
        print('# Relative shift along x-direction delta1 = '+str(delta1))
        print('# Relative shift along y-direction delta2 = '+str(delta2))
        print('# Interlayer distance dist = '+str(dist))

        ### Materials
        #   The material forming the structure 
        Mater = Si

        #   The environment
        Envir = PMMA

        ### Define the mode solver 
        ms = _2DSlab2LCircularHole(h,Lz,radius,dist,delta1,delta2,
                                   num_bands,1,resolution,kSpace,Mater,Envir)
    
        print(ms.k_points)

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

        Bands[i,:] = freqs 

    ##### Print the photonic bands into a file
    ###   Column 0: delta_Array 
    ###   Column j (1 <= j <= num_bands): Band j 
    namesave = '2p1D-CircularHole-'+kSpace+'-'+polarization

    with open(namesave+'.txt','w') as file:
        file.write('# Shift / a     Band 1      Band 2      ...')
        file.write('\n')

        for i in range(len(delta_Array)):
            file.write('%.8f           ' % delta_Array[i])

            for j in range(num_bands):
                file.writelines('%.8f           ' % Bands[i,j])

            file.write('\n')

    ##### Plot the figure
    PlotFigure2p1D(delta_Array,Bands,num_bands,namesave,show_fig)

##### Run the MAIN program
if __name__ == '__main__':
    main()