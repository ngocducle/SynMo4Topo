import numpy as np 
import meep as mp 
from meep import mpb 
import matplotlib.pyplot as plt 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab2LCircularHole
from Materials import * 

##############################################################################
#                                                                            #
#    The program scans the 2DSlab2L with varying interlayer distance         #
#    and calculate the energy eigenvalues at the M-point                     # 
#    The shifts deltax = deltay = 0                                          #
#                                                                            #
##############################################################################

##### FUNCTION: Plot the figure
def PlotFigureM(DistArray,Bands,num_bands,namesave,show_fig):
    fig,ax = plt.subplots()

    for j in range(num_bands):
        plt.plot(DistArray,Bands[:,j],'o',markerfacecolor='black',markersize=1)

    ax.set_xlabel('d/a',fontsize=14)
    ax.set_ylabel(r'$\omega a / (2 \pi c)$',fontsize=14)
    plt.savefig(namesave+'.png')
    if show_fig == 'Yes':
        plt.show()

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

    ### We focus on the M-point 
    kSpace = 'M'
    print('# We focus on the M-point')

    ### Resolution 
    resolution = mp.Vector3(32,32,32)   # pixels/a 
    print('# The resolution:'+str(resolution))

    ### Geometrical parameters 
    h = 0.3         # Thickness of the slab 
    radius = 0.4    # Radius of the hole 
    Lz = 10.0       # The height of the unit cell along the z-direction 

    print('# Thickness of the slab h = '+str(h))
    print('# Radius of the circular hole r = '+str(radius))
    print('# the height of the simulation cell Lz = '+str(Lz))

    ### Number of bands 
    num_bands = 20 
    print('# The number of bands to simulate: '+str(num_bands))

    ### Show figure (Yes/No)
    show_fig = 'Yes'
    print('# Show the figures: '+str(show_fig))

    ### Materials
    # The material forming the structure 
    Mater = Si 

    # The environment
    Envir = PMMA 

    ##########################################################################
    #                                                                        #
    #                   Here the simulation starts                           #
    #                                                                        #
    ##########################################################################

    ### Define the array for the interlayer distance
    Ndist = 101
    DistArray = np.linspace(0.0,2.0,Ndist)

    ### Initialize the arrays of bands
    Bands = np.zeros((Ndist,num_bands))

    for i in range(Ndist):

        ### Extract the interlayer distance 
        dist = DistArray[i]

        ### Define the mode solver 
        ms = _2DSlab2LCircularHole(h,Lz,radius,dist,0.0,0.0,
                                   num_bands,0,resolution,kSpace,Mater,Envir)
        
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
    ### Column 0: DistArray 
    ### Column j (1 <= j <= num_bands): Band j

    namesave = '2DSlab2L-CHole-DistScanM'+polarization  

    with open(namesave+'.txt','w') as file:
        file.write('# Interlayer distance / a  Band1    Band2   ...')
        file.write('\n')

        for i in range(Ndist):
            file.write('%.8f        ' % DistArray[i])

            for j in range(num_bands):
                file.writelines('%.8f       ' % Bands[i,j])

            file.write('\n')

    ##### Plot the figure
    PlotFigureM(DistArray,Bands,num_bands,namesave,show_fig)

##### Run the MAIN program 
if __name__ == "__main__":
    main()