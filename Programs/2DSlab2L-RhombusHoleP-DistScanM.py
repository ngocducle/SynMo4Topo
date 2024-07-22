import numpy as np 
import meep as mp 
from meep import mpb
import matplotlib.pyplot as plt 

import sys 
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab2LRhombusHole
from Materials import * 

#################################################################################
#                                                                               #
#       The program scans the 2DSlab2L with varying interlayer distance         #
#       and calculate the energy eigenvalues at the M-point                     #
#       The shifts deltax = deltay = delta                                      #
#                                                                               #
#################################################################################

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
    #############################################################################
    #                                                                           #
    #       Please define the parameters of the simulation here                 #
    #                                                                           #
    #############################################################################
    
    ### The light polarization
    # Choose between: 'all', 'zeven', 'zodd'
    polarization = 'all'
    print('# Polarization = '+polarization)

    ### We focus on the M-point
    kSpace = 'M'
    k_points = [mp.Vector3(0.5,0.5,0.0)]
    print('# We focus on the M-point')

    ### Resolution 
    resolution = mp.Vector3(16,16,16)   # pixels/a
    print('# The resolution:'+str(resolution))

    ### Geometrical parameters 
    # The upper layer 
    h1 = 0.35   # Thickness of the upper layer
    b1 = 0.38    # The edge length of the undeformed square hole 
    e1 = 0.10   # The anisotropy between the two diagonals of the upper layer
    
    print('\n# Upper slab:')
    print('# Thickness h1 = '+str(h1))
    print('# The edge length of the undeformed square hole b1 = '+str(b1))
    print('# The anisotropy between the two diagonals e1 = '+str(e1))

    # The lower layer 
    h2 = 0.35   # Thickness of the lower layer
    b2 = 0.38   # The edge length of the undeformed square hole    
    e2 = 0.10   # The anisotropy between the two diagonals of the lower layer

    # The interlayer parameters 
    dist = 0.10  # Distance between the two layers 
    deltax = 0.0 # Relative displacement along the x-direction 
    deltay = 0.0 # Relative displacement along the y-direction 
    Lz = 10.0    # The height of the unit cell along the z-direction 

    print('\n# The distance between the two layers d = '+str(dist))
    print('# Relative displacement along x-direction deltax = '+str(deltax))
    print('# Relative displacement along y-direction deltay = '+str(deltay))

    ### The vertices of the rhombus holes
    vertices1 = [
        mp.Vector3( 0.5*b1*(1+e1)/(1-e1)+0.5*deltax,  
                    0.5*b1*(1+e1)/(1-e1)+0.5*deltay, 0.0 ),
        mp.Vector3( 0.5*b1*(1-e1)/(1+e1)+0.5*deltax, 
                   -0.5*b1*(1-e1)/(1+e1)+0.5*deltay, 0.0 ),
        mp.Vector3(-0.5*b1*(1+e1)/(1-e1)+0.5*deltax, 
                   -0.5*b1*(1+e1)/(1-e1)+0.5*deltay, 0.0 ),
        mp.Vector3(-0.5*b1*(1-e1)/(1+e1)+0.5*deltax,  
                    0.5*b1*(1-e1)/(1+e1)+0.5*deltay, 0.0 ) 
    ]

    vertices2 = [
        mp.Vector3( 0.5*b2*(1+e2)/(1-e2)-0.5*deltax,
                    0.5*b2*(1+e2)/(1-e2)-0.5*deltay, 0.0 ),
        mp.Vector3( 0.5*b2*(1-e2)/(1+e2)-0.5*deltax, 
                   -0.5*b2*(1-e2)/(1+e2)-0.5*deltay, 0.0 ),
        mp.Vector3(-0.5*b2*(1+e2)/(1-e2)-0.5*deltax, 
                   -0.5*b2*(1+e2)/(1-e2)-0.5*deltay, 0.0 ),
        mp.Vector3(-0.5*b2*(1-e2)/(1+e2)-0.5*deltax,  
                    0.5*b2*(1-e2)/(1+e2)-0.5*deltay, 0.0 ) 
    ]


    ### Number of bands 
    num_bands = 30 
    print('# The number of bands to simulate: '+str(num_bands))

    ### Show figure (Yes/No)
    show_fig = 'Yes'
    print('\n# Show the figure: '+str(show_fig))

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
    Ndist = 21
    DistArray = np.linspace(0,0.2,Ndist)

    ### Initialize the arrays of bands
    Bands = np.zeros((len(DistArray),num_bands))

    for i in range(len(DistArray)):

        ### Extract the interlayer distance
        dist = DistArray[i]

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

        Bands[i,:] = freqs 

    ##### Print the photonic bands into a file
    ### Column 0: DistArray 
    ### Column j (1 <= j <= num_bands): Band j

    namesave = '2DSlab2L-RHoleP-DistScanM-'+polarization

    with open(namesave+'.txt','w') as file:
        file.write('# Interlayer distance / a  Band1    Band2   ...')
        file.write('\n')

        for i in range(len(DistArray)):
            file.write('%.8f        ' % DistArray[i])

            for j in range(num_bands):
                file.writelines('%.8f       ' % Bands[i,j])

            file.write('\n')

    ##### Plot the figure
    PlotFigureM(DistArray,Bands,num_bands,namesave,show_fig)

##### Run the MAIN program
if __name__ == "__main__":
    main()