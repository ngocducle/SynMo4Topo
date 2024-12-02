import numpy as np 
import meep as mp 
from meep import mpb 

import sys
sys.path.insert(0,'../src/')
from ModeSolvers import _2DSlab2LRhombusHole
from DielectricProfile import * 
from ExportData import *
from Materials import *
from LightCone import LightCone 

##### The main program goes here 
def main():
    #################################################################################
    #                                                                               #
    #       Band structure of 2D slab 2L- rhombus hole, case (1+1)D                 #
    #       This code studies a 2DSlab2L rhombus holes with:                        #
    #       Slab thickness h, the parent square holes has edge length b             #
    #       and the deformation parameter of the rhombus is e                       #
    #                                                                               #
    #       The hole is deformed from the C4 square, breaking the C4 symmetry       #
    #       The product of the diagonal remains unchanged, so we add the suffix     #
    #       P (product). It allows to keep the hole area, filling factor, and       #
    #       effective refractive index unchanged                                    # 
    #       That means let d be the diagonal of the square hole, the diagonals      #
    #       of the rhombi are:                                                      #
    #                   d1 =  d*(1+e)/(1-e)                                         #
    #                   d2 =  d*(1-e)/(1+e)                                         #
    #                                                                               #
    #       The two slabs are displaced by distances delta1 and delta2              #
    #       along the x and y directions, respectively                              #
    #       In the vicinity of the M point M(0.5,0.5)                               #
    #       From the M point:                                                       #
    #               Genuine momenta: kx = ky = k                                    #
    #               Synthetic momenta: qx = qy = q                                  #
    #                                                                               #
    #       Please define the parameters of the simulation here                     #
    #                                                                               #
    #################################################################################

    ### Resolution 
    resolution = mp.Vector3(16,16,16)      # pixels/a 
    print('# Resolution = '+str(resolution))

    ### Number of bands 
    num_bands = 25
    print('# The number of bands to simulate: '+str(num_bands))

    ### Geometrical parameters 
    # The upper layer 
    h1 = 0.35 # Thickness of the upper layer 
    b1 = 0.38  # The edge length of the undeformed square hole 
    e1 = 0.0  # The anisotropy between the two diagonals of the upper layer 

    print('\n# Upper slab:')
    print('# Thickness h1 = '+str(h1))
    print('# The edge length of the undeformed square hole b1 = '+str(b1))
    print('# The anisotropy between the two diagonals e1 = '+str(e1)) 

    # The lower layer 
    h2 = 0.35 # Thickness of the lower layer 
    b2 = 0.38  # The edge length of the undeformed square hole 
    e2 = 0.0  # The anisotropy between the two diagonals of the lower layer 

    print('\n# Lower slab:')
    print('# Thickness h2 = '+str(h2))
    print('# The edge length of the undeformed square hole b2 = '+str(b2))
    print('# The anisotropy between the two diagonals e2 = '+str(e2))

    # The interlayer parameters
    dist = 0.10   # Distance between the two layers 
    Lz = 10.0     # The height of the unit cell along the z-direction 

    print('\n# The distance between the two layers d = '+str(dist))
    print('# The height of the simulation cell Lz = '+str(Lz))

    ### Materials 
    # The material forming the structure 
    Mater = Si 

    # The environment 
    Envir = PMMA

    ### The number of k-points (genuine momenta) to interpolate 
    Nk = -1
    print('# The number of k-points to interpolate the high-symmetry line Nk = '+str(Nk))

    ### The set of k-points (MPB)
    Kmin = 0.00
    Kmax = 0.00 

    """k_points = [
        mp.Vector3(Kmin,Kmin,0.0),    
        mp.Vector3(Kmax,Kmax,0.0)     # M 
    ]

    k_points = mp.interpolate(Nk,k_points)"""
    k_points = [mp.Vector3(0.0,0.0,0.0)]

    ### The array of k-points (to print)
    k_array = np.linspace(Kmin,Kmax,Nk+2) 

    ### The number of q-points (synthetic momenta)
    Nq = 51

    ### The array of q-points (synthetic momenta)
    q_array = np.linspace(0.0,0.5,Nq)

    ### Show the figure 
    show_fig = 'No'

    ### The title and the name of the files 
    namesave = '1p1D-2DSlab2L-RHoleP-h1_'+str(h1)+'-b1_'+str(b1)+'-e1_'+str(e1) \
        + '-h2_'+str(h2)+'-b2_'+str(b2)+'-e2_'+str(e2) \
        + '-d_'+str(dist)+'_kq'
    
    ### Calculate the dispersion for the light cone with Envir and k_points
    print('# Calculate the light cone dispersion')
    lightcone = LightCone(Envir,k_points,resolution)


    ################################################################################
    #                                                                              #
    #                        Here the simulation starts                            #
    #                                                                              #
    ################################################################################

    ##### The array of frequencies 
    #Bands = np.zeros((Nk*Nq,num_bands))

    ##### We scan over the array of q-points
    for iq in range(Nq):
        ### The synthetic momentum 
        q = q_array[iq]

        ### The shift along the x-direction 
        deltax = q 

        ### The shift along the y-direction 
        deltay = q 

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

        ### Define the mode solver 
        ms = _2DSlab2LRhombusHole(h1,b1,e1,vertices1,
                                  h2,b2,e2,vertices2,
                                  Lz,dist,deltax,deltay,
                                  num_bands,k_points,resolution,Mater,Envir)
        
        ### Run the simulation 
        ms.run()

        ### Extract the frequencies of the modes from the ModeSolver 
        freqs = ms.all_freqs 
        #print(freqs)

        ### Save the frequencies to the array Bands 
        #for ik in range(Nk):
        #    Bands[iq*Nq+ik,:] = freqs[ik,:]

        ### Print the band structure 
        namesavetxt = namesave+'-Band'
        namesavetxt = namesavetxt + "-q__{0:.4f}".format(q) 
        namesavetxt = namesavetxt + '.txt'

        with open(namesavetxt,'w') as file:  
            for ik in range(Nk+2):
                file.write('%.8f    ' % q)
                file.writelines('%.8f    ' % k_array[ik])
                file.writelines('%.8f    ' % w for w in freqs[ik])
                file.write('\n')

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
        #x_plot,y_plot,z_array,epsilon_xy_array \
        #    = DielectricProfileXY(ms,Ncellx,Ncelly,zmin,zmax,Nx,Ny,Nz)
    
        # Plot the dielectric profile, the name of the figure is:
        #           namesave+'-z_'+str(k)+'.png'
        # where k is the number of the value of z in the array z_array
        #PlotDielectricProfileXY(x_plot,y_plot,z_array,epsilon_xy_array,namesave+str(iq),show_fig)
    
        # Print the dielectric profile to the file:
        #           namesave+'-epsilon-xy.txt'
        #PrintDielectricProfileXY(x_plot,y_plot,z_array,epsilon_xy_array,namesave+str(iq))


    ##### Print the band structure 
    #with open(namesave+'-Band.txt','w') as file: 
    #    for iq in range(Nq):
    #        for ik in range(Nk):
    #            file.write('%.8f        ' % q_array[iq])
    #            file.writelines('%.8f   ' % k_array[ik])
    #            file.writelines('%.8f   ' % w for w in Bands[iq*Nq+ik])
    #            file.write('\n')


if __name__ == '__main__':
    main()
