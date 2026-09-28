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

    #################################################################################
    #                                                                               #
    #       Band structure of 2D slab 2L- rhombus hole, case (1+1)D                 #
    #       This code studies a 2DSlab2L rhombus holes with:                        #
    #       Slab thickness h, the parent square holes have edge length b            #
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
    #       The hybrid Brillouin zone band structure is calculated along the path   #
    #           Point     (k,q)      (K,shift)                                      #
    #       ----------------------------------------------------------------------  # 
    #           Gamma   (-0.5,-0.5)   (0,0)                                         # 
    #           Q       (0,-0.5)      (0.5,0)                                       # 
    #           O       (0,0)         (0.5,0.5)                                     # 
    #           K       (-0.5,0)      (0,0.5)                                       #
    #                                                                               #
    #       The path is: Gamma-Q-O-K-Gamma                                          #
    #                                                                               #
    #       Please define the parameters of the simulation here                     #
    #                                                                               #
    #################################################################################

### Resolution 
resolution = mp.Vector3(16,16,16)       # pixels/a 
print('# Resolution = '+str(resolution))

### Number of bands 
num_bands = 12
print('# The number of bands to simulate: '+str(num_bands))

### Geometrical parameters 
# The upper layer 
h1 = 0.35   # Thickness of the upper layer 
b1 = 0.38   # The edge length of the undeformed square hole 
e1 = -0.1   # The anisotropy between the two diagonals of the upper layer 

print('\n# Upper slab:')
print('# Thickness h1 = '+str(h1))
print('# The edge length of the undeformed square hole b1 = '+str(b1))
print('# The anisotropy between the two diagonals e1 = '+str(e1))

# The lower layer 
h2 = 0.35   # Thickness of the lower layer 
b2 = 0.38   # The edge length of the undeformed square hole 
e2 = -0.1   # The anisotropy between the two diagonals of the lower layer 

print('\n# Lower slab:')
print('# Thickness h2 = '+str(h2))
print('# The edge length of the undeformed square hole b2 = '+str(b2))
print('# The anisotropy between the two diagonals e2 = '+str(e2))

# The interlayer parameters 
dist = 0.10     # Distance between the two layers 
Lz = 10.0       # The height of the unit cell along the z-direction 

print('\n# The distance between the two layers d = '+str(dist))
print('# The height of the simulation cell Lz = '+str(Lz))

### Materials 
# The material forming the structure 
Mater = Si 

# The environment 
Envir = PMMA 

##### The path along the Brillouin zone 
K1 = 0.0; delta1 = 0.0    # Gamma 
K2 = 0.5; delta2 = 0.0    # Q 
K3 = 0.5; delta3 = 0.5    # O 
K4 = 0.0; delta4 = 0.5    # K 

# Number of (k,q)-points to interpolate 
Nk = 19

K_array1 = np.linspace(K1,K2,Nk+2)
K_array2 = np.linspace(K2,K3,Nk+2)
K_array3 = np.linspace(K3,K4,Nk+2)

k_array = np.concatenate((K_array1,K_array2[1:Nk+2],K_array3[1:Nk+2]))

delta_array1 = np.linspace(delta1,delta2,Nk+2)
delta_array2 = np.linspace(delta2,delta3,Nk+2)
delta_array3 = np.linspace(delta3,delta4,Nk+2)

delta_array = np.concatenate((delta_array1,delta_array2[1:Nk+2],delta_array3[1:Nk+2]))

print('k_array =')
print(k_array)
print('delta_array = ')
print(delta_array)

# Length of the path 
Npath = len(k_array)
print('Path length Npath = '+str(Npath))

### Show the figure 
show_fig = 'No'

### The title and the name of the files 
namesave = '1p1D-2DSlab2L-RHoleP-h1_'+str(h1)+'-b1_'+str(b1)+'-e1_'+str(e1) \
        + '-h2_'+str(h2)+'-b2_'+str(b2)+'-e2_'+str(e2) \
        + '-d_'+str(dist)+'_kq'

### Calculate the dispersion for the light cone with Envir and k_points 
print('# Calculate the light cone dispersion')
k_points_LC = mp.interpolate(Nk,[mp.Vector3(K1,0.0,0.0),mp.Vector3(K2,0.0,0.0)])
lightcone = LightCone(Envir,k_points_LC,resolution)

###############################################################################
#                                                                             #
#                        Here the simulation starts                           #
#                                                                             #
###############################################################################

##### Inititate the array of Bands 
Bands = []

##### We scan over the path 
for ikq in range(Npath):
    ### The synthetic momentum 
    delta = delta_array[ikq]

    ### The genunime momentum 
    k = k_array[ikq]
    k_points = [mp.Vector3(k,k,0.0)]

    ### The vertices of the rhombus holes 
    vertices1 = [
            mp.Vector3( 0.5*b1*(1+e1)/(1-e1)+0.5*delta,  
                        0.5*b1*(1+e1)/(1-e1)+0.5*delta, 0.0 ),
            mp.Vector3( 0.5*b1*(1-e1)/(1+e1)+0.5*delta, 
                       -0.5*b1*(1-e1)/(1+e1)+0.5*delta, 0.0 ),
            mp.Vector3(-0.5*b1*(1+e1)/(1-e1)+0.5*delta, 
                       -0.5*b1*(1+e1)/(1-e1)+0.5*delta, 0.0 ),
            mp.Vector3(-0.5*b1*(1-e1)/(1+e1)+0.5*delta,  
                        0.5*b1*(1-e1)/(1+e1)+0.5*delta, 0.0 ) 
        ]

    vertices2 = [
            mp.Vector3( 0.5*b2*(1+e2)/(1-e2)-0.5*delta,
                        0.5*b2*(1+e2)/(1-e2)-0.5*delta, 0.0 ),
            mp.Vector3( 0.5*b2*(1-e2)/(1+e2)-0.5*delta, 
                       -0.5*b2*(1-e2)/(1+e2)-0.5*delta, 0.0 ),
            mp.Vector3(-0.5*b2*(1+e2)/(1-e2)-0.5*delta, 
                       -0.5*b2*(1+e2)/(1-e2)-0.5*delta, 0.0 ),
            mp.Vector3(-0.5*b2*(1-e2)/(1+e2)-0.5*delta,  
                        0.5*b2*(1-e2)/(1+e2)-0.5*delta, 0.0 ) 
        ]

    ### Define the mode solver 
    ms = _2DSlab2LRhombusHole(h1,b1,e1,vertices1,
                            h2,b2,e2,vertices2,
                            Lz,dist,delta,delta,
                            num_bands,k_points,resolution,Mater,Envir)
    
    ### Run the simulation 
    ms.run()

    ### Extract the frequencies of the modes from the ModeSolver 
    freqs = ms.all_freqs 

    ### Append freqs to Band 
    Bands.append(freqs[0])

### Convert Bands to array
Bands = np.array(Bands)
print(Bands)
print(np.shape(Bands))

### The data to be printed 
datasave = np.column_stack((k_array,delta_array,Bands))
print(np.shape(datasave))

### Print the band structure 
namesavetxt = namesave+'-Band.txt'
with open(namesavetxt,'w') as file:
    np.savetxt(namesavetxt,datasave,'%.8f')