import numpy as np
import scipy
import scipy.linalg as sla
import matplotlib.pyplot as plt

import sys 
sys.path.insert(0,'../src/')
from EffectiveModels import EffModel_2DSlab2L_M

### ================================================================================ ###
#                                                                                      #
#           Fit the effective model for 2D photonic crystal slab bilayer               #
#                                                                                      #
### ================================================================================ ###

### ====================================================================================
### Fit the modes arising from monolayer EVEN bands

# The shift along the x-direction 
deltax = 0.5

# The shift along the y-direction
deltay = 0.5

# Frequency at the M point 
omegaMe = 0.37773038
print('omegaMe = '+str(omegaMe))

# Group velocity  
ve = 0.35
print('Group velocity = '+str(ve))

# Counter-propagating coupling strength 
Ue = -0.02365418
print('Ue = '+str(Ue))

# Orthogonally-propagating coupling stength 
We = 0.00721296 
print('We = '+str(We))

# Interlayer coupling strength at zero interlayer distance 
V0e = 0.044 
print('V0e = '+str(V0e))

# Characteristic length for interlayer coupling 
d0e = 0.35 
print('d0e = '+str(d0e))

# Diagonalize the effective Hamiltonian 
Ndist = 121
DistArray = np.linspace(0.0,3.0,Ndist)

Eeven1 = np.empty(Ndist, dtype = complex) 
Eeven2 = np.empty(Ndist, dtype = complex)
Eeven3 = np.empty(Ndist, dtype = complex)
Eeven4 = np.empty(Ndist, dtype = complex)  
Eeven5 = np.empty(Ndist, dtype = complex)
Eeven6 = np.empty(Ndist, dtype = complex)
Eeven7 = np.empty(Ndist, dtype = complex)
Eeven8 = np.empty(Ndist, dtype = complex) 

# Calculate the energy bands arising from monolayer EVEN modes as function of the 
# distance 
for i in range(Ndist):
    qx = 0.0 
    qy = 0.0 

    d = DistArray[i]

    Ve = V0e*np.exp(-d/d0e)

    evalues, evectors = EffModel_2DSlab2L_M(qx,qy,deltax,deltay,omegaMe,ve,Ue,We,Ve)

    Eeven1[i] = evalues[0]
    Eeven2[i] = evalues[1]
    Eeven3[i] = evalues[2]
    Eeven4[i] = evalues[3] 
    Eeven5[i] = evalues[4]
    Eeven6[i] = evalues[5]
    Eeven7[i] = evalues[6]
    Eeven8[i] = evalues[7] 

### ============================================================================
### Load the data 
#DataEven = np.loadtxt('2DSlab2L-CHole-DistScanMzeven.txt')
#DataOdd = np.loadtxt('2DSlab2L-CHole-DistScanM-zodd.txt') 

#DistArrayMPB = DataEven[:,0]

Data = np.loadtxt('2DSlab2L-CHole-DistScanM-all.txt')
DistArrayMPB = Data[:,0]

### ==================================================================================
### Plot the figure for EVEN bands 
'''''
fig, ax = plt.subplots()
plt.plot(DistArrayMPB, DataEven[:,1],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1)
plt.plot(DistArrayMPB, DataEven[:,2],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1) 
plt.plot(DistArrayMPB, DataEven[:,3],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1)  
plt.plot(DistArrayMPB, DataEven[:,4],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1)  
plt.plot(DistArrayMPB, DataOdd[:,1],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1)
plt.plot(DistArrayMPB, DataOdd[:,2],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1) 
plt.plot(DistArrayMPB, DataOdd[:,3],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1)  
plt.plot(DistArrayMPB, DataOdd[:,4],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1) 


plt.plot(DistArray,Eeven1,linewidth=1,color='blue')
plt.plot(DistArray,Eeven2,linewidth=1,color='blue')
plt.plot(DistArray,Eeven3,linewidth=1,color='blue')
plt.plot(DistArray,Eeven4,linewidth=1,color='blue')
plt.plot(DistArray,Eeven5,linewidth=1,color='blue')
plt.plot(DistArray,Eeven6,linewidth=1,color='blue')
plt.plot(DistArray,Eeven7,linewidth=1,color='blue')
plt.plot(DistArray,Eeven8,linewidth=1,color='blue')

plt.xlim(0.0,2.0)
ax.set_xlabel(r'$d / a$', fontsize = 14)
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('z-even bands')
plt.savefig('Even-bands.png')
plt.show() 
'''

### ==================================================================================
### Fit the ODD bands
    
# Frequency at the M point
omegaMo = 0.435940955
print('omegaMo = '+str(omegaMo)) 

# Group velocity
vo = 0.41
print('Group velocity = '+str(vo))

# Counter-propagating coupling strength 
Uo =  0.008515315
print('Uo = '+str(Uo))

# Orthogonally-propagating coupling strength
Wo = 0.0154966775
print('Wo = '+str(Wo)) 

# Interlayer coupling strength amplitude
V0o = 0.049
print('V0o = '+str(V0o))

# Characteristic length for interlayer coupling
d0o = 0.42
print('d0o = '+str(d0o)) 

# Diagonalize the effective Hamiltonian
Ndist = 121 
DistArray = np.linspace(0.0,3.0,Ndist) 

Eodd1 = np.empty(Ndist, dtype = complex) 
Eodd2 = np.empty(Ndist, dtype = complex)
Eodd3 = np.empty(Ndist, dtype = complex)
Eodd4 = np.empty(Ndist, dtype = complex)  
Eodd5 = np.empty(Ndist, dtype = complex)
Eodd6 = np.empty(Ndist, dtype = complex)
Eodd7 = np.empty(Ndist, dtype = complex)
Eodd8 = np.empty(Ndist, dtype = complex) 

# Calculate the energy bands of odd modes as function of the distance 
for i in range(Ndist):
    qx = 0.0
    qy = 0.0 

    d = DistArray[i]

    Vo = V0o * np.exp(-d/d0o)

    evalues, evectors = EffModel_2DSlab2L_M(qx,qy,0.0,0.0,omegaMo,vo,Uo,Wo,Vo) 
    
    Eodd1[i] = evalues[0] 
    Eodd2[i] = evalues[1]
    Eodd3[i] = evalues[2]
    Eodd4[i] = evalues[3] 
    Eodd5[i] = evalues[4]
    Eodd6[i] = evalues[5]
    Eodd7[i] = evalues[6]
    Eodd8[i] = evalues[7] 

### ==================================================================================
### Plot the figure for ODD bands 
''''
fig, ax = plt.subplots()
plt.plot(DistArrayMPB, DataEven[:,5],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1)
plt.plot(DistArrayMPB, DataEven[:,6],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1) 
plt.plot(DistArrayMPB, DataEven[:,7],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1)  
plt.plot(DistArrayMPB, DataEven[:,8],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1)  
plt.plot(DistArrayMPB, DataOdd[:,5],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1)
plt.plot(DistArrayMPB, DataOdd[:,6],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1) 
plt.plot(DistArrayMPB, DataOdd[:,7],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1)  
plt.plot(DistArrayMPB, DataOdd[:,8],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1) 

plt.plot(DistArray,Eodd1,linewidth=1,color='red')
plt.plot(DistArray,Eodd2,linewidth=1,color='red')
plt.plot(DistArray,Eodd3,linewidth=1,color='red')
plt.plot(DistArray,Eodd4,linewidth=1,color='red')
plt.plot(DistArray,Eodd5,linewidth=1,color='red')
plt.plot(DistArray,Eodd6,linewidth=1,color='red')
plt.plot(DistArray,Eodd7,linewidth=1,color='red')
plt.plot(DistArray,Eodd8,linewidth=1,color='red')

plt.xlim(0.0,2.0)
ax.set_xlabel(r'$d / a$', fontsize = 14)
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('z-odd bands')
plt.savefig('Odd-bands.png')
plt.show() 
'''''
### ==================================================================
### Plot all the bands together
'''
fig, ax = plt.subplots()
plt.plot(DistArrayMPB, DataEven[:,1],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1)
plt.plot(DistArrayMPB, DataEven[:,2],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1) 
plt.plot(DistArrayMPB, DataEven[:,3],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1)  
plt.plot(DistArrayMPB, DataEven[:,4],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1)  
plt.plot(DistArrayMPB, DataOdd[:,1],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1)
plt.plot(DistArrayMPB, DataOdd[:,2],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1) 
plt.plot(DistArrayMPB, DataOdd[:,3],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1)  
plt.plot(DistArrayMPB, DataOdd[:,4],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1) 


#plt.plot(DistArray,Eeven1,linewidth=1,color='blue')
#plt.plot(DistArray,Eeven2,linewidth=1,color='blue')
#plt.plot(DistArray,Eeven3,linewidth=1,color='blue')
#plt.plot(DistArray,Eeven4,linewidth=1,color='blue')
#plt.plot(DistArray,Eeven5,linewidth=1,color='blue')
#plt.plot(DistArray,Eeven6,linewidth=1,color='blue')
#plt.plot(DistArray,Eeven7,linewidth=1,color='blue')
#plt.plot(DistArray,Eeven8,linewidth=1,color='blue')

plt.plot(DistArrayMPB, DataEven[:,5],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1)
plt.plot(DistArrayMPB, DataEven[:,6],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1) 
plt.plot(DistArrayMPB, DataEven[:,7],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1)  
plt.plot(DistArrayMPB, DataEven[:,8],'o',markerfacecolor='blue',\
    markeredgecolor='blue',markersize=1)  
plt.plot(DistArrayMPB, DataOdd[:,5],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1)
plt.plot(DistArrayMPB, DataOdd[:,6],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1) 
plt.plot(DistArrayMPB, DataOdd[:,7],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1)  
plt.plot(DistArrayMPB, DataOdd[:,8],'o',markerfacecolor='red',\
    markeredgecolor='red',markersize=1) 

#plt.plot(DistArray,Eodd1,linewidth=1,color='red')
#plt.plot(DistArray,Eodd2,linewidth=1,color='red')
#plt.plot(DistArray,Eodd3,linewidth=1,color='red')
#plt.plot(DistArray,Eodd4,linewidth=1,color='red')
#plt.plot(DistArray,Eodd5,linewidth=1,color='red')
#plt.plot(DistArray,Eodd6,linewidth=1,color='red')
#plt.plot(DistArray,Eodd7,linewidth=1,color='red')
#plt.plot(DistArray,Eodd8,linewidth=1,color='red')

plt.xlim(0.0,2.0)
ax.set_xlabel(r'$d / a$', fontsize = 14)
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('All bands')
plt.savefig('All-bands.png')
plt.show() 
''' 


fig, ax = plt.subplots()
plt.plot(DistArrayMPB, Data[:,1],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1)
plt.plot(DistArrayMPB, Data[:,2],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1) 
plt.plot(DistArrayMPB, Data[:,3],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1)  
plt.plot(DistArrayMPB, Data[:,4],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1)  
plt.plot(DistArrayMPB, Data[:,5],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1)
plt.plot(DistArrayMPB, Data[:,6],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1) 
plt.plot(DistArrayMPB, Data[:,7],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1)  
plt.plot(DistArrayMPB, Data[:,8],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1) 


plt.plot(DistArray,Eeven1,linewidth=1,color='blue')
plt.plot(DistArray,Eeven2,linewidth=1,color='blue')
plt.plot(DistArray,Eeven3,linewidth=1,color='blue')
plt.plot(DistArray,Eeven4,linewidth=1,color='blue')
plt.plot(DistArray,Eeven5,linewidth=1,color='blue')
plt.plot(DistArray,Eeven6,linewidth=1,color='blue')
plt.plot(DistArray,Eeven7,linewidth=1,color='blue')
plt.plot(DistArray,Eeven8,linewidth=1,color='blue')

plt.plot(DistArrayMPB, Data[:,9],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1)
plt.plot(DistArrayMPB, Data[:,10],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1) 
plt.plot(DistArrayMPB, Data[:,11],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1)  
plt.plot(DistArrayMPB, Data[:,12],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1)  
plt.plot(DistArrayMPB, Data[:,13],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1)
plt.plot(DistArrayMPB, Data[:,14],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1) 
plt.plot(DistArrayMPB, Data[:,15],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1)  
plt.plot(DistArrayMPB, Data[:,16],'o',markerfacecolor='green',\
    markeredgecolor='green',markersize=1) 

plt.plot(DistArray,Eodd1,linewidth=1,color='red')
plt.plot(DistArray,Eodd2,linewidth=1,color='red')
plt.plot(DistArray,Eodd3,linewidth=1,color='red')
plt.plot(DistArray,Eodd4,linewidth=1,color='red')
plt.plot(DistArray,Eodd5,linewidth=1,color='red')
plt.plot(DistArray,Eodd6,linewidth=1,color='red')
plt.plot(DistArray,Eodd7,linewidth=1,color='red')
plt.plot(DistArray,Eodd8,linewidth=1,color='red')

plt.xlim(0.0,2.0)
ax.set_xlabel(r'$d / a$', fontsize = 14)
ax.set_ylabel(r'$\omega a / (2 \pi c)$', fontsize = 14) 
plt.title('All bands')
plt.savefig('All-bands.png')
plt.show() 