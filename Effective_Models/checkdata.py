import numpy as np 
import scipy 

epsilon = 1e-3

datapy = np.loadtxt('output_python12')
datam = np.loadtxt('output_octave12')

#print(i for i in range(len(datapy)) if abs(datapy[i]-datam[i])<epsilon)

for i in range(len(datapy)):
    if (abs(datapy[i] - datam[i]) > epsilon):
        print(i)

for i in range(len(datam)):
    if (abs(datam[i]) < epsilon):
        print(str(i)+','+str(datam[i]))