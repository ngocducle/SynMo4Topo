import numpy as np 
import scipy 
import matplotlib.pyplot as plt 

data = np.loadtxt('q_0.3050.txt')

freq = data[:,0]
T = data[:,1]

fig,ax = plt.subplots()
ax.plot(freq,T)
plt.show() 
