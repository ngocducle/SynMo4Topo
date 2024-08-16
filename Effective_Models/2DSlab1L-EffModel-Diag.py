import numpy as np 
import scipy 
import sympy as sp 
from sympy.matrices import Matrix

v,W,U,u = sp.symbols('v W U u')

H = Matrix([[0,W,W,U],[W,0,u,W],[W,u,0,W],[U,W,W,0]])

evalues = H.eigenvals()
evectors = H.eigenvects()

print('The energy eigenvalues are: ')
print(evalues)

print('The corresponding eigenvectors: ')
print(evectors)