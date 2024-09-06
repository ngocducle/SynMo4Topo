import numpy as np 
import scipy 
import sympy as sp 
from sympy.matrices import Matrix 

U,V = sp.symbols('U V') 

H = Matrix([[0,U,-1j*V,0],[U,0,0,1j*V],[1j*V,0,0,U],[0,-1j*V,U,0]]) 

evalues = H.eigenvals() 

evectors = H.eigenvects() 

print('The energy eigenvalues are: ') 
print(evalues)  

print('The corresponding eigenvectors: ') 
print(evectors)

### The matrix form by columns of eigenvectors of H 
UMatrix = Matrix(
    [ 
        [ 1j*(U**2+V**2)**0.5/V , 1j*U/V, -1j*(U**2+V**2)**0.5/V, 1j*U/V ], 
    [ -1j*U/V, -1j*(U**2+V**2)**0.5/V, -1j*U/V , 1j*(U**2+V**2)**0.5/V ], 
     [ 1,0,1,0  ],
     [ 0,1,0,1 ] 
]
)

print('The U matrix')
print(UMatrix) 

### The inverse of UMatrix 
Uinv = UMatrix.inv() 
print('Uinv = ')
print(Uinv) 
