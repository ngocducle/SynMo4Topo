### Minh code for bilayer tight-binding 
### using pybinding 
### ATTENTION! The pybinding code actually does not work in my computers 
import numpy as np 
import pybinding as pb 
import matplotlib.pyplot as plt 
from math import pi,exp,sqrt 

D = 0.25 # interlayer separation 

def bilayer_square_lattice(a,d1,d2,alpha,shift,gate):
    # Hamiltonian around the AB stacking configuration 
    # a: lattice constant 
    # d1, d2: positions of 4 atoms within a unit cell; d1 = d2 preserves C4 
    # alpha: (around -0.1 to 0.1) deformation parameter 
    # shift: list or tuple 
    # We parametrize the hopping terms as follows 

    def interhop(xx,yy):
        # Interlayer hopping 
        c2 = D**2 / (xx**2 + yy**2 + D**2) 
        if sqrt(xx**2 + yy**2) < 4*a:
            res = exp(-(xx**2 + yy**2 + D**2)/a**2)*(2.5*c2-2)
        else:
            res = 0 
        return res 

    dx = shift[0]
    dy = shift[1]

    # Intralayer hopping 
    t1x = -exp(-d1**2/a**2)
    t1y = -exp(-d2**2/a**2)
    t2x = -exp(-(a-d1)**2/a**2)
    t2y = -exp(-(a-d2)**2/a**2)
    td1 = -0.5*(1-alpha)*exp( - (d1**2+d2**2)/a**2 )
    td2 = -0.5*(1+alpha)*exp( - (d1**2+d2**2)/a**2 )

    lat = pb.Lattice(a1 = [a,0], a2 = [0,a] ) 

    lat.add_sublattices( ('A1',[a/2+dx,a/2+dy,D/2],gate/2),
                         ('B1',[a/2+d1+dx,a/2+dy,D/2],gate/2),
                         ('C1',[a/2+dx,a/2+d2+dy,D/2],gate/2),
                         ('D1',[a/2+d1+dx,a/2+d2+dy,D/2],gate/2),
                         ('A2',[0,0,-D/2],-gate/2),
                         ('B2',[d1,0,-D/2],-gate/2),
                         ('C2',[0,d2,-D/2],-gate/2),
                         ('D2',[d1,d2,-D/2],-gate/2) )
                        
    lat.add_hopping( ( [0,0], 'A1','B1',t1x ),
                     ( [0,0], 'A1','C1',t1y ),
                     ( [0,0], 'B1','D1',t1y ),
                     ( [0,0], 'C1','D1',t1x ),
                     ( [0,0], 'A1','D1',td1 ),
                     ( [0,0], 'B1','C1',td2 ),
                     ( [1,0], 'B1','A1',t2x ),
                     ( [1,0], 'D1','C1',t2x ),
                     ( [0,1], 'C1','A1',t2y ),
                     ( [0,1], 'D1','B1',t2y ) )

    
    lat.add_hopping( ( [0,0], 'A2','B2',t1x ),
                     ( [0,0], 'A2','C2',t1y ),
                     ( [0,0], 'B2','D2',t1y ),
                     ( [0,0], 'C2','D2',t1x ),
                     ( [0,0], 'A2','D2',td1 ),
                     ( [0,0], 'B2','C2',td2 ),
                     ( [1,0], 'D2','C2',t2x ),
                     ( [1,0], 'D2','C2',t2x ),
                     ( [0,1], 'C2','A2',t2y ),
                     ( [0,1], 'D2','B2',t2y ) )

    # Interlayer hopping 
    atom2 = ['A2','B2','C2','D2']
    pos2 = np.array([[0,0],[d1,0],[0,d2],[d1,d2]])

    atom1 = ['A1','B1','C1','D1']
    pos1 = pos2 + a/2 + np.array([dx,dy])

    for nx in range(-3,4):
        for ny in range(-3,4):
            intercell = [nx,ny]
            for n1 in range(4):
                for n2 in range(4):
                    disp = pos1[n1] - np.array([nx*a,0]) - pos2[n2] - np.array([0,ny*a]) 
                    lat.add_one_hopping(intercell,atoms1[n1],atoms2[n2],interhop(disp[0],disp[1]))

    return lat
