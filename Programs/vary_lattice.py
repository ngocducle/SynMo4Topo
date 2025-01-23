import numpy as np
import meep as mp
from meep import mpb 
import matplotlib.pyplot as plt


##### The conclusion of this code is:
# The axes of the Block object and the vertices of the Prism object
# are defined in meep.Vector3, in terms of the basis lattice

num_bands = 8 

alpha = np.pi/4

geometry_lattice = mp.Lattice(
    size = mp.Vector3(1,1),
    basis1 = mp.Vector3(np.cos(alpha),np.sin(alpha)),
    basis2 = mp.Vector3(np.cos(alpha),-np.sin(alpha))
)

#geometry = [
#    mp.Cylinder(
#        radius = 0.2,
#        material = mp.Medium(epsilon=12)
#    )
#]

edge = 0.5

#geometry = [
#    mp.Block(
#        center = mp.Vector3(0.0,0.0,0.0),
#        size = mp.Vector3(edge,edge),
#        material = mp.Medium(epsilon=12)
#    )
#]

geometry = [
    mp.Prism(
        vertices = [mp.Vector3(-0.5,0.5,0),mp.Vector3(0.0,0.5,0),mp.Vector3(0.5,0,0)],
        height = 1,
        axis = mp.Vector3(0,0,1),
        center = mp.Vector3(0,0,0),
        material = mp.Medium(epsilon=12)
    )
]

k_points = [
    mp.Vector3(0,0,0),
    mp.Vector3(0.5,0,0)
]

resolution = 32

ms = mpb.ModeSolver(
    geometry = geometry,
    geometry_lattice = geometry_lattice,
    k_points = k_points,
    resolution = resolution,
    num_bands = num_bands
)

ms.run_te()

md = mpb.MPBData(rectify=True, periods=3, resolution=32)
eps = ms.get_epsilon()
converted_eps = md.convert(eps)

print(np.cos(alpha))
print(np.sin(alpha))

plt.imshow(converted_eps.T, interpolation='spline36', cmap='binary')
plt.axis('off')
plt.show()