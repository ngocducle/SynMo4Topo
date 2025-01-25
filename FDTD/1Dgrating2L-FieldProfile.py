import numpy as np
import scipy 
import meep as mp 
import matplotlib.pyplot as plt 

### Parameters 
sy = 18 # Size of computational cell along y
shift = 0.5
d = shift/2 # Layer 2 (upper) is shifted by d and layer 1 (lower) by -d 
comp = mp.Ez    # The field component

resolution =16 # pixels/micrometer 

InP = mp.Medium(index = 3.17)
Sb3S2 = mp.Medium(index = 2.73)

w0 = 0.3
w1 = 0.46   # width of layer 
w2 = 0.5
f0 = 0.8
f1 = 0.81
f2 = 0.77
dist = 0.1  # distance between layers

dpml = 2    # PML thickness
dabs = 6

N = 10  # Number of unit cell per side

sx = 2*(dabs+N) # size of cell in x direction 
print(sx)

cell = mp.Vector3(sx,sy,0)

geometry = []

for i in range(N+4):
    geometry.append(mp.Block(center = mp.Vector3(-0.5-i+d,(w0+dist)/2,0),
                             size = mp.Vector3(f0,w0,mp.inf),
                             material = InP))
    geometry.append(mp.Block(center = mp.Vector3(0.5+i+d,(w0+dist)/2,0),
                             size = mp.Vector3(f0,w0,mp.inf),
                             material = InP))
    geometry.append(mp.Block(center = mp.Vector3(-0.5-i-d,-(w1+dist)/2,0),
                             size = mp.Vector3(f1,w1,mp.inf),
                             material = Sb3S2))
    geometry.append(mp.Block(center = mp.Vector3(0.5+i+d,-(w2+dist)/2,0),
                             size = mp.Vector3(f2,w2,mp.inf),
                             material = Sb3S2))
    
sim = mp.Simulation(
    cell_size = cell,
    geometry = geometry,
    sources = [],
    boundary_layers = [mp.Absorber(dabs,mp.X),
                       mp.PML(dpml,mp.Y)],
    resolution = resolution
)

fcen = 0.2383
df = 0.0016

sim.sources.append(
    mp.Source(
        mp.GaussianSource(fcen,fwidth=df),
        component = comp,
        center = mp.Vector3(-d-0.23,-0.2)
    )
)

nonpml_vol = mp.Volume(
    center = mp.Vector3(),
    size = mp.Vector3(sx - 2*dabs,sy - 2*dpml)
)

#print(sx-2*dabs)
x_array = np.linspace(-N,N,2*N*resolution+2)

dft_obj = sim.add_dft_fields([comp],fcen,0,1,where=nonpml_vol)

sim.run(until_after_sources=mp.stop_when_fields_decayed(300,comp,mp.Vector3(-d-0.2,0.18),1e-7))

eps_data = sim.get_array(vol=nonpml_vol,component=mp.Dielectric)
ez_data = np.real(sim.get_dft_array(dft_obj,comp,0))

np.savetxt('dielectric.txt',eps_data)
np.savetxt('full_field.txt',ez_data)
np.savetxt('dft_field_along_x.txt',
           np.column_stack((x_array,np.sum(np.square(ez_data),axis=1)))
           )