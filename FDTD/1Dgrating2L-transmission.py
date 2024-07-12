import numpy as np 
import scipy 
import meep as mp 

##### The function to calculate the transmitted flux 
def trans_func(shift,fcen,df,nfreq):
    ### The geometrical configuration 
    N = 3   # Number of gratings per side
    sy = 10 # Size of computing cell along y 

    shift = 0.3
    d = shift/2 # Layer 2 (upper) is shifted by d and layer 1 (lower) by -d 
    comp = mp.Ez 

    resolution = 16 # pixels/micrometer 

    InP = mp.Medium(index=3.17)
    Sb3S2 = mp.Medium(index=2.73)

    ### Widths of the layers
    w0 = 0.3
    w1 = 0.45
    w2 = 0.5

    ### Filling fractions
    f0 = 0.8
    f1 = 0.85
    f2 = 0.75 

    ### Interlayer distance
    dist = 0.1 

    ### Position of the light source away from the absorber 
    pos = 1.8 

    ### PML thickness 
    dpml = 2 

    ### Padding block 
    pad = 2 

    ### Size of cell in x direction 
    sx = 2 *(dpml + N + pad) 

    ### Define the cell 
    cell = mp.Vector3(sx,sy,0)

    ### Define the geometry 
    geometry = []

    for i in range(N+2):
        geometry.append(mp.Block(
            center = mp.Vector3(-0.5-i+d,(w0+dist)/2,0),
            size = mp.Vector3(f0,w0,mp.inf),
            material = InP))
        geometry.append(mp.Block(
            center = mp.Vector3(0.5+i+d,(w0+dist)/2,0),
            size = mp.Vector3(f0,w0,mp.inf),
            material = InP))
        geometry.append(mp.Block(
            center = mp.Vector3(-0.5-i+d,-(w1+dist)/2,0),
            size = mp.Vector3(f1,w1,mp.inf),
            material = Sb3S2))
        geometry.append(mp.Block(
            center = mp.Vector3(0.5+i+d,-(w2+dist)/2,0),
            size = mp.Vector3(f2,w2,mp.inf),
            material = Sb3S2))
        
    # Padding blocks
    geometry.append(mp.Block(
        size=mp.Vector3(pad+dpml,w0,mp.inf),
        center = mp.Vector3((pad+dpml-sx)/2,(w0+dist)/2,0),
        material = InP))
    geometry.append(mp.Block(
        size=mp.Vector3(pad+dpml,w1,mp.inf),
        center = mp.Vector3((pad+dpml-sx)/2,-(w1+dist)/2,0),
        material = Sb3S2))
    geometry.append(mp.Block(
        size=mp.Vector3(pad+dpml,w0,mp.inf),
        center = mp.Vector3((sx-pad-dpml)/2,(w0+dist)/2,0),
        material = InP))
    geometry.append(mp.Block(
        size=mp.Vector3(pad+dpml,w2,mp.inf),
        center = mp.Vector3((sx-pad-dpml)/2,-(w2+dist)/2,0),
        material = Sb3S2))
    
    sim = mp.Simulation(cell_size = cell,
                        geometry = geometry,
                        sources = [],
                        boundary_layers = [mp.PML(dpml,mp.X),mp.PML(1,mp.Y)],
                        resolution=resolution)    

    sim.sources.append(mp.Source(
        mp.GaussianSource(fcen,fwidth=df),
        component = comp,
        center = mp.Vector3(-0.5*sx+dpml+pos,(w0-w1)/2,0),
        size = mp.Vector3(0,w0+w1+dist)
    ))

    freg = mp.FluxRegion(center = mp.Vector3(0.5*sx-dpml-pos,(w0-w2)/2,0),
                         size=mp.Vector3(0,w0+w2+dist))

    ### Transmitted flux 
    trans = sim.add_flux(fcen,df,nfreq,freg)

    sim.run(until_after_sources=mp.stop_when_fields_decayed(1000,
                                                            comp,
                                                            mp.Vector3(0.5*sx-dpml-pad/2),
                                                            1e-4))
    
    return np.array(mp.get_fluxes(trans))   # List of fluxes 

def norm_func(fcen,df,nfreq):
    sy = 10 # Size of computing cell along y 

    comp = mp.Ez 

    resolution = 32 # pixels/micrometers 

    InP = mp.Medium(index = 3.17)
    Sb3S2 = mp.Medium(index = 2.73)

    w0 = 0.3 
    w1 = 0.45       # layer thickness 
    w2 = 0.5        
    dist = 0.1      # distance between layers 

    pos = 1.8       # position of light sources away from the absorber 
    dpml = 2        # PML thickness 
    pad = 2         # Padding block

    sx = 2*(dpml+pad)   # size of cell in x direction 

    cell = mp.Vector3(sx,sy,0)

    geometry = []

    geometry.append(mp.Block(
        size = mp.Vector3(pad+dpml,w0,mp.inf),
        center = mp.Vector3((pad+dpml-sx)/2,(w0+dist)/2,0),
        material = InP 
    ))
    geometry.append(mp.Block(
        size = mp.Vector3(pad+dpml,w1,mp.inf),
        center = mp.Vector3((pad+dpml-sx)/2,-(w1+dist)/2,0),
        material = Sb3S2
    ))
    geometry.append(mp.Block(
        size = mp.Vector3(pad+dpml,w0,mp.inf),
        center = mp.Vector3((sx-pad-dpml)/2,(w0+dist)/2,0),
        material = InP
    ))
    geometry.append(mp.Block(
        size = mp.Vector3(pad+dpml,w2,mp.inf),
        center = mp.Vector3((sx-pad-dpml)/2,-(w2+dist)/2,0),
        material = Sb3S2
    ))

    sim = mp.Simulation(cell_size = cell,
                        geometry = geometry,
                        sources = [],
                        boundary_layers = [mp.PML(dpml,mp.X),mp.PML(1,mp.Y)],
                        resolution = resolution)
    
    sim.sources.append(mp.Source(mp.GaussianSource(fcen,fwidth=df),
                                 component = comp,
                                 center = mp.Vector3(-0.5*sx+dpml+pos,(w0-w1)/2,0),
                                 size = mp.Vector3(0,w0+w1+dist)))
    
    freg = mp.FluxRegion(center = mp.Vector3(0.5*sx-dpml-pos,(w0-w2)/2,0),
                         size = mp.Vector3(0,w0+w2+dist))
    
    # Transmitted flux 
    trans = sim.add_flux(fcen,df,nfreq,freg)

    sim.run(until_after_sources = mp.stop_when_fields_decayed(
        1000,comp,mp.Vector3(0.5*sx-dpml-pad/2),1e-4))
    
    # Return the list of flux 
    return np.array(mp.get_fluxes(trans)), np.array(mp.get_flux_freqs(trans))

##### Define the MAIN function 
def main():
    ### The parameter of the source 
    fcen = 0.2387 
    df = 0.01 
    nfreq = 500 


    shift_array = np.linspace(0.4,0.6,3)

    t2, w2 = norm_func(fcen,df,nfreq)

    with open('data.txt','w') as file:
        for shift in shift_array:
            T = trans_func(shift,fcen,df,nfreq)/t2 
            np.savetxt(file,np.vstack((shift*np.ones(nfreq),w2,T)).transpose(),
                       fmt='%.6f')


##### Run the MAIN function 
if __name__ == '__main__':
    main()