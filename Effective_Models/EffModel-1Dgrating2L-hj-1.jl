import Pkg 
Pkg.add("Plots")
Pkg.add("MAT")
using LinearAlgebra;
using Plots;

alpha = 0.0513 
U = 0.0208/alpha 
V = 2.73*U/alpha 
Delta = 0.3*U 

U = 1.0 
Del = -0.1 

V = 1.1 

i = im 
Nq = 1000 
qq = range(-3.0,stop=3.0,length=Nq)

Nph = 100 
pph = range(pi/2.0-0.5,stop=pi/2.0+0.5,length=Nph)

allEdgeStates = Array{Array{Float64,1},1}(UndefInitializer(),Nph);
allEmin = Array{Float64,1}(UndefInitializer(),Nph);
allEmax = Array{Float64.1}(UndefInitializer(),Nph);

