%%% Solve for the edge state of the analytical model 
%%% of 2D slab bilayer 
%%% OCTAVE version 

%%% Parameters 
omega = 0.2978 
eta = -0.003 
v = 0.317 
U = -0.01537
W = 0.001466
alpha = 0.05 

pomega = 0.0 
omega1 = omega*(1+pomega)
omega2 = omega*(1-pomega)

v1 = v 
v2 = v 

pU = 0.05
U1 = U*(1+pU)
U2 = U*(1-pU)

pW = -0.05
W1 = W*(1+pW)
W2 = W*(1-pW)

eta1 = eta 
eta2 = eta 
alpha1 = alpha 
alpha2 = alpha 

d0 = 0.35 
dist = 0.1 
V = 0.038 
beta = -0.3 

%%% The gap we calculate the edge state 
gap = 1 

%%% The array of genuine momenta 
Nk = 101 
k_array = linspace(-0.25,0.25,Nk)

%%% The array of synthetic momenta 
Nq = 101 
q_array = linspace(-0.5,0.5,Nq)

%%% Criterion for 0 
epsilon = 1e-3 

%%% Number of E values to scan 
NE = 501

%%% Small incriment in band edge 
epsilonE = 1e-4 

