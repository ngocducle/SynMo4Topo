%%% Solve for the dispersion of the Hamiltonian 
%%% of 2D slab bilayer 
%%% OCTAVE version 

%%%%% ===================================================================================
%%%%% FUNCTION: Ha 
%%%%% The zero order of the Hamiltonian 
function Ha = H_a(q,omega1,eta1,v1,U1,W1,alpha1,omega2,eta2,v2,U2,W2,alpha2,V,dist,d0)
    Ha = zeros(8,8);
    Vt = V*exp(-dist/d0);

    %%% Block (1,1)
    Ha(1,1) = omega1 + eta1;
    Ha(1,2) = W1; 
    Ha(1,3) = W1; 
    Ha(1,4) = U1*(1+alpha1);

    Ha(2,1) = W1; 
    Ha(2,2) = omega1 - eta1; 
    Ha(2,3) = U1*(1-alpha1);
    Ha(2,4) = W1;

    Ha(3,1) = W1; 
    Ha(3,2) = U1*(1-alpha1);
    Ha(3,3) = omega1 - eta1; 
    Ha(3,4) = W1; 

    Ha(4,1) = U1*(1+alpha1);
    Ha(4,2) = W1; 
    Ha(4,3) = W1; 
    Ha(4,4) = omega1 + eta1; 

    %%% Block (1,2)
    Ha(1,5) = -Vt*exp(-i*2*pi*q);
    Ha(2,6) = Vt;
    Ha(3,7) = Vt;
    Ha(4,8) = -Vt*exp(i*2*pi*q);

    %%% Block (2,1)
    Ha(5,1) = -Vt*exp(i*2*pi*q);
    Ha(6,2) = Vt; 
    Ha(7,3) = Vt;
    Ha(8,4) = -Vt*exp(-i*2*pi*q);

    %%% Block (2,2)
    Ha(5,5) = omega2 + eta2;
    Ha(5,6) = W2; 
    Ha(5,7) = W2; 
    Ha(5,8) = U2*(1+alpha2);

    Ha(6,5) = W2; 
    Ha(6,6) = omega2 - eta2; 
    Ha(6,7) = U2*(1-alpha2);
    Ha(6,8) = W2;

    Ha(7,5) = W2; 
    Ha(7,6) = U2*(1-alpha2);
    Ha(7,7) = omega2 - eta2; 
    Ha(7,8) = W2; 

    Ha(8,5) = U2*(1+alpha2);
    Ha(8,6) = W2; 
    Ha(8,7) = W2; 
    Ha(8,8) = omega2 + eta2;
end % function Ha 

%%%%% ===================================================================================
%%%%% FUNCTION: The coefficient of k
function H1 = H_1(q,v1,v2,beta,dist,d0)
    H1 = zeros(8,8);
    K = 2*pi; 
    
    H1(1,1) = v1;
    H1(4,4) = -v1;

    H1(1,5) = -beta*exp(-dist/d0-i*K*q);
    H1(4,8) = beta*exp(-dist/d0+i*K*q);

    H1(5,1) = -beta*exp(-dist/d0+i*K*q);
    H1(8,4) = beta*exp(-dist/d0-i*K*q);

    H1(5,5) = v2;
    H1(8,8) = -v2; 

end % function H1

%%%%% ===================================================================================
%%%%% FUNCTION: The coefficient of k^2
function H2 = H_2(q,v1,v2,beta,dist,d0)
    H2 = zeros(8,8);

    B = beta*exp(-dist/d0)/sqrt(2);

    K = 2*pi; 

    H2(1,1) = v1/sqrt(2);
    H2(2,2) = v1/sqrt(2);
    H2(3,3) = v1/sqrt(2);
    H2(4,4) = v1/sqrt(2);

    H2(1,5) = -beta*exp(-i*K*q-dist/d0)/sqrt(2);
    H2(2,6) = beta*exp(-dist/d0)/sqrt(2); 
    H2(3,7) = beta*exp(-dist/d0)/sqrt(2);
    H2(4,8) = -beta*exp(i*K*q-dist/d0)/sqrt(2);

    H2(5,1) = -beta*exp(i*K*q-dist/d0)/sqrt(2);
    H2(6,2) = beta*exp(-dist/d0)/sqrt(2); 
    H2(7,3) = beta*exp(-dist/d0)/sqrt(2);
    H2(8,4) = -beta*exp(-i*K*q-dist/d0)/sqrt(2);

    H2(5,5) = v2/sqrt(2);
    H2(6,6) = v2/sqrt(2);
    H2(7,7) = v2/sqrt(2);
    H2(8,8) = v2/sqrt(2);
end % function H2 

%%%%% ==============================================================================
%%%%% FUNCTION: Hamiltonian 
function H = Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,omega2,eta2,v2,U2,W2,alpha2,V,beta,dist,d0)

    Ha = H_a(q,omega1,eta1,v1,U1,W1,alpha1,omega2,eta2,v2,U2,W2,alpha2,V,dist,d0);
    H1 = H_1(q,v1,v2,beta,dist,d0);
    H2 = H_2(q,v1,v2,beta,dist,d0);

    H = Ha + H1*k + H2*k*k;

end % function Hamiltonian

%%%%% ====================================================================================
%%%%% Parameters 
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

pU = 0.0
U1 = U*(1+pU)
U2 = U*(1-pU)

pW = 0.0 
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

%%%%% ===================================================================================
%%%%% The arrays of synthetic and genuine momenta 
%%% The array of genuine momenta 
Nk = 100
Nq = 101 

Qmax = 0.5
Kmax = 0.05 

k_array = [linspace(0,0,Nq),linspace(Kmax/Nq,Kmax,Nk)]
q_array = [linspace(-Qmax,0,Nq),linspace(0,0,Nk)]

%%%%% ===================================================================================
%%%%% The array of energy 
E_array = zeros(Nq+Nk,8)

%%%%% ===================================================================================
%%%%% Scan the hybrid momenta
for i = 1:(Nq+Nk)
    %%% The genuine momentum 
    k = k_array(i)

    %%% The synthetic momentum 
    q = q_array(i)

    %%% The Hamiltonian 
    H = Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,omega2,eta2,v2,U2,W2,alpha2,V,beta,dist,d0);

    %%% Diagonalize the Hamiltonian 
    [Vec,D] = eig(H);
    [eigval,ind] = sort(diag(D));
    eigval = real(eigval);
    Ds = D(ind,ind);
    eigvec = Vec(:,ind);

    %%% Save the energy eigenvalues to the array E_array 
    E_array(i,:) = eigval; 

end % i-loop 

%%%%% ====================================================================================
%%%%% Plot the figure 
figure(1)
plot(E_array(:,1:3))
ylim([0.24,0.27])
saveas(1,'BandStructure.png')