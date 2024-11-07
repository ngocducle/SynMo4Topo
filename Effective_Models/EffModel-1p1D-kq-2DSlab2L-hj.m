%%% Solve for the edge state of the analytical model 
%%% of 2D slab bilayer 
%%% OCTAVE version 

%%%%% ================================================================================
%%%%% FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer with 
%%%%% kx = ky = k, qx = qy = q 
%%%%% The synthetic momentum: q = delta/(sqrt(2)*a) - 0.5
%%%%% that means q = 0 at shift = sqrt(2)*a/2 
%%%%% The model is derived directly from the Maxwell equation 
function H = Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,omega2,eta2,v2,U2,W2,alpha2,V,beta,dist,d0)
    H = zeros(8,8);
    K = 2*pi;

    %%% Block (1,1)
    H(1,1) = omega1 + eta1 + v1*k + v1*k*k/sqrt(2);
    H(1,2) = W1;
    H(1,3) = W1;
    H(1,4) = U1*(1+alpha1); 

    H(2,1) = W1; 
    H(2,2) = omega1 - eta1 + v1*k*k/sqrt(2);
    H(2,3) = U1*(1-alpha1);
    H(2,4) = W1;

    H(3,1) = W1;
    H(3,2) = U1*(1-alpha1);
    H(3,3) = omega1 - eta1 + v1*k*k/sqrt(2);
    H(3,4) = W1; 

    H(4,1) = U1*(1+alpha1);
    H(4,2) = W1;
    H(4,3) = W1;
    H(4,4) = omega1 + eta1 - v1*k + v1*k*k/sqrt(2);

    %%% Block (1,2)
    H(1,5) = -(V + beta*k + beta*k*k/sqrt(2))*exp(-i*K*q-dist/d0);
    H(2,6) = (V + beta*k*k/sqrt(2))*exp(-dist/d0);
    H(3,7) = (V + beta*k*k/sqrt(2))*exp(-dist/d0);
    H(4,8) = -(V - beta*k + beta*k*k/sqrt(2))*exp(i*K*q-dist/d0);

    %%% Block (2,1)
    H(5,1) = -(V + beta*k + beta*k*k/sqrt(2))*exp(i*K*q-dist/d0);
    H(6,2) = (V + beta*k*k/sqrt(2))*exp(-dist/d0);
    H(7,3) = (V + beta*k*k/sqrt(2))*exp(-dist/d0);
    H(8,4) = -(V - beta*k + beta*k*k/sqrt(2))*exp(-i*K*q-dist/d0);

    %%% Block (2,2)
    H(5,5) = omega2 + eta2 + v2*k + v2*k*k/sqrt(2);
    H(5,6) = W2;
    H(5,7) = W2; 
    H(5,8) = U2*(1+alpha2);

    H(6,5) = W2; 
    H(6,6) = omega2 - eta2 + v2*k*k/sqrt(2);
    H(6,7) = U2*(1-alpha2);
    H(6,8) = W2; 

    H(7,5) = W2;
    H(7,6) = U2*(1-alpha2);
    H(7,7) = omega2 - eta2 + v2*k*k/sqrt(2);
    H(7,8) = W2;

    H(8,5) = U2*(1+alpha2);
    H(8,6) = W2; 
    H(8,7) = W2;
    H(8,8) = omega2 + eta2 - v2*k + v2*k*k/sqrt(2);

end % function H 

%%%%% ================================================================================
%%%%% FUNCTION: the component Halpha 
function Halpha = H_alpha(q,omega1,eta1,W1,U1,alpha1,omega2,eta2,W2,U2,alpha2,V,dist,d0)
    Halpha = zeros(8,8);

    Vt = V*exp(-dist/d0);

    % Block (1,1)
    Halpha(1,1) = omega1 + eta1;
    Halpha(1,2) = W1;
    Halpha(1,3) = W1; 
    Halpha(1,4) = U1*(1+alpha1);

    Halpha(2,1) = W1;
    Halpha(2,2) = omega1 - eta1;
    Halpha(2,3) = U1*(1-alpha1);
    Halpha(2,4) = W1; 

    Halpha(3,1) = W1;
    Halpha(3,2) = U1*(1-alpha1);
    Halpha(3,3) = omega1 - eta1; 
    Halpha(3,4) = W1; 

    Halpha(4,1) = U1*(1+alpha1);
    Halpha(4,2) = W1;
    Halpha(4,3) = W1;
    Halpha(4,4) = omega1 + eta1; 

    % Block (1,2)
    Halpha(1,5) = -Vt*exp(-i*2*pi*q);
    Halpha(2,6) = Vt;
    Halpha(3,7) = Vt;
    Halpha(4,8) = -Vt*exp(i*2*pi*q);

    % Block (2,1)
    Halpha(5,1) = -Vt*exp(i*2*pi*q);
    Halpha(6,2) = Vt;
    Halpha(7,3) = Vt;
    Halpha(8,4) = -Vt*exp(-i*2*pi*q);

    % Block (2,2)
    Halpha(5,5) = omega2 + eta2;
    Halpha(5,6) = W2;
    Halpha(5,7) = W2; 
    Halpha(5,8) = U2*(1+alpha2);

    Halpha(6,5) = W2;
    Halpha(6,6) = omega2 - eta2;
    Halpha(6,7) = U2*(1-alpha2);
    Halpha(6,8) = W2; 

    Halpha(7,5) = W2;
    Halpha(7,6) = U2*(1-alpha2);
    Halpha(7,7) = omega2 - eta2; 
    Halpha(7,8) = W2; 

    Halpha(8,5) = U2*(1+alpha2);
    Halpha(8,6) = W2;
    Halpha(8,7) = W2;
    Halpha(8,8) = omega2 + eta2;
end % function Halpha 

%%%%% ================================================================================
%%%%% FUNCTION: Matrix H1 
function H1 = H_1(q,v1,v2,beta,dist,d0)
    H1 = zeros(8,8);

    B = beta*exp(-dist/d0);

    H1(1,1) = v1;
    H1(4,4) = -v1;
    H1(1,5) = -B*exp(-i*2*pi*q);
    H1(4,8) = B*exp(i*2*pi*q);
    H1(5,1) = -B*exp(i*2*pi*q);
    H1(8,4) = B*exp(-i*2*pi*q);
    H1(5,5) = v2;
    H1(8,8) = -v2;
end % function H1 

%%%%% ================================================================================
%%%%% FUNCTION: Matrix H2 
function H2 = H_2(q,v1,v2,beta,dist,d0)
    H2 = zeros(8,8);
    B = beta*exp(-dist/d0);

    H2(1,1) = v1;
    H2(2,2) = v1;
    H2(3,3) = v1;
    H2(4,4) = v1; 

    H2(1,5) = -B*exp(-i*2*pi*q);
    H2(2,6) = B;
    H2(3,7) = B;
    H2(4,8) = -B*exp(i*2*pi*q);

    H2(5,1) = -B*exp(i*2*pi*q);
    H2(6,2) = B;
    H2(7,3) = B;
    H2(8,4) = -B*exp(-i*2*pi*q);

    H2(5,5) = v2;
    H2(6,6) = v2;
    H2(7,7) = v2;
    H2(8,8) = v2;

    H2 = H2/sqrt(2);
end % function H2 

%%%%% ===============================================================================================
%%%%% FUNCTION: Hamiltonian of 2D photonic crystal slab bilayer with 
%%%%% kx = ky = k, qx = qy = q 
%%%%% The synthetic momentum: q = delta/(sqrt(2)*a) - 0.5
%%%%% that means q = 0 at shift = sqrt(2)*a/2 
%%%%% The model is derived directly from the Maxwell equation 
function H = Hamiltonian1(k,q,omega1,eta1,v1,U1,W1,alpha1,omega2,eta2,v2,U2,W2,alpha2,V,beta,dist,d0)
    Halpha = H_alpha(q,omega1,eta1,W1,U1,alpha1,omega2,eta2,W2,U2,alpha2,V,dist,d0);
    H1 = H_1(q,v1,v2,beta,dist,d0);
    H2 = H_2(q,v1,v2,beta,dist,d0);

    H = Halpha + H1*k + H2*k*k;

end % function Hamiltonian1

%%%%% ================================================================================
%%%%% FUNCTION: Solve the polynomial eigenvalue problem 
function [kvecs,kvals] = kPolyEig(E,q,omega1,eta1,v1,W1,U1,alpha1,omega2,eta2,v2,W2,U2,alpha2,V,beta,dist,d0)
    %%% The matrix Halpha 
    Halpha = H_alpha(q,omega1,eta1,W1,U1,alpha1,omega2,eta2,W2,U2,alpha2,V,dist,d0);
    H0 = Halpha - E*eye(8);

    %%% The matrix H1
    H1 = H_1(q,v1,v2,beta,dist,d0);

    %%% The matrix H2 
    H2 = H_2(q,v1,v2,beta,dist,d0);

    %%% Solve the polynomial eigenvalue problem 
    [kvecs,kvals] = polyeig(H0,H1,H2);

    %%% Move the arguments to the range 0 <= argument <= 2*pi 
    arg_array = zeros(8,1);

    for j = 1:8 
        arg_array(j) = arg(kvals(j));

        if (arg_array(j) < 0) 
            arg_array(j) = arg_array(j) + 2*pi; 
        end % IF 
    end % j-loop 

    %%% Rearrange kvecs and kvals in increasing order of argument
    [arg_list,ind] = sort(arg_array);
    kvals = kvals(ind);
    kvecs = kvecs(:,ind);
    kvecs = kvecs./norm(kvecs,'Fro','cols') % Frobenius norm summed over columns 
end % function kpolyEig 

%%%%% ================================================================================
%%% Parameters 
omega = 0.2978 
eta = -0.003
v = 0.317 
U = -0.01537
W = 0.001466
alpha = 0.05 

pomega = 0.015
omega1 = omega*(1+pomega)
omega2 = omega*(1-pomega)

v1 = v 
v2 = v 

pU = 0.1
U1 = U*(1+pU)
U2 = U*(1-pU)

pW = -0.1
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
Kmax = 0.05 
k_array = linspace(-Kmax,Kmax,Nk);

%%% The array of synthetic momenta 
Nq = 101 
Qmax = 0.3
q_array = linspace(-Qmax,Qmax,Nq);

%%% Criterion for 0 
epsilon = 1e-3 

%%% Number of E values to scan 
NE = 101

%%% Small incriment in band edge 
epsilonE = 1e-4 

%%% Initialize the edge states to be empty 
edge_state = [];

%%%%% ====================================================================================
%%%%% Scan over the q_array 
for iq = 1:Nq 
    %%%%% The synthetic momentum 
    q = q_array(iq);

    % The array of left and right bands 
    EL = zeros(Nk,8);
    ER = zeros(Nk,8);

    %%%%% =================================================================================
    for ik = 1:Nk 
        %%%%% The genuine momentum 
        k = k_array(ik);

        %%%%% =================================================================================
        %%%%% The left-hand side 
        %%% The Hamiltonian 
        H = Hamiltonian(k,q,omega1,eta1,v1,U1,W1,alpha1,omega2,eta2,v2,U2,W2,alpha2,V,beta,dist,d0);

        %%% Diagonalize the Hamiltonian 
        [Vec,D] = eig(H);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        %%% Save the energy eigenvalues to the array EL
        EL(ik,:) = eigval;
 
        %%%%% =================================================================================
        %%%%% The right-hand side 
        %%% The Hamiltonian 
        H = Hamiltonian(k,q,omega2,eta2,v2,U2,W2,alpha2,omega1,eta1,v1,U1,W1,alpha1,V,beta,dist,d0);

        %%% Diagonalize the Hamiltonian 
        [Vec,D] = eig(H);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        %%% Save the energy eigenvalues to the array ER 
        ER(ik,:) = eigval;
    end % ik-loop 

    %%%%% ==============================================================================
    %%% Calculate the obstructed bands 
    allEmax(iq) = max( min(EL(:,gap+1)), min(ER(:,gap+1)) );
    allEmin(iq) = min( max(EL(:,gap)), max(ER(:,gap)) );

    %%% Calculate the bulk bands 
    bulk1(iq) = min( min(EL(:,gap+1)), min(ER(:,gap+1)) );
    bulk2(iq) = max( max(EL(:,gap)), max(ER(:,gap)) );

    %%%%% ===============================================================================
    %%%%% Calculate the edge states at the synthetic momentum q 
    %%% Array of energy 
    E_array = linspace(bulk2(iq)+epsilonE,bulk1(iq)-epsilonE,NE);

    %%% Scan the energy array E_array 
    %for iE = 1:NE 
        % The energy 
    %    E = E_array(iE);

        % The matrix of eigenstates 
    %    W = zeros(8);

        % Left-hand Hamiltonian 
    %    printf("Left-hand Hamiltonian \n");
    %   [WL,kL] = kPolyEig(E,q,omega1,eta1,v1,W1,U1,alpha1,omega2,eta2,v2,W2,U2,alpha2,V,beta,dist,d0);

        % Right-hand Hamiltonian 
    %    printf("Right-hand Hamiltonian \n");
    %    [WR,kR] = kPolyEig(E,q,omega2,eta2,v2,W2,U2,alpha2,omega1,eta1,v1,W1,U1,alpha1,V,beta,dist,d0);

        %%% Combine WL ana WR to the matrix of eigenstates 
    %    W(:,1:4) = WL(:,5:8);
    %    W(:,5:8) = -WR(:,1:4);

        %%% The determinant of W 
    %    S = abs(det(W))

    %    if (S < epsilon) 
    %        edge_state = [edge_state;[q,E]];
    %    end % IF 

    %end % iE-loop 

end % iq-loop 

%%%%% ===================================================================================
%%%%% Plot the figure 
figure(1)
plot(q_array,allEmin,color='red'); hold on;
plot(q_array,allEmax,color='red');
plot(q_array,bulk1,color='red');
plot(q_array,bulk2,color='red');
%scatter(edge_state(:,1),edge_state(:,2));
hold off;
xlabel('q')
ylabel('E')
saveas(1,'transmission.png')

allEmax-bulk1 

allEmin-bulk2 