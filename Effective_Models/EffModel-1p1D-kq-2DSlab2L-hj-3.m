%%%%% Solve for the edge states of the effective model 
%%%%% of 2D slab 2L 
%%%%% OCTAVE version 

%%%%% ===================================================================================
%%%%% FUNCTION: Ha 
%%%%% The zero order of the Hamiltonian 
function Ha = H_a(q,omega,domega,eta,U,dU,W,dW,alpha,V,dist,d0)
    Ha = zeros(8);

    omega1 = omega + domega; 
    eta1 = eta; 
    U1 = U+dU;
    W1 = W+dW; 
    alpha1 = alpha;

    omega2 = omega - domega; 
    eta2 = eta; 
    U2 = U-dU;
    W2 = W-dW; 
    alpha2 = alpha;

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

end % function H_a 

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
function H = Hamiltonian(k,q,omega,domega,eta,v1,v2,U,dU,W,dW,alpha,V,beta,dist,d0)

    Ha = H_a(q,omega,domega,eta,U,dU,W,dW,alpha,V,dist,d0);
    H1 = H_1(q,v1,v2,beta,dist,d0);
    H2 = H_2(q,v1,v2,beta,dist,d0);

    H = Ha + H1*k + H2*k*k;

end % function Hamiltonian

%%%%% ==================================================================================
%%%%% FUNCTION: Solve the polynomial eigenvalue problem 
function [kvecs,kvals] = kPolyEig(E,q,omega,domega,eta,v1,v2,U,dU,W,dW,alpha,V,beta,dist,d0)

    Ha = H_a(q,omega,domega,eta,U,dU,W,dW,alpha,V,dist,d0);
    H0 = Ha-E*eye(8);
    H1 = H_1(q,v1,v2,beta,dist,d0);
    H2 = H_2(q,v1,v2,beta,dist,d0);

    %H = Ha + H1*k + H2*k*k;

    [kvecs,kvals] = polyeig(H0,H1,H2);

    %%% Move the arguments to the range 0 <= argument <= 2*pi 
    arg_array = zeros(16,1);

    for j = 1:16
        arg_array(j) = arg(kvals(j));

        if (arg_array(j) < 0) 
            arg_array(j) = arg_array(j) + 2*pi;
        end % IF 
    end % j-loop 

    %%% Rearrange kvecs and kvals in increasing order of argument 
    [arg_list,ind] = sort(arg_array);
    kvals = kvals(ind);
    kvecs = kvecs(:,ind);
    kvecs = kvecs./norm(kvecs,'Fro','cols'); % Frobenius norm summed over columns
    
    ind; 
    kvals;
    kvecs;

end % function kPolyEig 

%%%%% ==================================================================================
%%%%% FUNCTION: Derivative of electric field 
function DE = FieldDerivative(k)
    DE = diag([sqrt(2)+2*k, ... 
               2*k, ... 
               2*k, ... 
               -sqrt(2)+2*k, ... 
               -sqrt(2)-2*k, ... 
               2*k, ... 
               2*k, ... 
               sqrt(2)-2*k]);

    %DE = diag([k,k,k,k,-k,k,k,-k]);
end % function FieldDerivative 

%%%%% ==================================================================================
%%%%% FUNCTION: 16x16 Coefficient matrix 
function WW = CoefficientMatrix(WL,kL,WR,kR)
    WW = zeros(16);

    WW(1:8,1:8) = WL(:,9:16);

    WW(1:8,9:16) = -WR(:,1:8);

    for j = 1:8 
        WW(9:16,j) = FieldDerivative(kL(j+8))*WL(:,j+8);
    end % j-loop 

    for j = 1:8 
        WW(9:16,j+8) = -FieldDerivative(kR(j))*WR(:,j);
    end % j-loop 

end % function CoefficientMatrix  

%%%%% ==================================================================================
%%%%% Parameters 
omega = 0.2978
domega = 0.015*omega 
eta = 0.003 
v = 0.317 
U = -0.01537 
dU = 0.1*U 
W = 0.001466 
dW = -0.1*W 
alpha = -0.05 

v1 = v 
v2 = v 

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
Nq = 201 
Qmax = 0.5
q_array = linspace(-Qmax,Qmax,Nq);

%%% Criterion for 0 
epsilon = 1e-3; 

%%% Number of E values to scan 
NE = 501 

%%% Small increment in band edge 
epsilonE = 1e-4

%%% Initialize the edge states to be empty 
edge_state = [];

%%%%% ===============================================================================
%%%%% Scan the q_array 
for iq = 1:Nq 
    %%%%% The synthetic momenta 
    q = q_array(iq);

    %%%%% The arrays of left and right bands 
    EL = zeros(Nk,8);
    ER = zeros(Nk,8);

    %%%%% ============================================================================
    for ik = 1:Nk 
        %%%%% The genuine momenta 
        k = k_array(ik);

        %%%%% ========================================================================
        %%%%% The left-hand side 
        %%% The Hamiltonian 
        H = Hamiltonian(k,q,omega,domega,eta,v1,v2,U,dU,W,dW,alpha,V,beta,dist,d0);

        %%% Diagonalize the Hamiltonian 
        [Vec,D] = eig(H);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        %%% Save the energy eigenvalues to the array EL 
        EL(ik,:) = eigval;

        %%%%% ========================================================================
        %%%%% The right-hand side 
        %%% The Hamiltonian 
        H = Hamiltonian(k,q,omega,-domega,eta,v1,v2,U,-dU,W,-dW,alpha,V,beta,dist,d0);

        %%% Diagonalize the Hamiltonian 
        [Vec,D] = eig(H);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        %%% Save the energy eigenvalues to the array EL 
        ER(ik,:) = eigval;

    end % ik-loop   

    %%%%% =============================================================================
    %%% Calculate the obstructed bands 
    allEmax(iq) = max( min(EL(:,gap+1)), min(ER(:,gap+1)) );
    allEmin(iq) = min( max(EL(:,gap)), max(ER(:,gap)) );

    %%% Calculate the bulk bands 
    bulk1(iq) = min( min(EL(:,gap+1)), min(ER(:,gap+1)) );
    bulk2(iq) = max( max(EL(:,gap)), max(ER(:,gap)) );

    %%%%% ==============================================================================
    %%%%% Calculate the edge states at the synthetic momentum q 
    %%% Array of energy 
    E_array = linspace(bulk2(iq)+epsilonE, bulk1(iq)-epsilonE,NE);
    %bulk1(iq)
    %bulk2(iq)

    %%% Scan the energy array E_array 
    for iE = 1:NE
        % The energy 
        E = E_array(iE);

        % Left-hand Hamiltonian 
        printf("Left-hand Hamiltonian \n");
        [WL,kL] = kPolyEig(E,q,omega,domega,eta,v1,v2,U,dU,W,dW,alpha,V,beta,dist,d0);

        % Right-hand Hamiltonian 
        printf("Right-hand Hamiltonian \n");
        [WR,kR] = kPolyEig(E,q,omega,-domega,eta,v1,v2,U,-dU,W,-dW,alpha,V,beta,dist,d0);

        %%% Assembly to the coefficient matrix 
        WW = CoefficientMatrix(WL,kL,WR,kR);

        %%% Calculate the absolute value of the determinant of WW 
        S = abs(det(WW)) 

        %%% If S<epsilon then add (k,E) to the edge state 
        if (S < epsilon) 
            edge_state = [edge_state;[q,E]];
        end % IF 

    end % iE-loop 

end % iq-loop 

%edge_state

%%%%% ===================================================================================
%%%%% Plot the figure 
figure(1)
plot(q_array,allEmin,'k'); hold on; 
plot(q_array,allEmax,'m');
plot(q_array,bulk1,color='red');
plot(q_array,bulk2,color='blue');
scatter(edge_state(:,1),edge_state(:,2));
hold off; 
xlabel('q');
ylabel('E');
saveas(1,'transmission2.png')