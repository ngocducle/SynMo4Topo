%%%%% Solve for the edge states of the effective model
%%%%% of 2D slab 2L
%%%%% MATLAB version
%%%%%
%%%%% The inputs of the code are:
%%%%% The values of omega1,v1,U1,W1,omega2,v2,U2,W2,V,dist,d0
%%%%% alpha,eta
%%%%% 
%%%%% The structure is the swapped one: the RHS is swapped from the LHS

clc
clear all
close all

%% The functions
%%%%% ==========================================================
%%%%% FUNCTION: Ha
%%%%% The zeroth-order of the Hamiltonian
function Ha = H_a(q,omega1,v1,U1,W1,omega2,v2,U2,W2,alpha,eta,V,dist,d0)
    Ha = zeros(8);

    Vt = V*exp(-dist/d0);

    alpha1 = alpha; 
    alpha2 = alpha;
    eta1 = eta;
    eta2 = eta;

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

%%%%% ==================================================================
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
end % function H_1

%%%%% ==================================================================
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
end % function H_2

%%%%% =====================================================================
%%%%% FUNCTION: Hamiltonian 
function H = Hamiltonian(k,q,omega1,v1,U1,W1,omega2,v2,U2,W2,eta,alpha,V,beta,dist,d0)

    Ha = H_a(q,omega1,v1,U1,W1,omega2,v2,U2,W2,alpha,eta,V,dist,d0);
    H1 = H_1(q,v1,v2,beta,dist,d0);
    H2 = H_2(q,v1,v2,beta,dist,d0);

    H = Ha + H1*k + H2*k*k;

end % function Hamiltonian

%%%%% =====================================================================
%%%%% FUNCTION: Derivative of the Hamiltonian
function dHk = dH_k(k,q,v1,v2,beta,dist,d0)
    
    H1 = H_1(q,v1,v2,beta,dist,d0);
    H2 = H_2(q,v1,v2,beta,dist,d0);

    dHk = H1 + 2*H2*k;

end % function derivative of the Hamiltonian

%%%%% =====================================================================
%%%%% FUNCTION: 16x16 Coefficient matrix
function [WW,Nbvariables] = CoefficientMatrix(q,v1,v2,beta,dist,d0,WL,kL,WR,kR,Kmax)
    % WW = zeros(16);

    %WW(1:8,1:8) = WL(:,9:16);
    %WW(1:8,9:16) = -WR(:,1:8);

    %for j = 1:8
    %    WW(9:16,j) = dH_k(kL(j+8),q,v1,v2,beta,dist,d0)*WL(:,j+8);
    %end % j-loop

    %for j = 1:8
    %    WW(9:16,j+8) = -dH_k(kR(j),q,v1,v2,beta,dist,d0)*WR(:,j);
    %end % j-loop

    %%% Tolerance 
    tol = 1e-7; % to eliminate pure real wave vectors

    %%% Count the number of wavenumbers k such that |Re(k)| <= Kmax 
    count = 0;

    for j = 1:8 
        if ( (abs(imag(kL(j+8))) > tol) & (abs(real(kL(j+8))) < Kmax) ) 
            count = count + 1;
        end 
    end % j-loop 

    for j = 1:8 
        if ( (abs(imag(kR(j))) > tol) & (abs(real(kR(j))) < Kmax) )
            count = count + 1;
        end 
    end % j-loop 

    Nbvariables = count; 

    %%% Initialize the matrix WW
    WW = zeros(16,Nbvariables);

    % Count the columns from beginning
    count = 0;

    for j = 1:8
        if ( (abs(imag(kL(j+8))) > tol) & (abs(real(kL(j+8))) < Kmax) )
            % Increase count by 1 
            count = count + 1;

            % The continuity condition for wavefunction 
            WW(1:8,count) = WL(:,j+8);

            % The continuity condition for current density 
            WW(9:16,count) = dH_k(kL(j+8),q,v1,v2,beta,dist,d0)*WL(:,j+8);
        end 
    end 

    for j = 1:8
        if ( (abs(imag(kR(j))) > tol) & (abs(real(kR(j))) < Kmax) )
            % Increase count by 1 
            count = count + 1; 

            % The continuity condition for wavefunction 
            WW(1:8,count) = -WR(:,j);

            % The continuity condition for current density 
            WW(9:16,count) = -dH_k(kR(j),q,v1,v2,beta,dist,d0)*WR(:,j);
        end 
    end 

end % function CoefficientMatrix

%%%%% =====================================================================
%%%%% FUCNTION: Solve the polynomial eigenvalue problem 
function [kvecs,kvals] = kPolyEig(E,q,omega1,v1,U1,W1,omega2,v2,U2,W2,eta,alpha,V,beta,dist,d0)
    %%% Establish the matrices H0, H1, H2 
    Ha = H_a(q,omega1,v1,U1,W1,omega2,v2,U2,W2,alpha,eta,V,dist,d0);
    H0 = Ha - E*eye(8);
    H1 = H_1(q,v1,v2,beta,dist,d0);
    H2 = H_2(q,v1,v2,beta,dist,d0);

    %%% Polynomial diagonalization 
    [kvecs,kvals] = polyeig(H0,H1,H2);

    %%% We move the real eigenvalues to have infinitely small positive imaginary parts,
    %%% phase and argument, so we multiply kvals by exp(i*1e-9)
    %%% In MATLAB/OCTAVE, the POSITIVE REAL eigenvalues are given by 
    %%% a - {small value}*i and have phase 6.28... 
    %%% We need to convert their argument/phase to infinitisemal positive value 
    %%% Otherwise, the program will be wrong 
    kvals = kvals*exp(i*1e-9);

    %%% Move the arguments to the range 0 <= argument <= 2*pi
    arg_array = zeros(16,1);

    for j = 1:16 
        arg_array(j) = angle(kvals(j)); % The function angle gives values in the range [-pi,pi]
        
        if (arg_array(j)<0)
            arg_array(j) = arg_array(j) + 2*pi;
        end % IF-arg_array 
   
    end % j-loop 

    %%% Rearrange kvecs and kvals in increasing order of arguments
    [arg_list,ind] = sort(arg_array);
    kvals = kvals(ind);
    kvecs = kvecs(:,ind);
    kvecs = kvecs./vecnorm(kvecs); % normalize the columns

end % function kPolyEig 

%%%%% =====================================================================
%% The MAIN program goes here 
%%%%% =====================================================================
%%%%% Parameters 
% h = 0.35, b1 = 0.30 
omega1 = 0.28832692   
v1 = 0.30552363   
U1 = -0.01112476   
W1 = 0.00062376

% h = 0.35, b2 = 0.46 
omega2 = 0.31005337   
v2 = 0.32069812   
U2 = -0.01938797   
W2 = 0.00289523

% The hole parameter 
eta = 0.0032 
alpha = -0.062 

% The interlayer coupling 
d0 = 0.35
dist = 0.1 
V = 0.039 
beta = -0.3 

%%% The gap we calculate the edge state 
gap = 1 

%%% The array of genuine momenta 
Nk = 101 
Kmax = 0.1
k_array = linspace(-Kmax,Kmax,Nk);

%%% The array of synthetic momenta 
Nq = 101 
Qmax = 0.5
q_array = linspace(-Qmax,Qmax,Nq);

%%% Criterion for 0 
epsilon = 3e-3; % sqrt(sqrt(sqrt(sqrt(1e-55))));

%%% Number of E values to scan 
NE = 2001 

%%% Small increment in band edge 
epsilonE = 1e-5 

%%% Initialize the edge states to be empty 
edge_state = []; 

%%%%% =====================================================================
%%%%% Scan the q_array 
for iq = 1:Nq
    %%% The synthetic momentum 
    q = q_array(iq);

    %%% The arrays of left and right bands 
    EL = zeros(Nk,8);
    ER = zeros(Nk,8);

    %%%%% =================================================================
    for ik = 1:Nk 
        %%%%% The genuine momentum 
        k = k_array(ik);

        %%%%% =================================================================
        %%%%% The left-hand side 
        %%% The Hamiltonian
        H = Hamiltonian(k,q,omega1,v1,U1,W1,omega2,v2,U2,W2,eta,alpha,V,beta,dist,d0);

        %%% Diagonalize the Hamiltonian 
        [Vec,D] = eig(H);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        %%% Save the energy eigenvalues to the array EL
        EL(ik,:) = eigval;

        %%%%% =============================================================
        %%%%% The right-hand side 
        %%% The Hamiltonian 
        H = Hamiltonian(k,q,omega2,v2,U2,W2,omega1,v1,U1,W1,eta,alpha,V,beta,dist,d0);

        %%% Diagonalize the Hamiltonian 
        [Vec,D] = eig(H);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        %%% Save the energy eigenvalues to the array EL 
        ER(ik,:) = eigval;

    end % ik-loop

    %%%%% =================================================================
    %%% Calculate the obstructed bands 
    allEmax(iq) = max( min(EL(:,gap+1)), min(ER(:,gap+1)) );
    allEmin(iq) = min( max(EL(:,gap)), max(ER(:,gap)) );

    %%% Calculate the bulk bands 
    bulk1(iq) = min( min(EL(:,gap+1)), min(ER(:,gap+1)) );
    bulk2(iq) = max( max(EL(:,gap)), max(ER(:,gap)));

    %%%%% =================================================================
    %%%%% Calculate the edge states at the synthetic momentum q 
    %%% Array of energy 
    E_array = linspace(bulk2(iq)+epsilonE,bulk1(iq)-epsilonE,NE);

    S = zeros(NE,1);

    %%% Scan the energy array E_array 
    for iE = 1:NE
        % The energy 
        E = E_array(iE);

        % Left-hand Hamiltonian 
        [WL,kL] = kPolyEig(E,q,omega1,v1,U1,W1,omega2,v2,U2,W2,eta,alpha,V,beta,dist,d0)

        % Right-hand Hamiltonian
        [WR,kR] = kPolyEig(E,q,omega2,v2,U2,W2,omega1,v1,U1,W1,eta,alpha,V,beta,dist,d0)

        % Assembly to the coefficient matrix
        [WW,Nbvariables] = CoefficientMatrix(q,v1,v2,beta,dist,d0,WL,kL,WR,kR,0.7);
        Nbvariables

        % Calculate the absolute value of the determinant of WW 
        % S(iE) = abs(det(WW));
        % det(WW)
        % svd(WW)

        % If S < epsilon then add (k,E) to the edge state 
        % if (S(iE) < epsilon) 
        %    edge_state = [edge_state;[q,E,S(iE)]];
        % end % IF

        % Calculate the rank of WW
        rankW = rank(WW,epsilon)

        % The kernel of WW 
        % kernel = null(WW,epsilon)

        % If rank(WW) < Nbvariables (the number of variables) then add (k,E) to
        % the edge state
        if (rankW < Nbvariables) 
            edge_state = [edge_state;[q,E,rankW]];
        end % IF

    end % iE-loop

end % iq-loop

%%%%% =====================================================================
%% Plot the figure
figure(1)
plot(q_array,allEmin,'k'); hold on;
plot(q_array,allEmax,'m');
plot(q_array,bulk1,'r');
plot(q_array,bulk2,'b');
scatter(edge_state(:,1),edge_state(:,2),'green','filled');
hold off;
xlabel('q');
ylabel('E');
title(['tol = ',num2str(epsilon),', Nq = ',num2str(Nq),', NE = ',num2str(NE)]);
saveas(1,'transmission8.png');