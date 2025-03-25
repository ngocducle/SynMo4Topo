clc
clear all
close all

%% The FUNCTIONS
function H = Hamiltonian(k,q,omega1,v1,U1,W1,omega2,v2,U2,W2,V,dist,d0,beta,alpha,eta)
    % The Hamiltonian 
    H = zeros(8);

    % The hat-parameters 
    Vh = V*exp(-dist/d0);
    betah = beta*exp(-dist/d0);

    % The block (1,1)
    H(1,1) = eta + v1*k + omega1 + v1*k*k/sqrt(2);
    H(1,2) = W1;
    H(1,3) = W1;
    H(1,4) = U1*(1+alpha);

    H(2,1) = W1;
    H(2,2) = -eta + omega1 + v1*k*k/sqrt(2);
    H(2,3) = U1*(1-alpha);
    H(2,4) = W1;

    H(3,1) = W1;
    H(3,2) = U1*(1-alpha);
    H(3,3) = -eta + omega1 + v1*k*k/sqrt(2);
    H(3,4) = W1;

    H(4,1) = U1*(1+alpha);
    H(4,2) = W1;
    H(4,3) = W1;
    H(4,4) = eta - v1*k + omega1 + v1*k*k/sqrt(2);

    % Block (1,2)
    H(1,5) = -(Vh + betah*k + betah*k*k/sqrt(2))*exp(-i*2.0*pi*q);
    H(2,6) = Vh + betah*k*k/sqrt(2);
    H(3,7) = Vh + betah*k*k/sqrt(2);
    H(4,8) = -(Vh - betah*k + betah*k*k/sqrt(2))*exp(i*2.0*pi*q);

    % Block (2,1)
    H(5,1) = -(Vh + betah*k + betah*k*k/sqrt(2))*exp(i*2.0*pi*q);
    H(6,2) = Vh + betah*k*k/sqrt(2);
    H(7,3) = Vh + betah*k*k/sqrt(2);
    H(8,4) = -(Vh - betah*k + betah*k*k/sqrt(2))*exp(-i*2.0*pi*q);

    % Block (2,2)
    H(5,5) = eta + v2*k + omega2 + v2*k*k/sqrt(2);
    H(5,6) = W2;
    H(5,7) = W2;
    H(5,8) = U2*(1+alpha);

    H(6,5) = W2;
    H(6,6) = -eta + omega2 + v2*k*k/sqrt(2);
    H(6,7) = U2*(1-alpha);
    H(6,8) = W2;

    H(7,5) = W2;
    H(7,6) = U2*(1-alpha);
    H(7,7) = -eta + omega2 + v2*k*k/sqrt(2);
    H(7,8) = W2;

    H(8,5) = U2*(1+alpha);
    H(8,6) = W2;
    H(8,7) = W2;
    H(8,8) = eta - v2*k + omega2 + v2*k*k/sqrt(2);
end % function Hamiltonian

%% The MAIN program goes here 
% The parameters 
omega1 = 0.31005337   
v1 = 0.32069812   
U1 = -0.01938797   
W1 = 0.00289523

omega2 = 0.28832692   
v2 = 0.30552363   
U2 = -0.01112476   
W2 = 0.00062376

V = 0.039 
beta = -0.3
dist = 0.1
d0 = 0.35 

alpha = -0.062
eta = 0.0032 

% The arrays of genuine and synthetic momenta 
Nk = 31
Kmax = 0.1
k_array = linspace(-Kmax,Kmax,Nk);

Nq = 21 
Qmax = 0.5
q_array = linspace(-Qmax,Qmax,Nq);

%% Calculate the dispersion 
E1_array = zeros(Nk,Nq);
E2_array = zeros(Nk,Nq);
E3_array = zeros(Nk,Nq);
E4_array = zeros(Nk,Nq);
E5_array = zeros(Nk,Nq);
E6_array = zeros(Nk,Nq);
E7_array = zeros(Nk,Nq);
E8_array = zeros(Nk,Nq);

% Scan the arrays of genuine and synthetic momenta 
for ik = 1:Nk
    k = k_array(ik);

    for iq = 1:Nq 
        q = q_array(iq);

        % The Hamiltonian 
        H = Hamiltonian(k,q,omega1,v1,U1,W1,omega2,v2,U2,W2,V,dist,d0,beta,alpha,eta);

        % Diagonalize the Hamiltonian
        [Vec,D] = eig(H);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        % Save the energy eigenvalues to array 
        E1_array(ik,iq) = eigval(1);
        E2_array(ik,iq) = eigval(2);
        E3_array(ik,iq) = eigval(3);
        E4_array(ik,iq) = eigval(4);
        E5_array(ik,iq) = eigval(5);
        E6_array(ik,iq) = eigval(6);
        E7_array(ik,iq) = eigval(7);
        E8_array(ik,iq) = eigval(8);
        
    end % iq-loop 
end % ik-loop 

%% Plot the dispersion surfaces 
figure(1);
surf(k_array,q_array,E1_array');hold on;
surf(k_array,q_array,E2_array');
hold off;
colormap summer;
xticks([-Kmax,-0.5*Kmax,0,0.5*Kmax,Kmax]);
yticks([-0.5,-0.25,0,0.25,0.5]);
xlabel('k','FontSize',18);
ylabel('q','Fontsize',18);
zlabel('E','FontSize',18);

