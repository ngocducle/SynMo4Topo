%%% Solve for the edge states of the effective model 
%%% of 1D grating bilayer 
%%% OCTAVE version 

%%%%% ================================================================================
%%%%% FUNCTION: Ha 
%%%%% The zeroth degree of the Hamiltonian 
function Ha = H_a(q,omega1,U1,omega2,U2,V)
    Ha = [omega1,U1,V*exp(-i*pi*q),0;
          U1,omega1,0,V*exp(i*pi*q);
          V*exp(i*pi*q),0,omega2,U2;
          0,V*exp(-i*pi*q),U2,omega2];
end % function H_a 

%%%%% =================================================================================
%%%%% FUNCTION: H1 
%%%%% The first degree of the Hamiltoian 
function H1 = H_1(v1,v2)
    H1 = diag([v1,-v1,v2,-v2]);
end % function H_1 

%%%%% ==================================================================================
%%%%% FUNCTION: Hamiltonian 
%%%%% The Hamiltonian 
function H = Hamiltonian(k,q,omega1,v1,U1,omega2,v2,U2,V)
    Ha = H_a(q,omega1,U1,omega2,U2,V);
    H1 = H_1(v1,v2);
    H = Ha + H1*k;
end % function Hamiltonian 

%%%%% ===================================================================================
%%%%% FUNCTION: kPolyEig 
function [kvecs,kvals] = kPolyEig(E,q,omega1,v1,U1,omega2,v2,U2,V)
    %%% The matrix H0 
    Ha = H_a(q,omega1,U1,omega2,U2,V);
    H0 = Ha - E*eye(4);
    %H0 = [omega1-E,U1,V*exp(-i*pi*q),0;
    %      U1,omega1-E,0,V*exp(i*pi*q);
    %      V*exp(i*pi*q),0,omega2-E,U2;
    %      0,V*exp(-i*pi*q),U2,omega2-E];

    %%% The matrix H1 
    H1 = H_1(v1,v2);

    %%% Solve the polynomial eigenvalue problem 
    [kvecs,kvals] = polyeig(H0,H1); 
    kvals = kvals*exp(1e-9*i);

    %%% Rearrange kvals with increasing order of argument 
    args = zeros(4,1);

    for j = 1:4 
        args(j) = arg(kvals(j));

        if (args(j)<0) 
            args(j) = args(j) + 2*pi;
        end % IF 
    end % j-loop 

    args; 

    [arg_list,ind] = sort(args);
    arg_list; 
    ind; 
    kvals = kvals(ind);
    kvecs = kvecs(:,ind);
    kvecs = kvecs./norm(kvecs,'Fro','cols'); % Frobenius norm summed over columns 
end % function kPolyEig 

%%%%% ===================================================================================
%%%%% The parameters 
omega0 = 0.27815
v0 = 0.37602 
U0 = 0.02232 

omega1 = 0.27999
v1 = 0.40590
U1 = 0.02621 
V1 = 0.04511 

omega2 = 0.28010
v2 = 0.40454
U2 = 0.02142 
V2 = 0.04728 

%%% Gap where we calculate the edge states 
gap = 1 

%%% Criterion for 0 
epsilon = 1e-3

%%% Increment of band edge 
epsilonE = 1e-4

%%% Array of synthetic momenta 
Nq = 1001
q_array = linspace(0.0,1.0,Nq);

%%% Array of genuine momenta 
Nk = 101
k_array = linspace(-0.12,0.12,Nk);

%%% Number of energy values 
NE = 251

%%% Arrays of left and right bands 
EL = zeros(Nk,4);
ER = zeros(Nk,4);

%%% Arrays of bulk bands 
bulk1 = zeros(Nq,1);
bulk2 = zeros(Nq,1);

%%% Edge state 
Edge_state = []

count = 1; 

%%%%% =================================================================================
%%%%% Scan the array of synthetic momenta 
for iq = 1:Nq 
    %%% The synthetic momentum 
    q = q_array(iq);

    %%% Scan the k_array 
    for ik = 1:Nk 
        %%% The genuine momentum 
        k = k_array(ik);

        %%% ============================================================================
        %%% The left-hand side 
        %%% The Hamiltonian 
        HL = Hamiltonian(k,q,omega0,v0,U0,omega1,v1,U1,V1);

        %%% Diagonalize the Hamiltonian 
        [Vec,D] = eig(HL);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        %%% Save the energy eigenvalues to the array EL 
        EL(ik,:) = eigval; 

        %%% ============================================================================
        %%% The left-hand side 
        %%% The Hamiltonian 
        HR = Hamiltonian(k,q,omega0,v0,U0,omega2,v2,U2,V2);

        %%% Diagonalize the Hamiltonian 
        [Vec,D] = eig(HR);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        %%% Save the energy eigenvalues to the array EL 
        ER(ik,:) = eigval; 

    end % ik-loop 

    %%%%% ==============================================================================
    %%% Calculate the bulk bands 
    bulk1(iq) = min(min(EL(:,gap+1)),min(ER(:,gap+1)));
    bulk2(iq) = max(max(EL(:,gap)),max(ER(:,gap)));

    %%%%% ===============================================================================
    %%% Array of energies 
    E_array = linspace(bulk2(iq)-epsilonE,bulk1(iq)+epsilonE,NE);

    %%% Array of determinants 
    %S = zeros(NE,1);

    for iE = 1:NE 
        %count 

        % The value of the energy 
        E = E_array(iE);

        % The matrix of eigenstates 
        %WW = zeros(4);

        % The left-hand side 
        [WL,kL] = kPolyEig(E,q,omega0,v0,U0,omega1,v1,U1,V1);

        % The right-hand side 
        [WR,kR] = kPolyEig(E,q,omega0,v0,U0,omega2,v2,U2,V2);

        % Combine WL and WR to the matrix of eigenstates 
        WW(:,1:2) = WR(:,1:2);
        WW(:,3:4) = -WL(:,3:4);
        WW; 

        %WW = cat(2,WR(:,1:2),-WL(:,3:4))

        % The determinant 
        %S(iE) = abs(det(WW));

        S = abs(det(WW));

        if (S < epsilon)
            Edge_state = [Edge_state;[q,E]];
        end % IF 

        count = count + 1;

    end % iE-loop 

    %%% Scan E_array again 
    %for iE = 1:NE 
    %    if (S(iE)<epsilon)
    %        Edge_state = [Edge_state;[q,E_array(iE)]]; 
    %    end % IF 
    %end % iE-loop 

end % iq-loop 

%Edge_state
%S 

%%%%% ================================================================================
%%%%% Plot the figure 
figure(1)
plot(q_array,bulk1,color='red'); hold on; 
plot(q_array,bulk2,color='blue');
plot(Edge_state(:,1),Edge_state(:,2),'o');
hold off;
xlabel('q');
ylabel('E');
saveas(1,'transmission.png')