%%% Solve for the edge states of the effective model
%%% of 1D grating bilayer 
%%% OCTAVE version 

%%%%% ================================================================================
%%%%% FUNCTION: Ha 
%%%%% The zero order of the Hamiltonian 
function Ha = H_a(q,omega0,U,Delta,V)

    Ha = [omega0,U+Delta,V*exp(-i*pi*q),0;
        U+Delta,omega0,0,V*exp(i*pi*q);
        V*exp(i*pi*q),0,omega0,U-Delta;
        0,V*exp(-i*pi*q),U-Delta,omega0 
        ]; 

end % FUNCTION: Ha 

%%%%% ==================================================================================
%%%%% FUNCTION: H1 
%%%%% The linear order of the Hamiltonian 
function H1 = H_1(v)
    H1 = diag([v,-v,v,-v]);
end % FUNCTION: H1 

%%%%% ===================================================================================
%%%%% FUNCTION: Hamiltonian 
%%%%% The Hamiltonian 
function H = Hamiltonian(k,q,omega0,v,U,Delta,V)
    Ha = H_a(q,omega0,U,Delta,V);
    H1 = H_1(v);
    H = Ha + H1*k;
end % FUNCTION: Hamiltonian

%%%%% ====================================================================================
%%%%% FUNCTION: kPolyEig 
function [kvecs,kvals] = kPolyEig(E,q,omega0,v,U,Delta,V)
    %%% The matrix alpha 
    Ha = H_a(q,omega0,U,Delta,V); 
    H0 = Ha - E*eye(4);

    %%% The matrix H1 
    H1 = H_1(v); 

    %%% Solve the polynomial eigenvalue problem 
    [kvecs,kvals] = polyeig(H0,H1); 

    %%% Move the arguments to the range 0 <= argument <= 2*pi 
    arg_array = zeros(4,1);

    for j = 1:4 
        arg_array(j) = arg(kvals(j));

        if (arg_array(j) < 0)
            arg_array(j) = arg_array(j) + 2*pi;  
        end % IF 
    end % j-loop

    %%% Rearrange kvecs and kvals in increasing ordre of arguments 
    [arg_list,ind] = sort(arg_array);
    kvals = kvals(ind);
    kvecs = kvecs(:,ind);
    kvecs = kvecs./norm(kvecs,'Fro','cols'); % Frobenius norm summed over columns  
end % FUNCTION: kPolyEig 

%%%%% ====================================================================================
%%%%% Parameters 
omega0 = 0.25 
v = 0.30
U = 0.02
V = 0.05 
Delta = 0.2*U 

%%%%% Array of k 
Kmax = 0.05
Nk = 101 
k_array = linspace(-Kmax,Kmax,Nk);

%%%%% Array of q 
Qmax = 0.3
Nq = 101 
q_array = linspace(-Qmax,Qmax,Nq);

%%% The gap we calculate the edge states 
gap = 1

%%% Criterion for 0 
epsilon = 1e-3

%%% Number of E values to scan 
NE = 501 

%%% Small increment in band edge 
epsilonE = 1e-4 

%%% Intialize the edge states to be empty 
edge_state = []

%%% The arrays of obstructed and bulk bands 
allEmax = zeros(Nq,1);
allEmin = zeros(Nq,1);
bulk1 = zeros(Nq,1);
bulk2 = zeros(Nq,1);

%%%%% ===================================================================================
%%%%% Scan over the q_array 
for iq = 1:Nq 
    %%%%% The synthetic momentum 
    q = q_array(iq);

    %%%%% The arrays of left and right bands
    EL = zeros(Nk,4);
    ER = zeros(Nk,4);

    %%%%% ===============================================================================
    %%%%% Scan the k_array 
    for ik = 1:Nk 
        %%%%% The genuine momenta 
        k = k_array(ik);

        %%%%% ============================================================================
        %%%%% The left-hand side 
        %%% The Hamiltonian 
        H = Hamiltonian(k,q,omega0,v,U,Delta,V);

        %%% Diagonalize the Hamiltonian
        [Vec,D] = eig(H);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        %%% Save the energy eigenvalues to the array EL 
        EL(ik,:) = eigval;

        %%%%% ===========================================================================
        %%%%% The right-hand side 
        %%% The Hamiltonian 
        H = Hamiltonian(k,q,omega0,v,U,-Delta,V);

        %%% Diagonalize the Hamiltonian 
        [Vec,D] = eig(H);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        %%% Save the energy eigenvalues to the array ER
        ER(ik,:) = eigval;

    end % ik-loop

    %%%%% ===============================================================================
    %%% Calculate the obstructed bands 
    allEmax(iq) = max( min(EL(:,gap+1)), min(ER(:,gap+1)) );
    allEmin(iq) = min( max(EL(:,gap)), max(ER(:,gap)) );

    %%% Calculate the bulk bands 
    bulk1(iq) = min( min(EL(:,gap+1)), min(ER(:,gap+1)) );
    bulk2(iq) = max( max(EL(:,gap)), max(ER(:,gap)) );

    %%%%% ===============================================================================
    %%%%% Calculate the edge states at the synthetic momentum q 
    %%% Array of energies 
    E_array = linspace(bulk2(iq)+epsilonE,bulk1(iq)-epsilon,NE);
    S = zeros(NE);

    %%% Scan the energy array E_array 
    for iE = 1:NE 
        % The energy 
        E = E_array(iE);

        % The matrix of eigenstates 
        W = zeros(4);

        % Left-hand Hamiltonian 
        printf("Left-hand Hamiltonian \n");
        [WL,kL] = kPolyEig(E,q,omega0,v,U,Delta,V); 

        % Right-hand Hamiltonian 
        printf("Right-hand Hamiltonian \n");
        [WR,kR] = kPolyEig(E,q,omega0,v,U,-Delta,V); 

        %%% Combine WL and WR to the matrix of eigenstates 
        W(:,1:2) = WL(:,3:4);
        W(:,3:4) = -WR(:,1:2);

        %%% The determinant of W 
        S(iE) = abs(det(W));

        %if (S < epsilon)
        %    edge_state = [edge_state;[q,E]];
        %end % IF 

    end % iE-loop 

    %%% Minimum of the determinant 
    [minS,imin]=min(S);
    if (minS<epsilon)
        edge_state = [edge_state;[q,E_array(imin)]];
    end % IF 

end % iq-loop 

%%%%% ===================================================================================
%%%%% Plot the figure 
figure(1)
plot(q_array,allEmin,color='red'); hold on;
plot(q_array,allEmax,color='blue');
plot(q_array,bulk1,color='red');
plot(q_array,bulk2,color='blue');
scatter(edge_state(:,1),edge_state(:,2));
hold off; 
xlim([-Qmax,Qmax]);
xlabel('q');
ylabel('E');
saveas(1,'transmission.png')
