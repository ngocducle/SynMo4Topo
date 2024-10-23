% Parameters 
v = 0.3283 
U = 0.0207 
omega0 = 0.2413+0.0053 
domega = -omega0*0.55*(0.764-0.8)
dv = -v*0.37*(0.764-0.8)

Delta = 0.2*U 
V1 = U 
V2 = U 

% The gap we calculate the edge state 
gap = 3 

% The array of genuine momemta 
Nk = 101 
k_array = linspace(-0.25,0.25,Nk); 

% The array of synthetic momenta 
Nq = 101
q_array = linspace(-0.03*2*pi,0.03*2*pi,Nq);

% Criterion for 0 
epsilon = 1e-5

% Number of E values 
NE = 501 

% Small incriment in band edge 
epsilonE = 1e-2

% Diagonal matrix of v 
Hv = diag([v,-v,v+dv,-v-dv,v,-v]);

% Array of obstructed bands (corresponding to the transmission spectrum)
allEmax = zeros(Nq,1);
allEmin = zeros(Nq,1);

% Array of bulk bands (closest bands)
bulk1 = zeros(Nq,1);
bulk2 = zeros(Nq,1);

% Initialize the array of edge state (k,E)
edge_state = [];

%%%%% ===================================================================================
%%%%% FUNCTION: Ha 
function Ha = H_a(q1,q2,domega,U1,U2,U3,V1,V2)
    Ha = [0,U1*exp(i*q1),V1,0,0,0;
          U1*exp(-i*q1),0,0,V1,0,0;
          V1,0,domega,U2,V2,0;
          0,V1,U2,domega,0,V2;
          0,0,V2,0,0,U3*exp(-i*q2);
          0,0,0,V2,U3*exp(i*q2),0];
end 

%%%%% ===================================================================================
%%%%% FUNCTION: Calculate the k-eigenvalues and k-eigenvectors for one side 
function [kvecs,kvals] = kPolyEig(E,q1,q2,domega,U,Delta,V1,V2,Hv)
    %%% The one-side bulk Hamiltonian 
    Ha = H_a(q1,q2,domega,U,U+Delta,U,V1,V2);

    %%% Solve the generalized eigenvalue problem A*x=k*B*x 
    % Method 1: Using the eig function
    %A = Ha-E*eye(6);
    %B = -Hv;  
    %[kvecs,kvals] = eig(A,B)

    % Method 2: Using the polyeig function 
    A0 = Ha-E*eye(6);
    A1 = Hv;
    [kvecs,kvals] = eig(A0,A1);

    %%% Move the arguments to the range 0 <= argument < 2*pi
    arg_array = zeros(6,1); 
    
    for j=1:6 
        arg_array(j) = arg(kvals(j,j));

        if (arg_array(j) < 0)
            arg_array(j) = arg_array(j)+2*pi;
        end % IF
    end % j-loop      

    %%% Rearrange kvecs and kvals in increasing order of argument
    [arg_list,ind] = sort(arg_array)
    kvals = kvals(ind,ind) % ATTENTION! kvals is a 6x6 matrix 
    kvecs = kvecs(:,ind);
    kvecs = kvecs./norm(kvecs,'Fro','cols'); % Frobenius norm summed over columns 
end 

%%%%% ===================================================================================
%%%%% Scan over the q_array 
for iq = 1:Nq 
    %%%%% The synthetic momneta 
    q1 = q_array(iq);
    q2 = -q1;

    %%% ==================================================================================
    %%% The left-hand side

    % The matrix Ha 
    Ha = H_a(q1,q2,domega,U,U+Delta,U,V1,V2); 

    % The arrays of left bands 
    EL = zeros(Nk,6);

    % Scan over k_array and diagonalize the Hamiltonian 
    for ik = 1:Nk 
        % The genuine momentum 
        k = k_array(ik);

        % The Hamiltonian 
        H = Hv*k + Ha;

        % Diagonalize the Hamiltonian 
        [Vec,D] = eig(H);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        % Save the eigenvalues to the array EL 
        EL(ik,:) = eigval;
    end  % ik-loop 

    %%% ==================================================================================
    %%% The right-hand side

    % The matrix Ha 
    Ha = H_a(q1,q2,domega,U,U-Delta,U,V1,V2); 

    % The arrays of right bands 
    ER = zeros(Nk,6);

    % Scan over k_array and diagonalize the Hamiltonian 
    for ik = 1:Nk 
        % The genuine momentum 
        k = k_array(ik);

        % The Hamiltonian 
        H = Hv*k + Ha;

        % Diagonalize the Hamiltonian 
        [Vec,D] = eig(H);
        [eigval,ind] = sort(diag(D));
        eigval = real(eigval);
        Ds = D(ind,ind);
        eigvec = Vec(:,ind);

        % Save the eigenvalues to the array ER 
        ER(ik,:) = eigval;
    end  % ik-loop

    %%%%% ==============================================================================
    %%% Calculate the obstructed bands 
    allEmax(iq) = max( min(EL(:,gap+1)), min(ER(:,gap+1)) );
    allEmin(iq) = min( max(EL(:,gap)), max(ER(:,gap))  );

    %%% Calculate the bulk bands 
    bulk1(iq) = min( min(EL(:,gap+1)), min(ER(:,gap+1)) ); 
    bulk2(iq) = max( max(EL(:,gap)), max(ER(:,gap)) );

    %%%%% ===============================================================================
    %%%%% Calculate the edge states at the synthetic momentum q 
    %%% Array of energy 
    E_array = linspace(bulk2(iq)+epsilonE,bulk1(iq)-epsilonE,NE);

    %%%%% Scan the energy array E_array 
    for iE = 1:NE 
        % The energy 
        E = E_array(iE); 

        % The matrix of eigenstates 
        W = zeros(6);

        %%% Left-hand Hamiltonian 
        %Ha = H_a(q1,q2,domega,U,U+Delta,U,V1,V2);

        % Solve the generalized eigenvalue problem A*x=k*B*x 
        %A = Ha-E*eye(6);
        %B = -Hv; 
        %[kvecs,kvals] = eig(A,B); 

        %%% Left-hand Hamiltonian 
        printf("Left-hand Hamiltonian\n");
        [WL,kL] = kPolyEig(E,q1,q2,domega,U,Delta,V1,V2,Hv); 
        %kL
        %WL 

        %%% Right-hand Hamiltonian
        printf("Right-hand Hamiltonian\n");
        [WR,kR] = kPolyEig(E,q1,q2,domega,U,-Delta,V1,V2,Hv); 
        %kR  
        %WR  

        %%% Combine WL and WR to the matrix of eigenstates 
        W(:,1:3) = WL(:,1:3);
        W(:,4:6) = -WR(:,4:6);

        %%% The determinant of W 
        S = abs(det(W))

        if (S < epsilon)
            edge_state = [edge_state;[q1,E]];
        end % IF 

    end 
end % iq-loop 

%%%%% ===================================================================================
%%%%% Plot the figure 
figure(1)
plot(q_array,allEmin+omega0,color='red'); hold on;
plot(q_array,allEmax+omega0,color='red');
plot(q_array,bulk1+omega0,color='red');
plot(q_array,bulk2+omega0,color='red');
scatter(edge_state(:,1),edge_state(:,2)+omega0);
hold off;
saveas(1,'transmission.png')