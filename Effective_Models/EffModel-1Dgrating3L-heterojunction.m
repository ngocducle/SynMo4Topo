% Parameters 
v = 0.3283 
U = 0.0207 
omega0 = 0.2413+0.0053 
domega = -omega0*0.55*(0.764-0.8)
dv = -v*0.37*(0.764-0.8)

Delta = 0.2*U
V1 = U 
V2 = U 

gap = 3 

Nk = 101
k_array = linspace(-0.25,0.25,Nk);
 
Nq = 101 
q_array = linspace(-0.03*2*pi,0.03*2*pi,Nq);

% Criterion for 0 
epsilon = 1e-3 

% Number of E values 
NE = 501 

% Small increment in band edge 
epsilonE = 1e-4 

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

%%%%% =================================================================================
%%%%% Scan over the q_array 
for iq = 1:Nq 
    %%% The synthetic momenta 
    q1 = q_array(iq)
    q2 = -q1

    %%% ===============================================================================
    %%% The left-hand side 
    U1 = U; 
    U2 = U+Delta;
    U3 = U;  

    % The matrix Ha (no velocity)
    Ha = [0,U1*exp(i*q1),V1,0,0,0;
          U1*exp(-i*q1),0,0,V1,0,0;
          V1,0,domega,U2,V2,0;
          0,V1,U2,domega,0,V2;
          0,0,V2,0,0,U3*exp(-i*q2);
          0,0,0,V2,U3*exp(i*q2),0];

    % The arrays of left bands 
    EL = zeros(Nk,6); 

    %%% Scan over k_array and diagonalize the Hamiltonian 
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

	EL(ik,:) = eigval;
    end 

    %%% ================================================================================
    %%% The right-hand side 
    U1 = U; 
    U2 = U-Delta;
    U3 = U;
    
    % The matrix Ha (no velocity)
    Ha = [0,U1*exp(i*q1),V1,0,0,0;
          U1*exp(-i*q1),0,0,V1,0,0;
          V1,0,domega,U2,V2,0;
          0,V1,U2,domega,0,V2;
          0,0,V2,0,0,U3*exp(-i*q2);
          0,0,0,V2,U3*exp(i*q2),0];

    % The arrays of right bands 
    ER = zeros(Nk,6); 

    %%% Scan the k_array and diagonalize the Hamiltonian 
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
	
	ER(ik,:) = eigval;
    end 
    
    %%% ==============================================================================
    % Calculate the obstructed bands 
    allEmax(iq) = max(min(EL(:,gap+1),min(ER(:,gap+1))));
    allEmin(iq) = min(max(EL(:,gap),max(ER(:,gap))));

    % Calculate the bulk bands 
    bulk1(iq) = min(min(EL(:,gap+1),min(ER(:,gap+1))));
    bulk2(iq) = max(max(EL(:,gap),max(ER(:,gap)))); 
 
    % Array of energy 
    E_array = linspace(bulk2(iq)+epsilonE,bulk1(iq)-epsilonE,NE);

    %%% Scan the energy array E_array 
    for iE = 1:NE
        % The energy 
        E = E_array(iE)  

        % The matrix of eigenstate 
        W = zeros(6,6);

        %%% Left-hand Hamiltonian 
        Ha = [0,U*exp(i*q1), V1,0,0,0;
              U*exp(-i*q1),0,0,V1,0,0;
              V1,0,domega,U+Delta,V2,0;
              0,V1,U+Delta,domega,0,V2;
              0,0,V2,0,0,U*exp(-i*q2);
              0,0,0,V2,U*exp(i*q2),0];

        % Solve the generalized eigenvalue problem A*x = lambda*B*x 
        A = E*eye-Ha; 
        B = -Hv;
        [kvecs,kvals] = eig(A,B);

	% Calculate the argument of the eigenvalues
	% If a k-eigenvalue has negative argument, +2*pi to convert to 
	% the range [0,2*pi[
	arg_array=zeros(6,1);
	for jE = 1:6
		kvals(jE,jE) = kvals(jE,jE)*exp(i*10e-9);
		arg_array(jE) = arg(kvals(jE,jE));

		if (arg_array(jE)<0) 
			arg_array(jE) = arg_array(jE)+2*pi;
		end 
	end 

	% Rearrange the k-eigenvalues in increasing order of argument
	[arg_list,ind] = sort(arg_array)
	kvecs = kvecs(:,ind);
	W(:,1:3) = kvecs(:,1:3)./norm(kvecs(:,1:3),'Fro','cols'); % Frobenius norm summed over columns

        %%% Right-hand Hamiltonian 
        Ha = [0,U*exp(i*q1),V1,0,0,0;
              U*exp(-i*q1),0,0,V1,0,0;
              V1,0,domega,U-Delta,V2,0;
              0,V1,U-Delta,domega,0,V2;
              0,0,V2,0,0,U*exp(-i*q2);
              0,0,0,V2,U*exp(i*q2),0];

        
        % Solve the generalized eigenvalue problem A*x = lambda*B*x 
        A = E*eye-Ha; 
        B = -Hv;
        [kvecs,kvals] = eig(A,B);

	% Calculate the argument of the eigenvalues and gives value in the interval [0,2*pi[
	arg_array = zeros(6,1);
	for jE = 1:6
		kvals(jE,jE) = kvals(jE,jE)*exp(i*10e-9);
		arg_array(jE) = arg(kvals(jE,jE));

		if (arg_array(jE)<0)
			arg_array(jE) = arg_array(jE)+2*pi;
		end 
	end 

	% Rearrange the k-eigenvalues in increasing order of argument
	[arg_list,ind] = sort(arg_array)
	kvecs = kvecs(:,ind);
		W(:,4:6) = kvecs(:,4:6)./norm(kvecs(:,4:6),'Fro','cols'); % Frobenius norm summed over columns

	% Caculate the determinant of W 
	S = abs(det(W))

	if (S < epsilon) 
		edge_state = [edge_state;[q1,E]];
	end 

    end 

end

edge_state

figure(1)
plot(q_array,allEmin+omega0,color='red');hold on; 
plot(q_array,allEmax+omega0,color='red');
plot(q_array,bulk1+omega0,color='red'); 
plot(q_array,bulk2+omega0,color='red');
scatter(edge_state(:,1),edge_state(:,2)+omega0);
hold off; 
%pause;
saveas(1,'transmission.png')
