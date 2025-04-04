% The Hamiltonians 
function H = Hamiltonian(k,q,m,omega,v,U,W,ro,rv,rW,eta,alpha,V,beta,dist,d0)
    % Initialize the Hamiltonian 
    H = zeros(8);

    % The hat parameters 
    Vh = V*exp(-dist/d0);
    betah = beta*exp(-dist/d0);

    % The parameters 1 and 2 
    omega1 = omega*(1+ro*m);
    omega2 = omega*(1-ro*m);
    v1 = v*(1+rv*m);
    v2 = v*(1-rv*m);
    U1 = U*(1+m);
    U2 = U*(1-m);
    W1 = W*(1+rW*m);
    W2 = W*(1-rW*m);

    % The interlayer coupling functions
    V00 = -(Vh + betah*k + betah*k*k/sqrt(2));
    V01 = Vh + betah*k*k/sqrt(2);
    V10 = Vh + betah*k*k/sqrt(2);
    V11 = -(Vh - betah*k + betah*k*k/sqrt(2));

    % The block (1,1)
    H(1,1) = omega1 + eta + v1*k + v1*k*k/sqrt(2);
    H(1,2) = W1;
    H(1,3) = W1;
    H(1,4) = U1*(1+alpha);

    H(2,1) = W1; 
    H(2,2) = omega1 - eta + v1*k*k/sqrt(2);
    H(2,3) = U1*(1-alpha);
    H(2,4) = W1;

    H(3,1) = W1;
    H(3,2) = U1*(1-alpha);
    H(3,3) = omega1 - eta + v1*k*k/sqrt(2);
    H(3,4) = W1;

    H(4,1) = U1*(1+alpha);
    H(4,2) = W1;
    H(4,3) = W1;
    H(4,4) = omega1 + eta - v1*k + v1*k*k/sqrt(2);

    % The block (1,2)
    H(1,5) = V00*exp(-i*2*pi*q);
    H(2,6) = V01;
    H(3,7) = V10;
    H(4,8) = V11*exp(i*2*pi*q);

    % The block (2,1)
    H(5,1) = V00*exp(i*2*pi*q);
    H(6,2) = V01;
    H(7,3) = V10;
    H(8,4) = V11*exp(-i*2*pi*q);

    % The block (2,2)
    H(5,5) = omega2 + eta + v2*k + v2*k*k/sqrt(2);
    H(5,6) = W2;
    H(5,7) = W2;
    H(5,8) = U2*(1+alpha);

    H(6,5) = W2; 
    H(6,6) = omega2 - eta + v2*k*k/sqrt(2);
    H(6,7) = U2*(1-alpha);
    H(6,8) = W2;

    H(7,5) = W2;
    H(7,6) = U2*(1-alpha);
    H(7,7) = omega2 - eta + v2*k*k/sqrt(2);
    H(7,8) = W2;

    H(8,5) = U2*(1+alpha);
    H(8,6) = W2;
    H(8,7) = W2;
    H(8,8) = omega2 + eta - v2*k + v2*k*k/sqrt(2);
end % function Hamiltonian 