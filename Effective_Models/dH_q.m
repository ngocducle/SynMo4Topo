% The derivative of the Hamiltonian with respect to q
function dHq = dH_q(k,q,V,beta,dist,d0)
    dHq = zeros(8);

    % The hat parameters 
    Vh = V*exp(-dist/d0);
    betah = beta*exp(-dist/d0);

    % The interlayer coupling 
    V00 = -(Vh + betah*k + betah*k*k/sqrt(2));
    V01 = Vh + betah*k*k/sqrt(2);
    V10 = Vh + betah*k*k/sqrt(2);
    V11 = -(Vh - betah*k + betah*k*k/sqrt(2));

    % The matrix elements
    dHq(1,5) = -i*2*pi*V00*exp(-i*2*pi*q);
    dHq(2,6) = 0;
    dHq(3,7) = 0;
    dHq(4,8) = i*2*pi*V11*exp(i*2*pi*q);

    dHq(5,1) = i*2*pi*V00*exp(i*2*pi*q);
    dHq(6,2) = 0;
    dHq(7,3) = 0;
    dHq(8,4) = -i*2*pi*V11*exp(-i*2*pi*q);

end % function dHq 