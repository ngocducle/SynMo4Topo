% The derivative of the Hamiltonian with respect to k 
function dHk = dH_k(k,q,m,v,rv,V,beta,dist,d0)
    dHk = zeros(8);

    % The hat parameters 
    Vh = V*exp(-dist/d0);
    betah = beta*exp(-dist/d0);

    % Calculate v1 and v2 
    v1 = v*(1+rv*m);
    v2 = v*(1-rv*m);

    % The derivatives of the interlayer coupling 
    dV00 = -betah*(1+sqrt(2)*k);
    dV01 = sqrt(2)*betah*k;
    dV10 = sqrt(2)*betah*k;
    dV11 = betah*(1-sqrt(2)*k);

    % Block (1,1)
    dHk(1,1) = v1*(1+sqrt(2)*k);
    dHk(2,2) = sqrt(2)*v1*k;
    dHk(3,3) = sqrt(2)*v1*k;
    dHk(4,4) = v1*(-1+sqrt(2)*k);

    % Block (1,2)
    dHk(1,5) = dV00*exp(-i*2*pi*q);
    dHk(2,6) = dV01;
    dHk(3,7) = dV10;
    dHk(4,8) = dV11*exp(i*2*pi*q);

    % Block (2,1)
    dHk(5,1) = dV00*exp(i*2*pi*q);
    dHk(6,2) = dV01;
    dHk(7,3) = dV10;
    dHk(8,4) = dV11*exp(-i*2*pi*q);

    % Block (2,2)
    dHk(5,5) = v2*(1+sqrt(2)*k);
    dHk(6,6) = sqrt(2)*v2*k;
    dHk(7,7) = sqrt(2)*v2*k;
    dHk(8,8) = v2*(-1+sqrt(2)*k);

end % function dHk  