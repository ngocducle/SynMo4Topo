% The derivative of the Hamiltonian with respect m 
function dHm = dH_m(k,m,omega,v,U,W,ro,rv,rW,alpha,eta)
    dHm = zeros(8);

    % Block (1,1)
    dHm(1,1) = omega*ro + v*rv*k + v*rv*k*k/sqrt(2);
    dHm(1,2) = W*rW;
    dHm(1,3) = W*rW;
    dHm(1,4) = U*(1+alpha);

    dHm(2,1) = W*rW;
    dHm(2,2) = omega*ro + v*rv*k*k/sqrt(2);
    dHm(2,3) = U*(1-alpha);
    dHm(2,4) = W*rW;

    dHm(3,1) = W*rW;
    dHm(3,2) = U*(1-alpha);
    dHm(3,3) = omega*ro + v*rv*k*k/sqrt(2);
    dHm(3,4) = W*rW;

    dHm(4,1) = U*(1+alpha);
    dHm(4,2) = W*rW;
    dHm(4,3) = W*rW;
    dHm(4,4) = omega*ro - v*rv*k + v*rv*k*k/sqrt(2);

    % Block (2,2)
    for row = 5:8
        for column = 5:8 
            dHm(row,column) = -dHm(row-4,column-4);
        end % column-loop 
    end % row-loop 
    
end % function dHm 