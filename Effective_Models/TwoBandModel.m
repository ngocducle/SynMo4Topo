% Reduce to 2 band model 
function H = Ham(k,q,omega,eta,v,U,W,alpha,V,beta0,dist,d0)
    Vh = V*exp(-dist/d0);
    betah = beta0*exp(-dist/d0);
    alphah = U*alpha; 
    qh = 2*pi*q; 

    dd = sqrt(Vh^2+4*W^2);
    cc = Vh/dd;
    ss = 2*W/dd; 

    d0 = omega + U - dd + v*k*k/sqrt(2) + Vh*sin(0.5*qh)^2 - betah*k*k*cos(0.5*qh)/sqrt(2);
    d1 = ss*v*k*k/sqrt(2) + ss*cc*Vh*sin(0.5*qh)^2 - betah*ss*cc*k*k*cos(0.5*qh)^2;
    d2 = 0.5*betah*ss*cc*k*sin(qh);
    d3 = (eta+alphah)*cc*cc + Vh*sin(0.5*qh)^2 + betah*k*k*sin(0.5*qh)^2/sqrt(2); 

    sigma0 = [1,0;0,1];
    sigmax = [0,1;1,0];
    sigmay = [0,-i;i,0];
    sigmaz = [1,0;0,-1];

    H = d0*sigma0 + d1*sigmax + d2*sigmay + d3*sigmaz;

end 

% Array of genuine momentum 
Nk = 101; 
k_array = linspace(-0.2,0.2,Nk);

% Arrya of synthetic momentum 
Nq = 101;
q_array = linspace(0,0.5,Nq);

%%% Parameters 
omega = 0.29780940 
eta = -0.0032 
v = 0.317 
U = -0.01503117
W = 0.001466 
alpha = 0.062 

d0 = 0.35 
dist = 0.1 
V = 0.038
beta = -0.3

for iq = 1:Nq
    for ik = 1:Nk  
        q = q_array(iq)
        k = k_array(ik)

        H = Ham(k,q,omega,eta,v,U,W,alpha,V,beta,dist,d0)

    end % ik-loop 
end % iq-loop 