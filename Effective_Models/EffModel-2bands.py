import numpy as np
import scipy 
import scipy.linalg as sla
import cmath 
import matplotlib.pyplot as plt

##### ======================================================================
##### FUNCTION: 2-band model
def TwoBandModel(omega,eta,v,U,W,alpha,d0,dist,V,beta,k,q):
    ### The Hamiltonian
    H = np.zeros((2,2),dtype=complex)

    Vh = V*np.exp(-dist/d0)
    betah = beta*np.exp(-dist/d0)

    ### The energy eigenvalues of the unperturbed Hamiltonian
    E = np.zeros(8)
    E[0] = -eta - Vh - U*(1-alpha)
    E[1] = -eta + Vh - U*(1-alpha)
    E[2] = eta - Vh - U*(1+alpha)
    E[3] = eta + Vh - U*(1+alpha)
    E[4] = U - np.sqrt( (U*alpha-Vh+eta)*(U*alpha-Vh+eta) + 4*W*W )
    E[5] = U + np.sqrt( (U*alpha-Vh+eta)*(U*alpha-Vh+eta) + 4*W*W )
    E[6] = U - np.sqrt( (U*alpha+Vh+eta)*(U*alpha+Vh+eta) + 4*W*W )
    E[7] = U + np.sqrt( (U*alpha+Vh+eta)*(U*alpha+Vh+eta) + 4*W*W )

    ### The functions of k and q 
    hp = v*k + v*k*k/np.sqrt(2)
    h0 = v*k*k/np.sqrt(2)
    hm = -v*k + v*k*k/np.sqrt(2)

    omegap = -Vh*(cmath.exp(-1j*q)-1) - (betah*k+betah*k*k/np.sqrt(2))*cmath.exp(-1j*q)
    omega0 = betah*k*k/np.sqrt(2)
    omegam = -Vh*(cmath.exp(1j*q)-1) - (-betah*k+betah*k*k/np.sqrt(2))*cmath.exp(1j*q)

    ### The parameters A, B, C, D
    A = (alpha*U - Vh + eta + np.sqrt( (alpha*U - Vh + eta)**2 + 4*W**2 )) / (2*W)
    B = (alpha*U - Vh + eta - np.sqrt( (alpha*U - Vh + eta)**2 + 4*W**2 )) / (2*W)
    C = (alpha*U + Vh + eta + np.sqrt( (alpha*U + Vh + eta)**2 + 4*W**2 )) / (2*W)
    D = (alpha*U + Vh + eta - np.sqrt( (alpha*U + Vh + eta)**2 + 4*W**2 )) / (2*W)

    ### The cos and sin 
    ca = 1/np.sqrt(1+A*A)
    sa = A/np.sqrt(1+A*A)
    cb = 1/np.sqrt(1+B*B)
    sb = B/np.sqrt(1+B*B)
    cc = 1/np.sqrt(1+C*C)
    sc = C/np.sqrt(1+C*C)
    cd = 1/np.sqrt(1+D*D)
    sd = D/np.sqrt(1+D*D)

    ### The matrix elements
    Hp = np.zeros((8,8),dtype=complex)

    Hp[4,0] = 0
    Hp[4,1] = 0
    Hp[4,2] = 0.5*ca*(hm-hp+np.real(omegam)-np.real(omegap))
    Hp[4,3] = 0.5*ca*(omegam-omegap-np.real(omegam)+np.real(omegap))
    Hp[4,4] = 0.5*( ca*ca*(hp+hm) + 2*(h0+omega0)*sa*sa + ca*ca*np.real(omegap+omegam) )
    Hp[4,5] = 0.5*( ca*cb*(hp+hm) + 2*(h0+omega0)*sa*sb + ca*cb*np.real(omegap+omegam) )
    Hp[4,6] = 0.5*1j*ca*cc*np.imag(omegap+omegam)
    Hp[4,7] = 0.5*1j*ca*cd*np.imag(omegap+omegam)

    #Hp7 = np.zeros(8,dtype=complex)
    Hp[6,0] = 0
    Hp[6,1] = 0 
    Hp[6,2] = -0.5*1j*cc*np.imag(omegam-omegap)
    Hp[6,3] = 0.5*cc*(hm-hp-np.real(omegam)+np.real(omegap))
    Hp[6,4] = -0.5*1j*ca*cc*np.imag(omegap+omegam)
    Hp[6,5] = -0.5*1j*cb*cc*np.imag(omegap+omegam)
    Hp[6,6] = 0.5*( cc*cc*(hp+hm) + 2*(h0-omega0)*sc*sc - cc*cc*np.real(omegap+omegam) )
    Hp[6,7] = 0.5*(cc*cd*(hp+hm) + 2*(h0-omega0)*sc*sd - cc*cd*np.real(omegap+omegam) )

    for i in range(8):
        Hp[i,4] = np.conjugate(Hp[4,i])
        Hp[i,6] = np.conjugate(Hp[6,i])


    ### The Hamiltonian
    H[0,0] = omega + E[4] + Hp[4,4]
    H[0,1] = Hp[4,6]
    H[1,0] = Hp[6,4] 
    H[1,1] = omega + E[6] + Hp[6,6]

    ### The second order 
    for l in [0,1,2,3,5,7]:
        H[0,0] = H[0,0] + np.abs(Hp[4,l])**2 / (E[4]-E[l])
        H[0,1] = H[0,1] + Hp[4,l]*Hp[l,6] / (E[4]-E[l])
        H[1,0] = H[1,0] + Hp[6,l]*Hp[l,4] / (E[6]-E[l])
        H[1,1] = H[1,1] + np.abs(Hp[6,l])**2 / (E[6]-E[l])

    ### The third order 
    #for l in [0,1,2,3,5,7]:
    #    for m in [0,1,2,3,5,7]:
    #        H[0,0] = H[0,0] - 0.5*Hp[4,l]*Hp[l,m]*Hp[m,4] / ((E[4]-E[l])*(E[m]-E[l])) 
    #        H[0,0] = H[0,0] - 0.5*Hp[4,m]*Hp[m,l]*Hp[l,4] / ((E[4]-E[l])*(E[m]-E[l])) 

    return H 


##### ======================================================================
##### The parameters of the effective model 
omega = 0.29780940 
eta = 0.0032
v = 0.317 
U = -0.01503117
W = 0.001466 
alpha = -0.062

pomega = 0.0
omega1 = omega*(1+pomega)
omega2 = omega*(1-pomega)

v1 = v
v2 = v 

pU = 0
U1 = U*(1+pU)
U2 = U*(1-pU)

pW = 0
W1 = W*(1+pW)
W2 = W*(1-pW)

eta1 = eta 
eta2 = eta 
alpha1 = alpha 
alpha2 = alpha 

d0 = 0.35 
dist = 0.1 
V = 0.038
beta = -0.3

##### ===================================================================================
##### The arrays of k and q 
### The array of intrinsic momenta k
Nk = 201 
Kmax = 0.05
k_array = np.linspace(0,Kmax,Nk)

### The array of synthetic momenta q
Nq = 201 
Qmax = 0.1
q_array = np.linspace(0,Qmax,Nq)

##### ===================================================================================
### The array of energy 
Energy_array = np.zeros((Nk,Nq,2))

### We scan over the intrinsic and the synthetic momenta
for ik in range(Nk):
    # The genuine momentum
    k = k_array[ik]

    for iq in range(Nq):
        # The synthetic momentum
        q = q_array[iq]

        # The Hamiltonian
        H = TwoBandModel(omega,eta,v,U,W,alpha,d0,dist,V,beta,k,q)

        # Diagonalize the Hamiltonian
        evals,evecs = sla.eigh(H)

        # We save the energy eigenvalues to the array Energy_array 
        Energy_array[ik,iq,:] = evals 

##### ===================================================================================
### Plot the energy bands 

### The array of X and Y
X,Y = np.meshgrid(k_array,q_array)

### Figure 
fig,ax = plt.subplots(subplot_kw={'projection':'3d'},figsize=(12,10))
ax.plot_surface(X,Y,Energy_array[:,:,0].T)
ax.plot_surface(X,Y,Energy_array[:,:,1].T)
ax.set_xlabel('k',fontsize=14)
ax.set_ylabel('q',fontsize=14)
plt.show()