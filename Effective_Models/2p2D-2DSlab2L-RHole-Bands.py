import numpy as np 
import scipy 
import scipy.linalg as sla 
import cmath 
import matplotlib.pyplot as plt 
from matplotlib import cm,colors 

import sys 
sys.path.insert(0,'../src')
from EffectiveModels import EffHam_2DSlab2L_RHole_M 

##### The MAIN program
def main():
    ### The arrays of genuine momenta 
    Nk = 201 
    kx_array = np.linspace(0.0,0.1,Nk)
    ky_array = np.linspace(0.0,0.1,Nk)
    kd_array = np.linspace(0.0,0.1,Nk) # kx = ky = kd
    ka_array = np.linspace(0.0,0.1,Nk) # kx = -ky = ka 

    ### The arrays of synthetic momenta 
    Nq = 201 
    qx_array = np.linspace(0.0,0.5,Nq)
    qy_array = np.linspace(0.0,0.5,Nq)
    qd_array = np.linspace(0.0,0.5,Nq) # qx = qy = qd 
    qa_array = np.linspace(0.0,0.5,Nq) # qx = -qy = qa 

    ### The parameters 
    omega1 = 0.29780940
    v1 = 0.317
    U1 = -0.01537
    W1 = 0.001466
    alpha1 = 0.1

    omega2 = 0.29780940
    v2 = 0.317
    U2 = -0.01537
    W2 = 0.001466
    alpha2 = 0.1

    V = 0.025

    delta = 0.05

    ##### ===============================================================================
    ### Define the shift of the plane passing pass (0,0,0.5,0.5) and 
    ### spanned by the vectors (1,-1,0,0) and (0,0,1,1)
    ### The plane has the normal vectors:
    ###         n1 = (1,1,0,0)
    ###         n2 = (0,0,1,-1)
    ### We shift the M (0,0,0.5,0.5) point by a vector n = a*n1 + b*n2 
    ### to the point M'(a,a,0.5+b,0.5-b)
    ### The plane passing through M' and parallel to the initial plane 
    ### has the following parameterization:
    ###
    ###         kx = a + ka 
    ###         ky = a - ka 
    ###         qx = 0.5 + b + qd 
    ###         qy = 0.5 + b + qd 
    a = 0.05
    b = 0.15 

    ##### ===============================================================================
    ##### Case kx=ky=kd, qx=qy=qd 
    Band1 = np.zeros((Nk,Nq))
    Band2 = np.zeros((Nk,Nq))
    Band3 = np.zeros((Nk,Nq))
    Band4 = np.zeros((Nk,Nq))

    for ik in range(Nk):
        for iq in range(Nq):
            # The momenta 
            kd = kd_array[ik]
            qd = qd_array[iq]

            # The Hamiltonian 
            Hamiltonian = EffHam_2DSlab2L_RHole_M(kd,kd,qd+delta,qd-delta,
                                                  omega1,v1,U1,W1,alpha1,
                                                  omega2,v2,U2,W2,alpha2,V)

            # Diagonalize the Hamiltonian to get the eigenvalues and eigenstates
            evalues,evectors = sla.eigh(Hamiltonian)

            # Save the values to bands 
            Band1[ik,iq] = evalues[0]
            Band2[ik,iq] = evalues[1]
            Band3[ik,iq] = evalues[2]
            Band4[ik,iq] = evalues[3] 


    ### Plot the band structure for kx = ky = kd, qx = qy = qd 
    namesave = "kd-qd-alpha_{0:.4f}".format(alpha1)+"-delta_{0:.4f}".format(delta)

    X,Y = np.meshgrid(kd_array,qd_array)

    fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                          figsize = (12,10))
    ax.plot_surface(X,Y,Band1.T,antialiased=False,cmap='plasma')
    ax.plot_surface(X,Y,Band2.T,antialiased=False,cmap='plasma')
    #ax.plot_surface(X,Y,Band3.T,antialiased=False,cmap='plasma')
    #ax.plot_surface(X,Y,Band4.T,antialiased=False,cmap='plasma')
    ax.set_xlabel(r'$k_d$',fontsize=14)
    ax.set_ylabel(r'$q_d$',fontsize=14)
    ax.set_title(r'$k_x = k_y, q_x = q_y, \alpha = $'+str(alpha1)+
                 ', a = '+str(a)+', b = '+str(b),fontsize=14)
    ax.view_init(elev=5,azim=135,roll=0)
    plt.savefig(namesave+'.png')
    plt.show()
    plt.close()

    ##### ===============================================================================
    ##### Case kx=-ky=ka, qx=qy=qd 
    Band1 = np.zeros((Nk,Nq))
    Band2 = np.zeros((Nk,Nq))
    Band3 = np.zeros((Nk,Nq))
    Band4 = np.zeros((Nk,Nq))

    for ik in range(Nk):
        for iq in range(Nq):
            # The momenta 
            ka = ka_array[ik]
            qd = qd_array[iq]

            # The Hamiltonian 
            Hamiltonian = EffHam_2DSlab2L_RHole_M(ka+a,-ka+a,qd+b,qd+b,
                                                  omega1,v1,U1,W1,alpha1,
                                                  omega2,v2,U2,W2,alpha2,V)

            # Diagonalize the Hamiltonian to get the eigenvalues and eigenstates
            evalues,evectors = sla.eigh(Hamiltonian)

            # Save the values to bands 
            Band1[ik,iq] = evalues[0]
            Band2[ik,iq] = evalues[1]
            Band3[ik,iq] = evalues[2]
            Band4[ik,iq] = evalues[3] 


    ### Plot the band structure for kx = -ky = ka, qx = qy = qd 
    namesave = "ka-qd-alpha_{0:.4f}".format(alpha1)+"-delta_{0:.4f}".format(delta)

    X,Y = np.meshgrid(ka_array,qd_array)

    fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                          figsize = (12,10))
    ax.plot_surface(X,Y,Band1.T,antialiased=False,cmap='plasma')
    ax.plot_surface(X,Y,Band2.T,antialiased=False,cmap='plasma')
    #ax.plot_surface(X,Y,Band3.T,antialiased=False,cmap='plasma')
    #ax.plot_surface(X,Y,Band4.T,antialiased=False,cmap='plasma')
    ax.set_xlabel(r'$k_a$',fontsize=14)
    ax.set_ylabel(r'$q_d$',fontsize=14)
    ax.set_title(r'$k_x = -k_y, q_x = q_y, \alpha = $'+str(alpha1)+
                 ', a = '+str(a)+', b = '+str(b),fontsize=14)
    ax.view_init(elev=5,azim=135,roll=0)
    plt.savefig(namesave+'.png')
    plt.show()
    plt.close()

    ##### ===============================================================================
    ##### Case kx = -ky = ka, qx = -qy = qa 
    Band1 = np.zeros((Nk,Nq))
    Band2 = np.zeros((Nk,Nq))
    Band3 = np.zeros((Nk,Nq))
    Band4 = np.zeros((Nk,Nq))

    for ik in range(Nk):
        for iq in range(Nq):
            # The momenta 
            ka = ka_array[ik]
            qa = qa_array[iq]

            # The Hamiltonian 
            Hamiltonian = EffHam_2DSlab2L_RHole_M(ka,-ka,qa,-qa,
                                                  omega1,v1,U1,W1,alpha1,
                                                  omega2,v2,U2,W2,alpha2,V)

            # Diagonalize the Hamiltonian to get the eigenvalues and eigenstates
            evalues,evectors = sla.eigh(Hamiltonian)

            # Save the values to bands 
            Band1[ik,iq] = evalues[0]
            Band2[ik,iq] = evalues[1]
            Band3[ik,iq] = evalues[2]
            Band4[ik,iq] = evalues[3] 


    ### Plot the band structure for kx = -ky = ka, qx = -qy = qa 
    namesave = "ka-qa-alpha_{0:.4f}".format(alpha1)+"-delta_{0:.4f}".format(delta)

    X,Y = np.meshgrid(ka_array,qa_array)

    fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                          figsize = (12,10))
    ax.plot_surface(X,Y,Band1.T,antialiased=False,cmap='plasma')
    ax.plot_surface(X,Y,Band2.T,antialiased=False,cmap='plasma')
    ax.plot_surface(X,Y,Band3.T,antialiased=False,cmap='plasma')
    ax.plot_surface(X,Y,Band4.T,antialiased=False,cmap='plasma')
    ax.set_xlabel(r'$k_a$',fontsize=14)
    ax.set_ylabel(r'$q_a$',fontsize=14)
    ax.set_title(r'$k_x = -k_y, q_x = -q_y, \alpha = $'+str(alpha1),fontsize=14)
    ax.view_init(elev=5,azim=135,roll=0)
    plt.savefig(namesave+'.png')
    #plt.show()
    plt.close()

    ##### ===============================================================================
    ##### Case kx = ky = kd, qx = -qy = qa 
    Band1 = np.zeros((Nk,Nq))
    Band2 = np.zeros((Nk,Nq))
    Band3 = np.zeros((Nk,Nq))
    Band4 = np.zeros((Nk,Nq))

    for ik in range(Nk):
        for iq in range(Nq):
            # The momenta 
            kd = kd_array[ik]
            qa = qa_array[iq]

            # The Hamiltonian 
            Hamiltonian = EffHam_2DSlab2L_RHole_M(kd,kd,qa,-qa,
                                                  omega1,v1,U1,W1,alpha1,
                                                  omega2,v2,U2,W2,alpha2,V)

            # Diagonalize the Hamiltonian to get the eigenvalues and eigenstates
            evalues,evectors = sla.eigh(Hamiltonian)

            # Save the values to bands 
            Band1[ik,iq] = evalues[0]
            Band2[ik,iq] = evalues[1]
            Band3[ik,iq] = evalues[2]
            Band4[ik,iq] = evalues[3] 


    ### Plot the band structure for kx = ky = kd, qx = -qy = qa 
    namesave = "kd-qa-alpha_{0:.4f}".format(alpha1)+"-delta_{0:.4f}".format(delta)

    X,Y = np.meshgrid(kd_array,qa_array)

    fig,ax = plt.subplots(subplot_kw = {'projection':'3d'},
                          figsize = (12,10))
    ax.plot_surface(X,Y,Band1.T,antialiased=False,cmap='plasma')
    ax.plot_surface(X,Y,Band2.T,antialiased=False,cmap='plasma')
    ax.plot_surface(X,Y,Band3.T,antialiased=False,cmap='plasma')
    ax.plot_surface(X,Y,Band4.T,antialiased=False,cmap='plasma')
    ax.set_xlabel(r'$kd$',fontsize=14)
    ax.set_ylabel(r'$q_a$',fontsize=14)
    ax.set_title(r'$k_x = k_y, q_x = -q_y, \alpha = $'+str(alpha1),fontsize=14)
    ax.view_init(elev=5,azim=135,roll=0)
    plt.savefig(namesave+'.png')
    #plt.show()
    plt.close()

##### Run the MAIN program
if __name__ == '__main__':
    main()