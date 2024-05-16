#include <iostream>  
#include <complex>

using namespace std;

// zheev_ is a symbol in the LAPACK library files 
extern "C" {
    extern void zheev_(char* jobz,char* uplo,int* lapackdim,complex<double>* A,
                       int* lda, double* w, complex<double>* work,
                       int* lwork,double* rwork,int* info);
}

int main()
{
    // The velocity
    float v = 0.5;
    std::cout << "v = " << v << "\n"; 

    // The parameter U 
    float U = 0.02;
    std::cout << "U = " << U << "\n";

    // The parameter V 
    float V = 1.5*U;
    std::cout << "V = " << V << "\n"; 
    //std::cout << "v*V = " << v*V << "\n";

    // The parameter Delta 
    float Delta = 0.1*U; 
    std::cout << "Delta = " << Delta << "\n";

    // Definition of the imaginary unit 
    complex <float> ic1(0.0,1.0);

    // The array of kx 
    int Nx = 201;

    float kx_array[Nx];
    for (int i = 0; i<Nx; i++)
    {
        kx_array[i] = -0.5+float(i)*1.0/float(Nx-1); 
        /* std::cout << Kx[i] << "\n"; */
    }

    // The array of ky 
    int Ny = 201;

    float ky_array[Ny];
    for (int i = 0; i<Ny;i++)
    {
        ky_array[i] = -0.5+float(i)*1.0/float(Ny-1); 
    }

    // Allocate the arrays of energy and Berry curvature 
    float E[Nx][Ny][2];
    float F_array[Nx][Ny][2]; 

    // Allocate the Hamiltonian and its derivatives with respect to x and y 
    complex <double> Hamiltonian[2][2];
    complex <double> dHx[2][2];
    complex <double> dHy[2][2];

    for (int i = 0; i<Nx; i++)
    {
        for (int j = 0; j<Ny; j++)
        {
            float kx = kx_array[i];
            float ky = ky_array[j];

            // The Hamiltonian H(kx,ky)
            Hamiltonian[0][0] = U-V+0.5*pow(v,2);
            Hamiltonian[0][1] = -v*Delta*kx/U-V*ky*ic1; 
            Hamiltonian[1][0] = -v*Delta*kx/U+V*ky*ic1;
            Hamiltonian[1][1] = -Hamiltonian[0][0];

            /*cout << Hamiltonian[0][0] << ' ' << Hamiltonian[0][1] << '\n';
            cout << Hamiltonian[1][0] << ' ' << Hamiltonian[1][1] << '\n';
            cout << '\n';*/

            // The derivative of the Hamiltonian with respect to x 
            dHx[0][0] = pow(v,2)*kx/U;
            dHx[0][1] = -v*Delta/U;
            dHx[1][0] = -v*Delta/U;
            dHx[1][1] = -dHx[0][0];

            // The derivative of the Hamiltonian with respect to y 
            dHy[0][0] = 0.0;
            dHy[0][1] = -V*ic1;
            dHy[1][0] = V*ic1;
            dHy[1][1] = 0.0;  

            // LAPACK diagonalization 
            complex <double> *States[2][2];
            
            States[0][0] = Hamiltonian[0][0];
            States[0][1] = Hamiltonian[0][1];
            States[1][0] = Hamiltonian[1][0];
            States[1][1] = Hamiltonian[1][1];

            int lapackdim = 2,lda = 2,lwork=128;
            int info; 
            double w[2],rwork[6];
            complex <double> work[lwork];

            zheev_("V","U",&lapackdim,&States,&lda,w,&work,&lwork,rwork,&info);

            // Check for convergence 
            if (info>0){
                cout << "Error \n";
            }
        }
    }
     

    return 0; 
}