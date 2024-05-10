#include <iostream>  

int main()
{
    float v = 0.5;
    std::cout << "v = " << v << "\n"; 

    float U = 0.02;
    std::cout << "U = " << U << "\n";

    float V = 1.5*U;
    std::cout << "V = " << V << "\n"; 

    float Delta = 0.1*U; 
    std::cout << "Delta = " << Delta << "\n";

    int Nx = 201;
    int Ny = 201;

    float Kx[Nx];
    for (int i = 0; i<Nx; i++)
    {
        Kx[i] = -0.5+float(i)*1.0/float(Nx-1); 
        /* std::cout << Kx[i] << "\n"; */
    }

    float Ky[Ny];
    for (int i = 0; i<Ny;i++)
    {
        Ky[i] = -0.5+float(i)*1.0/float(Ny-1); 
    }

    float Energy[Nx][Ny][2];
    float F_array[Nx][Ny][2]; 

    float Hamiltonian[2][2];
    float dHx[2][2];
    float dHy[2][2];

    
     

    return 0; 
}