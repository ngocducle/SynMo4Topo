program Dispersion2DSlab1L
implicit none 

real(8),parameter :: math_pi = 3.141592653589793d0

integer(8) :: Nk,i  

real(8) :: Kmax,omegaM,v,U,W,alpha  
real(8) :: kx,ky 

complex(8), dimension(4,4) :: Hamiltonian,States 
real(8),dimension(:),allocatable :: kx_array,ky_array 

!!! ==========================================================================
!!! The variables for LAPACK routine
!!! ATTENTION! Please never change the type of the variables
integer :: info !,lwork,lapackdim 
double precision, dimension(12) :: rwork 
complex(8), dimension(128) :: work 
real(8),dimension(4) :: eigenvalues 

!lapackdim = 2 
!lwork = 128 

omegaM = 0.297 
v = 0.317 
U = -0.0154
W = 0.0014 
alpha = 0.1 

! Number of k-points to interpolate 
Nk = 39 

! Maximum value of k 
Kmax = 0.1 

allocate(kx_array(2*Nk+1))
allocate(ky_array(2*Nk+1))

do i = 1,Nk+1  
    kx_array(i) = -Kmax*(i-Nk-1)/(Nk+1) 
    ky_array(i) = -kx_array(i) 
    print*, kx_array(i),ky_array(i)
end do ! i-loop 

do i = Nk+2,2*Nk+2 
    kx_array(i) = Kmax*(i-Nk-1)/(Nk+1)
    ky_array(i) = kx_array(i)
    print*, kx_array(i),ky_array(i) 
end do ! i-loop 

do i = 1,2*Nk+2 
    kx = kx_array(i) 
    ky = ky_array(i) 

    Hamiltonian(1,1) = omegaM + v*(kx+ky)/sqrt(2.0) 
    Hamiltonian(1,2) = W 
    Hamiltonian(1,3) = W 
    Hamiltonian(1,4) = U*(1+alpha) 

    Hamiltonian(2,1) = W 
    Hamiltonian(2,2) = omegaM + v*(kx-ky)/sqrt(2.0) 
    Hamiltonian(2,3) = U*(1-alpha)
    Hamiltonian(2,4) = W 

    Hamiltonian(3,1) = W 
    Hamiltonian(3,2) = U*(1-alpha)
    Hamiltonian(3,3) = omegaM + v*(-kx+ky)/sqrt(2.0)
    Hamiltonian(3,4) = W 

    Hamiltonian(4,1) = U*(1+alpha)
    Hamiltonian(4,2) = W 
    Hamiltonian(4,3) = W 
    Hamiltonian(4,4) = omegaM + v*(-kx-ky)/sqrt(2.0)  

    States = Hamiltonian
    call zheev('V','U',4,States,2,eigenvalues,work,128,rwork,info) 
end do ! i-loop 

end program Dispersion2DSlab1L 