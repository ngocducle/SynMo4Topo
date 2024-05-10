program ChernNumber 
implicit none 

integer(8) :: Nx,Ny
integer(8) :: i,j  
real(8) :: velocity,U,V,Delta 
real(8) :: x,y 

real(8),dimension(:), allocatable :: Kx,Ky 
real(8),dimension(:,:), allocatable :: Energy1,Energy2

complex(8),dimension(2) :: State1,State2 
complex(8),dimension(2,2) :: Hamiltonian,dHx,dHy 
complex(8),dimension(:,:), allocatable :: Berry_curvature1,Berry_curvature2

!!! ===============================================================================
!!! The variables for LAPACK routine 
!!! ATTENTION! Please never change the type of the variables 
integer :: info,lwork,lapackdim  
double precision, dimension(:),allocatable :: rwork 
complex(8), dimension(:),allocatable :: work 
real(8), dimension(:),allocatable :: eigenvalues

velocity = 0.5
print*,'# velocity = ',velocity 

U = 0.02
print*,'# U = ', U 

V = 1.5*U 
print*,'# V = ', V 

Delta = 0.1*U 
print*,'# Delta = ', Delta 

dHy(1,1) = 0.0
dHy(1,2) = cmplx(0.0,-V,8)
dHy(2,1) = -dHy(1,2)
dHy(2,2) = 0.0 

Nx = 2001
Ny = 2001  

allocate(Kx(Nx))
allocate(Ky(Ny))
allocate(Energy1(Nx,Ny))
allocate(Energy2(Nx,Ny))
allocate(Berry_curvature1(Nx,Ny))
allocate(Berry_curvature2(Nx,Ny))

do i = 1,Nx 
    Kx(i) = -0.5+(i-1)*1.0/(Nx-1) 
    !print*, Kx(i) 
end do ! i-loop 

do j = 1,Ny 
    Ky(j) = -0.5+(j-1)*1.0/(Ny-1) 
    !print*, Ky(j)
end do ! j-loop 

lapackdim = 2 
lwork = 128 
allocate(work(lwork))
allocate(rwork(6))
allocate(eigenvalues(2))

do j = 1,Ny 
do i = 1,Nx 
    !x = Kx(i)
    !print*, 'x = ', x 

    !y = Ky(j) 

    !print*,U - V + 0.5*velocity**2 * x**2/U

    Hamiltonian(1,1) = cmplx(U - V + 0.5*velocity**2 * Kx(i)**2/U,0.0d0,8)  
    Hamiltonian(1,2) = cmplx(-velocity*Delta/U * Kx(i), -V*Ky(j),8) 
    Hamiltonian(2,1) = cmplx(-velocity*Delta/U * Kx(i),  V*Ky(j),8)
    Hamiltonian(2,2) = -Hamiltonian(1,1)

    !print*,Hamiltonian 

    call zheev('V','U',lapackdim,Hamiltonian,2,eigenvalues,work,128,rwork,info)
    !print*,'# LAPACK info',info 
    !print*, eigenvalues(1),eigenvalues(2)

    Energy1(i,j) = eigenvalues(1)
    Energy2(i,j) = eigenvalues(2) 

    State1(1) = Hamiltonian(1,1)
    State1(2) = Hamiltonian(2,1)
    State2(1) = Hamiltonian(1,2)
    State2(2) = Hamiltonian(2,2) 

    !print*,State1(1),State1(2)
    !print*,State2(1),State2(2)

    dHx(1,1) = velocity**2 * Kx(i)/U 
    dHx(1,2) = -velocity*Delta/U 
    dHx(2,1) = -velocity*Delta/U 
    dHx(2,2) = -dHx(1,1) 

    Berry_curvature1(i,j) = -2.0*aimag(matmul(State1,matmul(dHx,State2)))
end do ! i-loop 
end do ! j-loop 

deallocate(work)
deallocate(rwork)
deallocate(eigenvalues)

end program ChernNumber 