program ChernNumber 
implicit none 

real(8),parameter :: math_pi = 3.141592653589793d0 

integer(8) :: i,j,m,n 
integer(8) :: Nx,Ny 

real(8) :: velocity,U,V,Delta 
real(8) :: kx_max,dkx,ky_max,dky,kx,ky,val 
real(8) :: C1,C2   

real(8), dimension(:), allocatable :: kx_array,ky_array 
real(8), dimension(:,:,:), allocatable :: F_array,E ! array of Berry curvature and energy 

complex(8),dimension(2,2) :: Hamiltonian,dHx,dHy,dHxe,dHye 
complex(8),dimension(2,2) :: States,StatesH
!complex(8),dimension(2,2) :: ProdL,ProdR 
!complex(8),dimension(:,:,:), allocatable :: F_array ! array of Berry curvature 

!!! ==========================================================================
!!! The variables for LAPACK routine
!!! ATTENTION! Please never change the type of the variables
integer :: info !,lwork,lapackdim 
double precision, dimension(6) :: rwork 
complex(8), dimension(128) :: work 
real(8),dimension(2) :: eigenvalues 

!lapackdim = 2 
!lwork = 128 

!!! The velocity 
velocity = 0.5
print*,'v = ',velocity 

!!! The parameter U  
U = 0.02 
print*,'U = ',U 

!!! The parameter V
V = 1.5*U 
print*,'V = ',V 

!!! The parameter Delta 
Delta = 0.1*U 
print*,'Delta = ',Delta 

!!! The array of kx 
Nx = 201
kx_max = 0.5d0 
dkx = 2*kx_max/(Nx-1) 

allocate(kx_array(Nx))

do i = 1,Nx 
    kx_array(i) = -kx_max + (i-1)*dkx 
end do ! i-loop 

!!! The array of ky 
Ny = 201 
ky_max = 0.5d0 
dky = 2*ky_max/(Ny-1)

allocate(ky_array(Ny))

do i = 1,Ny 
    ky_array(i) = -ky_max + (i-1)*dky 
end do ! i-loop 

!!! Allocate the arrays of energy and Berry curvature 
allocate(E(Nx,Ny,2))
allocate(F_array(Nx,Ny,2))

!!!!! We scan over the momenta 
do i = 1,Nx 
do j = 1,Ny 

    kx = kx_array(i)
    ky = ky_array(j)

    !!! The Hamiltonian H(kx,ky)
    Hamiltonian(1,1) = cmplx(U-V+0.5*velocity**2*kx**2/U,0.0d0,8)
    Hamiltonian(1,2) = cmplx(-velocity*Delta*kx/U,-V*ky,8)
    Hamiltonian(2,1) = cmplx(-velocity*Delta*kx/U,V*ky,8)
    Hamiltonian(2,2) = -Hamiltonian(1,1)

    !print*,'Hamiltonian = '
    !print*,Hamiltonian(1,:)
    !print*,Hamiltonian(2,:) 

    !!! The derivative of the Hamiltonian with respect to x 
    dHx(1,1) = cmplx(velocity**2 * kx / U, 0.0d0,8) 
    dHx(1,2) = cmplx(-velocity * Delta / U,0.0d0,8)  
    dHx(2,1) = cmplx(-velocity * Delta / U,0.0d0,8) 
    dHx(2,2) = -dHx(1,1)  

    !print*,'dH/dx = '
    !print*,dHx(1,:)
    !print*,dHx(2,:)

    !!! The derivative of the Hamiltonian with respect to y 
    dHy(1,1) = cmplx(0.0d0,0.0d0,8)  
    dHy(1,2) = cmplx(0.0d0,-V,8)
    dHy(2,1) = cmplx(0.0d0,V,8)
    dHy(2,2) = cmplx(0.0d0,0.0d0,8)

    !print*,'dH/dy = '
    !print*,dHy(1,:)
    !print*,dHy(2,:)

    !!! LAPACK diagonalization 
    States = Hamiltonian
    call zheev('V','U',2,States,2,eigenvalues,work,128,rwork,info) 
    !print*,'# LAPACK info', info 

    E(i,j,1) = eigenvalues(1)
    E(i,j,2) = eigenvalues(2)

    !print*,'States = ' 
    !print*,States(1,:)
    !print*,States(2,:) 

    ! The matrix StatesH = States^{\dagger}
    StatesH = transpose(conjg(States))
    !print*,'States^{\dagger| = '
    !print*,StatesH(1,:)
    !print*,StatesH(2,:)

    !!! ATTENTION! The formula to evaluate the Berry curvature is:
    !
    ! F_{xy}^n = \sum_{m \ne n| (-2)*<n|dHx|m><m|dHy|n> / (En-Em)^2 
    !
    ! In fact: <n|dHx|m> and <m|dHy|n> are the matrix elements of 
    ! the operators dHx and dHy in the basis of the energy eigenstates 
    ! of the Hamiltonian 
    ! 
    ! Therefore, we reexpress the matrices dHx and dHy in the basis of 
    ! the eigenstates. The transformation is done by the formula:
    ! 
    !   A' = StatesH*A*States 
    !
    ! here A = dHx or dHy 
    ! and the j-th column of States is the eigenvector corresponding 
    ! to the j-th eigenvalue (j is not the loop index in this program)
    dHxe = matmul(StatesH,matmul(dHx,States))
    dHye = matmul(StatesH,matmul(dHy,States))

    !!! We calculate the Berry curvature 
    do n = 1,2
        F_array(i,j,n) = 0.0d0 

        do m = 1,2 
            if (m /= n) then 
                val = -2.0d0*aimag(dHxe(n,m)*dHye(m,n))/(eigenvalues(n)-eigenvalues(m))**2
                F_array(i,j,n) = F_array(i,j,n) + val 
            end if ! m/=n IF
        end do ! m-loop 
    end do ! n-loop 

end do ! j-loop 
end do ! i-loop 

!!! Calculate  the Chern numbers 
C1 = sum(F_array(:,:,1))*dkx*dky/(2.0*math_pi)
print*,'Chern number C1 = ', C1 

C2 = sum(F_array(:,:,2))*dkx*dky/(2.0*math_pi)
print*,'Chern number C2 = ', C2 

!!!!! Print the data to file
open(unit = 1,file='1p1D-2levelHamiltonian.txt') 
100 format (' ',f12.8,' ',f12.8,' ',f12.8,' ',f12.8)

do n=1,2 
do j=1,Ny 
do i=1,Nx 
    write(unit=1,fmt=100) kx_array(i),ky_array(j),E(i,j,n),F_array(i,j,n)
end do ! i-loop 
end do ! j-loop 
end do ! n-loop 

close(1)

!!! Deallocate the arrays 
deallocate(E)
deallocate(F_array)

end program ChernNumber 