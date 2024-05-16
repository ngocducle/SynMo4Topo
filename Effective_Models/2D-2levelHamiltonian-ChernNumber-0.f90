program ChernNumber 

!!!!! ============================================================================== !!!!! 
!!!                                                                                    !!!
!!!         Calculate the Chern number for the Hamiltonian:                            !!!
!!!                                                                                    !!!
!!!      H = [ U - V + (v^2/2U)*kx^2, -(v*Delta/U)*kx - i*V*ky;                        !!! 
!!!           -(v*Delta/U)*kx + i*V*ky, -U + V - (v^2/2U)*kx^2 ]                       !!!
!!!                                                                                    !!!
!!!!! ============================================================================== !!!!! 

!!! To compile:
!   gfortran -llapack 1p1D-2levelHamiltonian-ChernNumber.f90

implicit none 

real(8),parameter :: math_pi = 3.141592653589793d0

integer(8) :: Nx,Ny
integer(8) :: i,j,k,l,m,n 
real(8) :: velocity,U,V,Delta 
real(8) :: kx,ky,kx_min,kx_max,ky_min,ky_max,dkx,dky
!real(8) :: Chern1,Chern2 
complex(8) :: F1,F2

real(8),dimension(2) :: Chern 
real(8),dimension(:), allocatable :: kx_scan,ky_scan 
real(8),dimension(:,:), allocatable :: Energy1,Energy2

complex(8),dimension(2) :: State1,State2 
complex(8),dimension(2,2) :: Hamiltonian,A,dHx,dHy,dHxe,dHye,State,StateH 
!complex(8),dimension(:,:), allocatable :: Berry_curvature1,Berry_curvature2
complex(8),dimension(:,:,:), allocatable :: F_array 

!!! ===============================================================================
!!! The variables for LAPACK routine 
!!! ATTENTION! Please never change the type of the variables 
integer :: info,lwork,lapackdim  
double precision, dimension(:),allocatable :: rwork 
complex(8), dimension(:),allocatable :: work 
real(8), dimension(:),allocatable :: eigenvalues

!!!!! The parameters of the model 
velocity = 0.5
print*,'# velocity = ',velocity 

U = 0.02
print*,'# U = ', U 

V = 1.5*U 
print*,'# V = ', V 

Delta = 0.1*U 
print*,'# Delta = ', Delta 

!!!!! The matrix dH/dy 
dHy(1,1) = 0.0
dHy(1,2) = cmplx(0.0,-V,8)
dHy(2,1) = -dHy(1,2)
dHy(2,2) = 0.0 

!!!!!  Allocate the arrays of momenta,energy and Berry curvature
Nx = 2001
Ny = 2001 
allocate(kx_scan(Nx))
allocate(ky_scan(Ny))
allocate(Energy1(Nx,Ny))
allocate(Energy2(Nx,Ny))
!allocate(Berry_curvature1(Nx,Ny))
!allocate(Berry_curvature2(Nx,Ny))
allocate(F_array(Nx,Ny,2))

!!!!! The array of kx 
kx_min = -3.0d0 
kx_max = 3.0d0  
dkx = (kx_max-kx_min)/(Nx-1)
do i = 1,Nx 
    kx_scan(i) = kx_min+(i-1)*dkx  
    !print*, kx_scan(i) 
end do ! i-loop 

!!!!! The array of ky 
ky_min = -30.0d0 
ky_max = 30.0d0 
dky = (ky_max-ky_min)/(Ny-1) 
do j = 1,Ny 
    ky_scan(j) = -0.5+(j-1)*dky  
    !print*, ky_scan(j)
end do ! j-loop 

!!!!! Initiate the Chern numbers of  band 1 and band 2 
!Chern1 = 0.0d0 
!Chern2 = 0.0d0 
Chern(1) = 0.0d0 
Chern(2) = 0.0d0 

!!!!! LAPACK allocations 
lapackdim = 2 
lwork = 128 
allocate(work(lwork))
allocate(rwork(6))
allocate(eigenvalues(2))

!!!!! We scan over the momenta kx and ky 
do j = 1,Ny 
do i = 1,Nx 
    !!! Take the values of kx and ky 
    kx = kx_scan(i)
    ky = ky_scan(j) 

    !!! The Hamiltonian H(kx,ky)
    Hamiltonian(1,1) = cmplx(U - V + 0.5*velocity**2 * kx**2/U,0.0d0,8)  
    Hamiltonian(1,2) = cmplx(-velocity*Delta/U * kx, -V*ky,8) 
    Hamiltonian(2,1) = cmplx(-velocity*Delta/U * kx,  V*ky,8)
    Hamiltonian(2,2) = -Hamiltonian(1,1)

    !print*,'Hamiltonian = '
    !print*,Hamiltonian(1,:)
    !print*,Hamiltonian(2,:)

    A = Hamiltonian 

    !!! Diagonalization 
    call zheev('V','U',lapackdim,A,2,eigenvalues,work,128,rwork,info)
    !print*,'# LAPACK info',info 
    !print*, eigenvalues(1),eigenvalues(2)

    !!! The energies of the bands 
    Energy1(i,j) = eigenvalues(1)
    Energy2(i,j) = eigenvalues(2) 

    !!! The eigenstates 
    !State1(1) = Hamiltonian(1,1)
    !State1(2) = Hamiltonian(2,1)
    !State2(1) = Hamiltonian(1,2)
    !State2(2) = Hamiltonian(2,2) 

    !print*,State1(1),State1(2)
    !print*,State2(1),State2(2)

    !!! The matrix of energy eigenstates
    State = Hamiltonian 

    !print*, 'Hamiltonian = ' 
    !print*, Hamiltonian(1,:)
    !print*, Hamiltonian(2,:)

    !print*, 'State = ' 
    !print*, State(1,:)
    !print*, State(2,:) 

    StateH = transpose(State)
    do l=1,2
    do k=1,2 
        StateH(k,l) = conjg(StateH(k,l))
    end do ! k-loop 
    end do ! l-loop 
    !print*, 'StateH = '
    !print*, StateH(1,:)
    !print*, StateH(2,:)

    !!! The matrix dH/dx 
    dHx(1,1) = velocity**2 * kx/U 
    dHx(1,2) = -velocity*Delta/U 
    dHx(2,1) = -velocity*Delta/U 
    dHx(2,2) = -dHx(1,1) 

    !!! The matrix dHxe represents the operator dHx in the basis of energy eigenstates
    dHxe = matmul(StateH,matmul(dHx,State)) 

    !!! The matrix dHdeltae represents the operator dHdelta in the basis of energy eigenstates
    dHye = matmul(StateH,matmul(dHy,State))

    !!! Berry curvature and Chern number of band 1 
    !F1 = cmplx(0.0,0.0)
    !F2 = cmplx(0.0,0.0)

    !do l = 1,2 
    !do k = 1,2 
    !    F1 = F1 + conjg(State1(k))*dHx(k,l)*State2(l)
    !    F2 = F2 + conjg(State2(k))*dHy(k,l)*State1(l)
    !end do ! k-loop 
    !end do ! l-loop 

    !Berry_curvature1(i,j) = Berry_curvature1(i,j) &
    !    & -2*aimag(F1*F2)/(eigenvalues(1)-eigenvalues(2))**2

    !Chern1 = Chern1 + Berry_curvature1(i,j)

    !!! Berry curvature and Chern number of band 2 
    !F1 = cmplx(0.0,0.0)
    !F2 = cmplx(0.0,0.0)

    !do l = 1,2 
    !do k = 1,2 
    !    F1 = F1 + conjg(State2(k))*dHx(k,l)*State1(l)
    !    F2 = F2 + conjg(State1(k))*dHy(k,l)*State2(l) 
    !end do ! k-loop 
    !end do ! l-loop 

    !Berry_curvature2(i,j) = Berry_curvature2(i,j) &
    !    & -2*aimag(F1*F2)/(eigenvalues(1)-eigenvalues(2))**2

    !Chern2 = Chern2 + Berry_curvature2(i,j)

    do n=1,2
        F_array(i,j,n) = cmplx(0.0,0.0,8) 
    
        do m=1,2 
            if (m /= n) then 
                F_array(i,j,n) = F_array(i,j,n) & 
                    &-2.0*aimag(dHxe(n,m)*dHye(m,n))/(eigenvalues(n)-eigenvalues(m))**2
            end if ! m /= n -IF 
        end do ! m-loop

        Chern(n) = Chern(n) + F_array(i,j,n) 
    end do ! n-loop 

end do ! i-loop 
end do ! j-loop 

!!!!! Chern number of band 1 
Chern(1) = Chern(1)*dkx*dky/(2.0*math_pi)
print*,'# Chern1 = ', Chern(1) 

!!!!! Chern number of band 2 
Chern(2) = Chern(2)*dkx*dky/(2.0*math_pi)
print*,'# Chern2 = ', Chern(2) 

deallocate(work)
deallocate(rwork)
deallocate(eigenvalues)
deallocate(F_array)

end program ChernNumber 