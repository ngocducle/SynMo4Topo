module matrix_diagonalization

implicit none 

contains 

! =========================================================================
! Lapack_diagonalization: Diagonalization using LAPACK library 
subroutine Lapack_diagonalization(HamiltonianDim,Hamiltonian,MatrixEigenvalues,MatrixEigenvectors)
    
integer(8) :: HamiltonianDim 
complex(8),dimension(:) :: MatrixEigenvalues
complex(8),dimension(:,:) :: Hamiltonian,MatrixEigenvectors 

!!! The variables for LAPACK routine
!!! ATTENTION! Please never change the type of the variables! 
integer :: lapackdim,lwork,info
real, dimension(:), allocatable :: rwork
complex, dimension(:), allocatable :: work  

real, dimension(:), pointer :: eigenvalues
complex, dimension(:,:), pointer :: Matrix

print*,'LAPACK diagonalization'

allocate(Matrix(lapackdim,lapackdim))
allocate(eigenvalues(lapackdim))

Matrix = Hamiltonian

! This paragraph is used for cheev
lwork = 64*lapackdim 
allocate(work(lwork))
allocate(rwork(3*lapackdim))

call cheev('N','L',lapackdim,Matrix,lapackdim,eigenvalues,work,lwork,rwork,info)

print*,'# LAPACK info',info

! Deallocate the arrays 
deallocate(Matrix)
deallocate(eigenvalues)
deallocate(rwork)
deallocate(work)

end subroutine Lapack_diagonalization

end module matrix_diagonalization