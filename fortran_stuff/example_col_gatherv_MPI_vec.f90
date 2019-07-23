program main
  !
  ! TODO:
  !   1. Make sure using offset and block counts with mpi_type_vector
  !     works correctly. 
  !   2. See if we can exploit Column-major ordering and contiguous memory 
  !     in Fortran to easily count and send the columns to be stacked on the
  !     receiving end. 
  !
  !
  use mpi

  implicit none
  integer :: ierr, procno, numprocs, root, comm, vec_type !! MPI values
  integer :: iProcs !! loop iterators

  real(kind=8), allocatable :: vec(:), Gvec(:)
  integer, allocatable :: displs(:), rcounts(:)

!======================================

  call mpi_init(ierr)
  comm = mpi_comm_world
  call mpi_comm_rank(comm, procno, ierr)
  call mpi_comm_size(comm, numprocs, ierr)

  root = 0

!======================================

  call mpi_type_vector(4, 4, 4, mpi_real8, vec_type, ierr)
  call mpi_type_commit(vec_type, ierr)

!======================================

  allocate(vec(8))
  vec(:) = real(procno, 8)

  !print *, vec

  if (procno .eq. root) then
    allocate(Gvec(numprocs * 8))
    allocate(displs(numprocs))
    allocate(rcounts(numprocs))

    rcounts = 2
    !print *, rcounts

    ! calculate displacements
    offset = 0
    do iProcs = 1,nProcs

    enddo
  endif 

!======================================

  !call mpi_gatherv(vec, 2, vec_type)

!======================================
  if (allocated(vec)) deallocate(vec)
  if (allocated(Gvec)) deallocate(Gvec)
  if (allocated(displs)) deallocate(displs)
  if (allocated(rcounts)) deallocate(rcounts)

  call mpi_type_free(vec_type, ierr)
  call mpi_finalize(ierr)
end program main