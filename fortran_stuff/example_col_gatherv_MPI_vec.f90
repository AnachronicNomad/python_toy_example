program main
  !
  ! TODO:
  !   1. DONE, blocks and offset size working correctly.
  !   2. ~~See if we can exploit Column-major ordering and contiguous memory 
  !     in Fortran to easily count and send the columns to be stacked on the
  !     receiving end. ~~ <- DONE, it works!!
  !   3. Have arbitrary process ordering (rcount)
  !
  !
  use mpi

  implicit none
  integer :: ierr, procno, numprocs, root, comm, vec_type !! MPI values
  integer :: iProcs, offset, ncols, i, j !! loop iterators, indices
  integer, parameter :: nrows = 5

  real(kind=4), allocatable :: vec(:,:), Gvec(:,:)
  integer, allocatable :: displs(:), rcounts(:)

!======================================

  call mpi_init(ierr)
  comm = mpi_comm_world
  call mpi_comm_rank(comm, procno, ierr)
  call mpi_comm_size(comm, numprocs, ierr)

  root = 0

  ncols = procno

!======================================

  call mpi_type_vector(1, nrows, nrows, mpi_real4, vec_type, ierr)
  call mpi_type_commit(vec_type, ierr)

!======================================

  allocate(vec(nrows,ncols))
  do i=1,nrows
    do j=1,ncols
      vec(i,j) = procno + (0.1 * j)
    enddo
  enddo

  call mpi_barrier(comm, ierr)

!======================================

  if (procno .eq. root) then
    allocate(Gvec(nrows, ncols*numprocs))
    allocate(displs(numprocs))
    allocate(rcounts(numprocs))

    rcounts = (ncols)
    !print *, rcounts

    ! calculate displacements
    offset = 0
    do iProcs = 1,numprocs
      displs(iProcs) = offset
      offset = offset + rcounts(iProcs)
    enddo
    !print *, displs
  endif 

!======================================

  call mpi_gatherv(vec, ncols, vec_type, &
                   Gvec, rcounts, displs, vec_type, &
                   root, comm, ierr)

  if (procno .eq. root) then
    print *, shape(Gvec)
    do i=1,nrows
      print *, Gvec(i,:)
    enddo
  endif

!======================================
  if (allocated(vec)) deallocate(vec)
  if (allocated(Gvec)) deallocate(Gvec)
  if (allocated(displs)) deallocate(displs)
  if (allocated(rcounts)) deallocate(rcounts)

  call mpi_type_free(vec_type, ierr)
  call mpi_finalize(ierr)
end program main