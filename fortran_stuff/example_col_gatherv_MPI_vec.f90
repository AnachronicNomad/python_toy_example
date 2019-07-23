program main
  !
  ! TODO:
  !   1. DONE, blocks and offset size working correctly.
  !   2. ~~See if we can exploit Column-major ordering and contiguous memory 
  !     in Fortran to easily count and send the columns to be stacked on the
  !     receiving end. ~~ <- DONE, it works!!
  !   3. ~~Have arbitrary process ordering (rcount)~~ <- DONE
  !
  !
  use mpi

  implicit none
  integer :: ierr, procno, numprocs, root, comm, vec_type !! MPI values
  integer :: iProcs, offset, ncols, i, j !! loop iterators, indices
  integer, parameter :: nrows = 5

  real(kind=4), allocatable :: vec(:,:), Gmat(:,:)
  integer, allocatable :: displs(:), rcounts(:)

!======================================

  call mpi_init(ierr)
  comm = mpi_comm_world
  call mpi_comm_rank(comm, procno, ierr)
  call mpi_comm_size(comm, numprocs, ierr)

  root = 0

  ncols = procno+1

!======================================
!! Define the vector type that we're going to use to transfer columns 
!! of 5x{1->procno} sub-matrix on each process

  call mpi_type_vector(1, nrows, nrows, mpi_real4, vec_type, ierr)
  call mpi_type_commit(vec_type, ierr)

!======================================
!! Have each process create its number of columns, then 
!! gather the column count per process on the main communicator

  !
  ! For easy printout, leading digit is proc_id in comm, following digit is 
  ! the column that the value is supposed to have.
  !
  allocate(vec(nrows,ncols))
  do i=1,nrows
    do j=1,ncols
      vec(i,j) = real(10 * procno + j,kind=4)
    enddo
  enddo

  !print *, 'proc', procno, 'completed filling array with nrows',  &
  !         nrows, 'ncols', ncols
  !print *, shape(vec)
  !do i=1,nrows
  !  print *, vec(i,:)
  !enddo

  if (procno .eq. root) then
    allocate(rcounts(numprocs))
  endif

  call mpi_gather(ncols, 1, mpi_integer, &    !! Sending (sendbuf, num_elements, send_type)
                  rcounts, 1, mpi_integer, &  !! Receiving (recvbuf, num_elements, recv_type)
                  root, comm, ierr)           !! (receiving root, communicator, error handler)

  !if (procno .eq. root) then
  !  print *, rcounts 
  !endif

!======================================
!! Root (receving process) will use the total column counts (rcounts)
!! to allocate the Global matrix. 
!! Then, the displacement/offset of memory region inside the 
!! global matrix is just multiples of data type that we defined

  if (procno .eq. root) then
    allocate(Gmat(nrows, sum(rcounts)))
    allocate(displs(numprocs))

    !rcounts = (ncols)
    !print *, rcounts

    ! calculate displacements of each vector
    offset = 0
    do iProcs = 1,numprocs
      displs(iProcs) = offset
      offset = offset + rcounts(iProcs)
    enddo
    !print *, displs
  endif 

!======================================
!! Use "variable-length" gather to place all the matrices in order
!! onto the root process, which will do the indexing. 

  call mpi_gatherv(vec, ncols, vec_type, &
                   Gmat, rcounts, displs, vec_type, &
                   root, comm, ierr)

  !if (procno .eq. root) then
  !  print *, shape(Gmat)
  !  do i=1,nrows
  !    print *, Gmat(i,:)
  !  enddo
  !endif

!======================================
!! Do teardown and free

  if (allocated(vec)) deallocate(vec)
  if (allocated(Gmat)) deallocate(Gmat)
  if (allocated(displs)) deallocate(displs)
  if (allocated(rcounts)) deallocate(rcounts)

  call mpi_type_free(vec_type, ierr)
  call mpi_finalize(ierr)
end program main