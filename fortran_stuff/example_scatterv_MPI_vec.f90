program main
  use mpi

  implicit none
  integer :: ierr, procno, numprocs, comm, vec_type !! MPI values
  integer :: iproc, offset, ncols, i, j !! loop iterators, indices
  integer, parameter :: root = 0

  integer, allocatable :: senddata(:), recvdata(:), displs(:), rcounts(:)

!======================================

  call mpi_init(ierr)
  comm = mpi_comm_world
  call mpi_comm_rank(comm, procno, ierr)
  call mpi_comm_size(comm, numprocs, ierr)

  ncols = procno+1

  if (procno .eq. root) then
    allocate(rcounts(numprocs))
    allocate(displs(numprocs))
  endif

  call mpi_gather(ncols, 1, mpi_integer, & 
                  rcounts, 1, mpi_integer, &
                  root, comm, ierr)

  if (procno .eq. root) then
    allocate(senddata(sum(rcounts)))
    senddata = 1

    !print *, size(senddata), ",", shape(senddata)
    offset = 0
    do iproc = 1,numprocs
      displs(iproc) = offset
      offset = offset + rcounts(iproc)
    enddo
  endif

!======================================

  allocate(recvdata(ncols))

  call mpi_scatterv(senddata, rcounts, displs, mpi_integer, &
                    recvdata, ncols, mpi_integer, &
                    root, comm, ierr)

  print *, procno, ": ", recvdata

!======================================
!! Do teardown and free


  call mpi_finalize(ierr)
end program main