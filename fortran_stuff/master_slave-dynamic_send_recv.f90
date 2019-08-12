program main
  use mpi

  implicit none
  !! MPI specific
  integer :: ierr, procno, nprocs, comm, vec_type
  integer :: worker_status(MPI_STATUS_SIZE)
  
  !! loop iterators, indices
  integer :: iProcs, offset, ncols, i, j, prob_index
  
  integer, parameter :: nrows = 5, root = 0, num_bundles = 3
  integer, parameter :: NUM_ELEMENTS = 100
  integer, parameter :: DONE = 0
  integer :: target_proc, complete_flag

  ! Operating data
  integer, allocatable :: root_data(:), alt_data(:), local_data(:)
  integer :: recv_data, other_recv_data, recv_tag, temp

!======================================

  call mpi_init(ierr)
  comm = mpi_comm_world
  call mpi_comm_rank(comm, procno, ierr)
  call mpi_comm_size(comm, nprocs, ierr)

!======================================

  if (procno .eq. root) then
    allocate(root_data(1:NUM_ELEMENTS))
    allocate(alt_data(1:NUM_ELEMENTS))
    ! TODO: USE alt_data TO TEST SENDING WITH THE SAME TAG
    do i=1,NUM_ELEMENTS
      root_data(i) = i
      alt_data(i) = i
    enddo
  endif

  if (procno .eq. root) then
    do i=1,NUM_ELEMENTS
      target_proc = modulo(i, (nprocs-1))+1
      call mpi_send(root_data(i), 1, MPI_INTEGER, &
                    target_proc, i, comm, ierr)
      call mpi_send(alt_data(i), 1, MPI_INTEGER, &
                    target_proc, i, comm, ierr)
    enddo

    do i=1,nprocs-1
      call mpi_send(0, 1, MPI_INTEGER, &
                    i, DONE, comm, ierr)
    enddo
  else
    allocate(local_data(NUM_ELEMENTS))
    local_data = -1

    do while(0 .eq. 0)
      call mpi_probe(root, MPI_ANY_TAG, comm, worker_status, ierr)
      recv_tag = worker_status(MPI_TAG)
      ! Will have to use MPI_get_count here when sending variable sized arrays

      if (recv_tag .eq. DONE) then
        exit
      endif

      call mpi_recv(recv_data, 1, MPI_INTEGER, &
                    root, MPI_ANY_TAG, comm, worker_status, ierr) 
      call mpi_recv(other_recv_data, 1, MPI_INTEGER, &
                    root, MPI_ANY_TAG, comm, worker_status, ierr) 

      local_data(recv_tag) = 

      print *, "proc", procno, "result", temp, "tag", recv_tag

    enddo

  endif

!======================================

  call mpi_finalize(ierr)
end program main