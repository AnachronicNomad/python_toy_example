program main
  use mpi

  implicit none
  !! MPI specific
  integer :: ierr, procno, nprocs, comm, vec_type
  integer :: worker_status(MPI_STATUS_SIZE), root_status(MPI_STATUS_SIZE)
  
  !! loop iterators, indices
  integer :: iProcs, offset, ncols, i, j, prob_index
  
  integer, parameter :: nrows = 5, root = 0, num_bundles = 3
  integer, parameter :: NUM_ELEMENTS = 10
  integer, parameter :: DONE = 0
  integer :: target_proc, complete_flag, recv_length

  ! Operating data
  integer, allocatable :: root_data(:), worker_data(:), new_data(:)
  integer, allocatable :: worker_buffer
  integer :: recv_tag, worker_item, send_tag, root_item

  type bundle_type
    integer, allocatable :: data1(:)
    integer, allocatable :: data2(:)
  end type bundle_type

  type(bundle_type), allocatable :: root_bundles(:), worker_bundles(:)

!======================================

  call mpi_init(ierr)
  comm = mpi_comm_world
  call mpi_comm_rank(comm, procno, ierr)
  call mpi_comm_size(comm, nprocs, ierr)

!======================================

  if (procno .eq. root) then
    allocate(root_bundles(1:NUM_ELEMENTS))

    do i=1,NUM_ELEMENTS
      allocate(root_data(i))
      root_data = i

      root_bundles(i)%data1 = root_data

      deallocate(root_data)      
    enddo
  endif

  if (procno .eq. root) then
    do i=1,NUM_ELEMENTS
      target_proc = modulo(i, (nprocs-1))+1
      send_tag = i

      call mpi_send(root_bundles(i)%data1, i, MPI_INTEGER, &
                    target_proc, send_tag, comm, ierr)
    enddo

    do i=1,nprocs-1
      call mpi_send(0, 1, MPI_INTEGER, &
                    i, DONE, comm, ierr)
    enddo

    !!!!!
    do i=1, NUM_ELEMENTS
      call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, root_status, ierr)
      recv_tag = root_status(MPI_TAG)
      call mpi_get_count(root_status, MPI_INTEGER, recv_length, ierr)

      allocate(new_data(recv_length))
      call mpi_recv(new_data, recv_length, MPI_INTEGER, &
                    MPI_ANY_SOURCE, MPI_ANY_TAG, comm, root_status, ierr)

      root_bundles(recv_tag)%data2 = new_data

      if (allocated(new_data)) deallocate(new_data)
    enddo

    do i=1, NUM_ELEMENTS
      print *, i
      print *, root_bundles(i)%data1
      print *, root_bundles(i)%data2

      print *, ":::"
    enddo

    deallocate(root_bundles)
  else
    allocate(worker_bundles(1:NUM_ELEMENTS))

    do while(0 .eq. 0)
      call mpi_probe(root, MPI_ANY_TAG, comm, worker_status, ierr)
      recv_tag = worker_status(MPI_TAG)

      if (recv_tag .eq. DONE) then
        exit
      endif

      call mpi_get_count(worker_status, MPI_INTEGER, recv_length, ierr)
      !print *, "proc", procno, "recv_length", recv_length
      allocate(worker_data(1:recv_length))
      
      call mpi_recv(worker_data, recv_length, MPI_INTEGER, &
                    root, MPI_ANY_TAG, comm, worker_status, ierr) 

      do j=1,recv_length
        worker_data(j) = worker_data(j) + 1
      enddo

      worker_bundles(recv_tag)%data2 = worker_data 

      if (allocated(worker_data)) deallocate(worker_data)
    enddo

    !print *, "proc", procno, ":", worker_data

    !!!!!
    do i=1, NUM_ELEMENTS
      if (modulo(i,nprocs-1)+1-procno .eq. 0) then
        send_tag = i
        call mpi_send(worker_bundles(i)%data2, i, MPI_INTEGER, &
                      root, send_tag, comm, ierr)
      endif
    enddo

    deallocate(worker_bundles)
  endif

!======================================

  call mpi_finalize(ierr)
end program main