program main
  !
  ! Answer taken from https://stackoverflow.com/a/38302499, 
  ! commentary and some refactoring for understanding done by anachronicnomad. 
  !
  use mpi

  implicit none
  integer :: ierr, myRank, iProcs, nProcs, master, newtype !! MPI values
  integer :: ix, iy, ip !! loop counters
  integer :: nxSum, offset

  integer :: nx
  integer, parameter :: ny=5 ! Each process has a matrix with 5 rows

  integer, allocatable :: vec(:), vecG(:), nxAll(:), displs(:), rcounts(:), matG(:,:)

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world, myRank, ierr)
  call mpi_comm_size(mpi_comm_world, nProcs, ierr)

!
! Have each process build its local row collection in a 1D array
!

  master = 0

  nx = myRank+1
  allocate(vec(nx*ny))

  do ix = 1,nx
  ! For each row
     do iy = 1,ny
     ! For each column
        ! Update every index of local vec array with this processes' rank; 
        ! index = (row_count) + col_offset
        ip = (ix-1)*ny + iy
        vec(ip) = myRank
     enddo
  enddo

  call mpi_barrier(mpi_comm_world,ierr)

!
! Calculate/Allocate the necessary counts for the later gatherv op for root receiver
!
  
  ! Gather the row count for all processes
  allocate(nxAll(nProcs))
  call mpi_gather(nx, 1, mpi_integer, &         ! sendbuf, sendcount (per process), sendtype
                  nxAll, 1, mpi_integer, &      ! recvbuf, recvcount (per process), rectype
                  master, mpi_comm_world, ierr) ! recv_proc, comm_group, ierr

  ! Do collection on root process
  if (myRank == master) then
     ! print *, 'nxAll = ', nxAll, 'sum(nxAll) = ',sum(nxAll)

     nxSum = sum(nxAll)        ! Total number of rows 
     allocate(vecG(nxSum*ny))  ! Allocate 1D buffer on root process to hold the 
                               ! 1D array representation from each sender

     allocate(displs(nProcs))
     allocate(rcounts(nProcs))
     offset = 0
     do iProcs = 1,nProcs
        displs(iProcs) = offset
        rcounts(iProcs) = nxAll(iProcs)*ny
        offset = offset + rcounts(iProcs)
        ! print *,'iProcs',iProcs,'displs = ',displs(iProcs),'rcounts',rcounts(iProcs)
     enddo

  endif

!
! Type definition for transfer, and the gathering mechanism
!

  call mpi_type_vector(nx*ny, & ! count - # blocks
                       1, & ! block length - # elements (oldtype) in each block
                       1, & ! stride - num elements between start of each block
                       mpi_integer, & ! old_type (what this vector is constituent of)
                       newtype, & ! type reference for MPI_Object for rest of prog
                       ierr) ! err_flag for MPI
  ! Add the MPI_Object reference type for MPI Implementation
  call mpi_type_commit(newtype,ierr)

  ! Gather, with variable displacement, the values to place. 
  call mpi_gatherv(vec, & ! sendbuf
                   1, & ! send count - num_elements in send buffer
                   newtype, & ! send type - handle for MPI_Object 
                   vecG, & ! receive buffer on root process
                   rcounts, ! recvcounts - int array of num_elements recv from each proc
                   displs, ! integer array (length grp size); entry i specifies displacement rel. to recvbuf at which to place incoming data
                   mpi_integer, & ! rec_type
                   master, & !recv process_id
                   mpi_comm_world, & ! communicator to send on
                   ierr) & ! MPI Error flag handle


  if (myRank == master) then
     print *, 'Global vector, vecG = ',vecG

     ! Reshape into matrix
     print *, 'Global matrix'
     allocate(matG(nxSum,ny))
     do ix = 1,nxSum
        do iy = 1,ny
           ip = (ix-1)*ny + iy
           matG(ix,iy) = vecG(ip)
        enddo
        print *, (matG(ix,iy),iy=1,ny)
     enddo

  endif

  call mpi_finalize(ierr)

end program main