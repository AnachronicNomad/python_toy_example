!
! The purpose of this module is to show an `nconstraint` row, with `npart`
! columns matrix, in addition to a summed constraint vector; shared 
! "axisymmetrically" 
!
module example_mod
  use mpi

  implicit none
  ! MPI 
  integer :: procno, numprocs, sml_comm, mpi_err
  integer :: axi_color
  integer :: axi_comm

  ! Sim data
  integer, private :: resamp_musize             !< Number of sampling bins in v_perp/mu direction
  integer, private :: resamp_vpsize             !< Number of sampling bins in v_para direction
  integer, private :: resamp_inode1             !< Index of the first node of local patch of the configuration space mesh (=f0_inode1)
  integer, private :: resamp_inode2             !< Index of the last node of local patch of the configuration space mesh (=f0_inode2)

  type resamp_bin_type
   integer :: npart                          !< particles in bin(old)
   integer :: bin_id(3)                      !< corresponding node number for Voronoi cell,velocity
   !type(ptl_type), allocatable :: new_ptl(:) !< new particles generated in the bin
   real (8), allocatable:: constraint(:)     !< bin moments
   real (8), allocatable:: Cmat(:,:)         !< constraint matrix for new particles
   integer :: nconstraint                    !< Number of constraints for optimization

  end type resamp_bin_type

  type(resamp_bin_type), allocatable, private :: bins(:,:,:)


contains
  subroutine simulation_init(inode1, inode2, vp_max, mu_max)
    !
    implicit none
    integer, intent(IN) :: inode1, inode2, vp_max, mu_max
    integer :: node,bmu,bvp,i                             !< loop counters

    resamp_inode1 = inode1
    resamp_inode2 = inode2
    resamp_vpsize = vp_max
    resamp_musize = mu_max

    allocate(bins(resamp_inode1:resamp_inode2,resamp_musize,resamp_vpsize))

    !
    ! Put some values into the bins
    !
    do bvp=1,vp_max
      do bmu=1,mu_max
        do node=inode1,inode2
          bins(node,bmu,bvp)%npart = modulo(procno, numprocs) + 1
          bins(node,bmu,bvp)%nconstraint = 5
          allocate(bins(node,bmu,bvp)%constraint(1:bins(node,bmu,bvp)%nconstraint))
          allocate(bins(node,bmu,bvp)%Cmat( 1:bins(node,bmu,bvp)%nconstraint, 1:bins(node,bmu,bvp)%npart) )

          do i=1,bins(node,bmu,bvp)%nconstraint
            bins(node,bmu,bvp)%constraint(i) = real(i+bvp+bmu, 8)
          enddo

          do i=1,bins(node,bmu,bvp)%npart
            bins(node,bmu,bvp)%Cmat(1,i)=bins(node,bmu,bvp)%npart
            bins(node,bmu,bvp)%Cmat(2,i)=bins(node,bmu,bvp)%npart
            bins(node,bmu,bvp)%Cmat(3,i)=bins(node,bmu,bvp)%npart
            bins(node,bmu,bvp)%Cmat(4,i)=bins(node,bmu,bvp)%npart
            bins(node,bmu,bvp)%Cmat(5,i)=bins(node,bmu,bvp)%npart
          enddo
        enddo
      enddo
    enddo
    !
  end subroutine

  subroutine share_mats(bins, nconstraint, num_part, Axi_Mat)
    !
    implicit none
    type (resamp_bin_type), intent(INOUT) :: bins(:,:,:)
    real (8), allocatable, intent(INOUT) :: Axi_Mat
    integer, intent(IN) :: nconstraint

    integer :: node,bmu,bvp,iProcs,offsets                 !< Loop counter
    integer :: axi_pe, axi_pe_size, root, vec_type  !< MPI_Object reference handles
    integer :: nparts
    real (8), allocatable :: constraint_sum(:)
    integer, allocatable :: nparts_bins(:), displs(:)

    ! 
    !  Create an axisymmetric communicator. 
    ! Chances are, this can actually be handled by sml_intpl_comm, available from sml_module
    !

    axi_color = 1 !modulo((resamp_inode2 - resamp_inode1), procno+1)
    ! ^^^ this can also be done by modulo(grid%nnode, inode1) as long as node numbering is 1-indexed
    call MPI_Comm_split(sml_comm, axi_color, procno, axi_comm, mpi_err)
    call MPI_COMM_RANK(axi_comm, axi_pe, mpi_err)
    call MPI_COMM_SIZE(axi_comm, axi_pe_size, mpi_err)

    root = 0

    !======================================
    !! Define the vector type that we're going to use to transfer columns 
    !! of {nconstraint}x{1->procno} sub-matrix on each process

    call mpi_type_vector(1, nrows, nrows, mpi_real8, vec_type, mpi_err)
    call mpi_type_commit(vec_type, mpi_err)

    !======================================

    if (axi_pe .eq. root) then
      allocate(constraint_sum(1:nconstraint))
      allocate(nparts_bins(1:axi_pe_size))
      allocate(displs(axi_pe_size))
    endif

    do bvp=1,vp_max
      do bmu=1,mu_max
        do node=inode1,inode2
          call MPI_Barrier(axi_comm, mpi_err)

          !print *, 'proc', procno, 'node', node, 'bmu', bmu, 'bvp', bvp

          !
          ! Reduce the constraint vector onto root process of the axisymmetric communicator
          !
          ! The sum of constraint vectors in the axisymmetric communicator is now there. 
          !

          call MPI_Reduce(bins(node,bmu,bvp)%constraint, nconstraint, MPI_REAL8, &
                          MPI_SUM, root, axi_comm, mpi_err)

          nparts = bins(node,bmu,bvp)%npart ! + size(bins(node,bmu,bvp)%new_ptl)

          ! count how many particles there are per bin 
          call MPI_Gather(nparts, 1, MPI_INT, & ! sendbuf, sendcount, sendtype
                          nparts_bins, 1, MPI_INT, &              ! recvbuf, recvcount, recvtype
                          root, axi_comm, mpi_err)                   ! root, comm, error_flag


          !======================================
          !! Root (receving process) will use the total column counts (rcounts)
          !! to allocate the Global matrix. 
          !! Then, the displacement/offset of memory region inside the 
          !! global matrix is just multiples of data type that we defined
          if (axi_pe .eq. root) then
            allocate(Axi_Mat(nconstraint, sum(nparts_bins)))
            allocate(displs(axi_pe_size))

            ! calculate displacements of each vector
            offset = 0
            do iProcs = 1,axi_pe_size
              displs(iProcs) = offset
              offset = offset + nparts_bins(iProcs)
            enddo
            !print *, displs
          endif

          !======================================
          !! Use "variable-length" gather to place all the matrices in order
          !! onto the root process, which will do the indexing. 

          call mpi_gatherv(bins(node,bmu,bvp)%Cmat, nparts, vec_type, &
                           Axi_Mat, nparts_bins, displs, vec_type, &
                           root, axi_comm, mpi_err)

          !if (axi_pe .eq. root) then
          !  print *, 'shape: ', shape(Gmat)
          !  do i=1,nrows
          !    print *, Gmat(i,:)
          !  enddo
          !endif
          
          !
          ! Do something
          !
          if (axi_pe .eq. root) then
            ! solve optimization here
          endif

          !
          ! Broadcast back to axisymmetric communicator
          !
          
          call MPI_Barrier(axi_comm, mpi_err)
        enddo
      enddo
    enddo

    !
  end subroutine

  subroutine sim()
    !
      implicit none

      call share_mats()
    !
  end subroutine

  subroutine teardown()
    !
    deallocate(bins)

    call MPI_FINALIZE(mpi_err)
    !
  end subroutine


end module example_mod

program example_prog
  !use mpi
  use example_mod

  implicit none

  call MPI_INIT(mpi_err)
  call MPI_Comm_dup(MPI_COMM_WORLD, sml_comm, mpi_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD, procno, mpi_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, mpi_err)

  call simulation_init(3*procno+1, 3*procno+3, 2, 2)

  call sim()

  call teardown()

end program example_prog

