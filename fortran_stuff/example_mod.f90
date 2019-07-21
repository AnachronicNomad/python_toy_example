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
    do node=inode1,inode2
      do bmu=1,mu_max
        do bvp=1,vp_max
          bins(node,bmu,bvp)%npart = modulo(procno, numprocs) + 1
          bins(node,bmu,bvp)%nconstraint = 5
          allocate(bins(node,bmu,bvp)%constraint(1:bins(node,bmu,bvp)%nconstraint))
          allocate(bins(node,bmu,bvp)%Cmat(bins(node,bmu,bvp)%nconstraint, bins(node,bmu,bvp)%npart))

          do i=1,bins(node,bmu,bvp)%nconstraint
            bins(node,bmu,bvp)%constraint = real(1, 8)
          enddo

          bins(node,bmu,bvp)%Cmat = real(procno, 8) 
        enddo
      enddo
    enddo
    !
  end subroutine

  subroutine share_mats()
    !
    implicit none
    !type (resamp_bin_type), allocatable, intent(INOUT) :: bins(:,:,:)
    integer :: node,bmu,bvp

    axi_color = (resamp_inode2 - resamp_inode1)
    call MPI_Comm_split(sml_comm, axi_color, procno, axi_comm, mpi_err)

    do node=resamp_inode1,resamp_inode2
      do bmu=1,resamp_mu_max
        do bvp=1,resamp_vp_max

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

  call simulation_init(3*procno+1, 3*procno+3, 4, 4)

  call sim()

  call teardown()

end program example_prog
