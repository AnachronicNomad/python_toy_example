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
  type(resamp_bin_type)

  type resamp_bin_type
   integer :: npart                          !< particles in bin(old)
   real (8), allocatable:: constraint(:)     !< bin moments
   real (8), allocatable:: Cmat(:,:)         !< constraint matrix for new particles
   integer :: nconstraint                    !< Number of constraints for optimization
  end type resamp_bin_type

contains
  subroutine simulation_init(inode1, inode2, vp_max, mu_max)
    !
    call MPI_INIT(mpi_err)
    call MPI_Comm_dup(MPI_COMM_WORLD, sml_comm, mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, procno, mpi_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, mpi_err)
    !
  end subroutine

  subroutine share_mats()
    !
    type (resamp_bin_type), allocatable :: bins(:,:,:)
    integer :: node,bmu,bvp



    !
  end subroutine

  subroutine sim()
    !

    !
  end subroutine

  subroutine teardown()
    !


    call MPI_FINALIZE(mpi_err)
    !
  end subroutine


end module example_mod

program example_prog
   !use mpi
   use example_mod

   implicit none

   call simulation_init()

   call sim()

   call teardown()

end program example_prog

