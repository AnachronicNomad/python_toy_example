!!  Fortran example  
module example_mod
  use mpi

  implicit none
  ! MPI 
  integer :: procno, numprocs, sml_comm, mpi_err
  integer :: win
  integer, dimension(:), allocatable :: shared_data

  ! Sim data
  integer :: loop_cntr
  integer, dimension(:), allocatable :: local_data

contains
  subroutine simulation_init()
    !
    call MPI_INIT(mpi_err)
    call MPI_Comm_dup(MPI_COMM_WORLD, sml_comm, mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, procno, mpi_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, mpi_err)

    allocate(local_data(1:10))

    call MPI_Win_create_dynamic(MPI_INFO_NULL, sml_comm, win, mpi_err)
    !
  end subroutine

  subroutine sim()
    !
    do loop_cntr=1,10
      local_data(loop_cntr) = procno * loop_cntr
    enddo

    print *, 'proc', procno, 'has', local_data

    call MPI_Win_attach(win, local_data, 10, mpi_err)
    !
  end subroutine

  subroutine teardown()
    !
    call MPI_Win_detach(win, local_data)
    deallocate(local_data)

    call MPI_Win_free(win, mpi_err)

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

