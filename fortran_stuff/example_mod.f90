!!  Fortran example  
module example_mod
  use mpi

  implicit none
  ! MPI 
  integer :: procno, numprocs, mpi_err

contains
  subroutine simulation_init()
    call MPI_INIT(mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, procno, mpi_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, mpi_err)

    print *, 'proc',  procno

  end subroutine

  subroutine teardown()
    call MPI_FINALIZE(mpi_err)
  end subroutine


end module example_mod

program example_prog
   !use mpi
   use example_mod

   implicit none

   call simulation_init()

   call teardown()

end program example_prog

