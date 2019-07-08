!!  Fortran example  
program hello
   use mpi

   implicit none
   integer world_rank, world_size, ierror, tag, status(MPI_STATUS_SIZE)
   
   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD, world_rank, ierror)

   if (modulo(world_rank, 2) .eq. 0) then
        print *, 'node', world_rank, ': I am even.'
   else
        print *, 'node', world_rank, ': I am odd.'
   end if

   call MPI_FINALIZE(ierror)
end program hello

