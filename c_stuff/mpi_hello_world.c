#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Takes a number of processes, divides them into the supplied
 * number of 'groups'.  
 * In this example, groups aren't actually created. It just uses
 * Comm_split, using the process rank to determine 'color' for 
 * the shared intra-communicator to use. 
 */
int main(int argc, char* argv[]) {
	if(argc != 2)
	{
		printf("usage: mpirun -n {num_procs} ./hello {num_groups}\n");
		exit(1);
	}

	// Initialize and get process rank, size of MPI_COMM_WORLD (# created processes)
	MPI_Init(NULL, NULL);

	int world_size, world_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	// capture command line arg, check if valid
	int mod_factor = atoi(argv[1]);
	if (mod_factor > world_size)
	{
		printf("Number of groups larger than number of processes, abort.\n");
		exit(1);
	}

	// Make a copy of COMM_WORLD to make edits, different 'context' than COMM_WORLD
	// which is mutable
	MPI_Comm sml_world;
	MPI_Comm_dup(MPI_COMM_WORLD, &sml_world);

	// Create the new communicator
	int color = world_rank % mod_factor;

	MPI_Comm new_comm;
	MPI_Comm_split(sml_world, color, world_rank, &new_comm);

	// Determine local rank in the communicator
	int local_rank;
	MPI_Comm_rank(new_comm, &local_rank);

	printf("Process %d on comm %d with local rank %d\n", world_rank, color, local_rank);

	// If we're the lowest ranked on the communicator, broadcast our global rank to the 
	// other processes in our communicators, so they can send their values back to us
	
	int sum;
	MPI_Reduce(&world_rank, &sum, 1, MPI_INT, MPI_SUM, 0, new_comm);

	if (local_rank == 0)
	{
		printf("The sum of process ranks on comm %d is %d\n", color, sum);
	}

	MPI_Comm_free(&sml_world);
	MPI_Comm_free(&new_comm);
	MPI_Finalize();
}
