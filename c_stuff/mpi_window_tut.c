#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char* argv[]) {
	MPI_Init(NULL, NULL);
	int procno, numproc;
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &procno);

	int *a;
	MPI_Win win;

	/* create private memory */ 
	MPI_Alloc_mem(10*sizeof(int), MPI_INFO_NULL, &a);     
	/* use private memory like you normally would */     
	for(int i=0; i<10; i++) {
		a[i] = procno;

	/* collectively declare memory as remotely accessible */ 
	MPI_Win_create(a, 10*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);  

	/* Array ‘a’is now accessible by all processes in MPI_COMM_WORLD */ 
	
	MPI_Win_free(&win);
       	MPI_Free_mem(a);

	MPI_Finalize();
	return 0;
}
