static char help[] = "Hello world program.\n\n";
#include <petsc.h>

int main( int arc, char **argv)
{
	PetscMPIInt rank;
	PetscErrorCode ierr;

	PetscInitialize(&argc, &argv, (char*) 0, help);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_SELF, "  Hello from %d\n", rank);CHKERRQ(ierr);

	PetscFinalize();

	return 0;
}
