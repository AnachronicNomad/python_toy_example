This is a toy example to make sure I know what I'm talking about.

# Installation & Runtime

This assumes you have a working MPI implementation. The author created this
with MPICH. 

There is also a requirements.txt file which contains the package list for the 
python3 virtualenv used to run this code. 

## Running/Usage

In order to run the program on the command line, use
```
mpiexec -n {_numprocs_} python -m mpi4py {_pyfile_} [arg]
```

The author decided to use main.py as the entry point for the simulation.
