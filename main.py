# file: main.py
from mpi4py import MPI
from data_structures.sim_structures import *

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    print(f"""Hello, I'm process {rank} of {size}""")


if __name__ == '__main__':
    main()
