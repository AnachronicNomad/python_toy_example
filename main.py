# file: main.py
from mpi4py import MPI
from data_structures.sim_structures import *
from simulation import sim

def main():
    sim(4, 3, debug=True)


if __name__ == '__main__':
    main()
