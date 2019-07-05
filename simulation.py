# file simulation.py 
from mpi4py import MPI
from data_structures.sim_structures import *

def sim(num_planes, num_bins_per_plane, max_particles_per_bin, 
        r=(0.5, 1.5), z=(0,1), v=(0,1000),
        debug=False):
    """
        Every process will be running this function. As a result, this func
        needs to check for its rank and determine what it is actually doing.
    """
    comm_world = MPI.COMM_WORLD

    num_procs = comm_world.Get_size()
    rank = comm_world.Get_rank()

    print(f"""Hello, I'm process {rank} of {num_procs}""")
