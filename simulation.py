# file simulation.py 
import numpy as np
from mpi4py import MPI
from data_structures.sim_structures import *

def sim(num_planes, max_particles_per_bin, 
        r=(0.5, 1.5), z=(0.0, 1.0), v=(0.0, 1.0),
        v_grid_size=2,
        debug=False):
    """
        Every process will be running this function. As a result, this func
        needs to check for its rank and determine what it is actually doing.

        The assumption currently held by the program is that each bin in 
        the phase space φ,r,z,v is represented by an MPI process. Care 
        needs to be taken in how the arguments are chosen when running the
        simulation given how the MPI processes are created, then. 

        The default args for r,z,v represent domain of values for those
        args. 
    """
    # boilerplate
    comm_world = MPI.COMM_WORLD
    num_procs = comm_world.Get_size() # Remember, this is num. of bins
    rank = comm_world.Get_rank() # what bin 
    
    # Here, we're going to use divisions of ranks to determine what 
    # spatial values each bin is tied to, and what velocity it will
    # represent, and which phi-slice (plane) it's in
    
    plane = np.int(np.floor(rank / num_planes))
    num_bins_per_plane = np.int(np.floor(num_procs / num_planes))

    # Now, based on its rank, this bin needs an r,z,(v0,v1) assignment. 
    per_plane_rank = rank % num_bins_per_plane
    
