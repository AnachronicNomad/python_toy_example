# file: simulation.py 
import numpy as np
from mpi4py import MPI
from data_structures.sim_structures import *


def sim(num_planes, max_particles_per_bin, 
        r=(0.5, 1.5), z=(0.0, 1.0), v=(0.0, 1.0),
        vgrid_size=2,
        debug=False):
    """
        Every process will be running this function. As a result, this func
        needs to check for its rank and determine what it is actually doing.

        The assumption currently held by the program is that each bin in 
        the phase space φ,r,z,v is represented by an MPI process. Care 
        needs to be taken in how the arguments are chosen when running the
        simulation given how the MPI processes are created, then. 

        There was an alternative here but I forgot what it was right now

        The default args for r,z,v represent domain of values for those
        args. 
    """
    # boilerplate
    comm_world = MPI.COMM_WORLD # can be an arbitrary communicator
    num_procs = comm_world.Get_size() # Remember, this is num. of bins
    rank = comm_world.Get_rank() # what bin 
    
    # Here, we're going to use divisions of ranks to determine what 
    # spatial values each bin is tied to, and what velocity it will
    # represent, and which phi-slice (plane) it's in
    
    plane = np.int(np.floor(rank / num_planes))
    num_bins_per_plane = np.int(np.floor(num_procs / num_planes))

    # how we might rank the bins per plane, this process' rank inside the 
    # plane we're currently assigning it to
    per_plane_rank = rank % num_bins_per_plane

    """
        Assume the bins are ordered correctly and have been assigned their
        r,z,v values here.

        Going to start with just having spatial locations realized, not
        going to split up velocity space at all.
    """
    r_point = None
    z_point = None
    bin_ = None

    # if the number of bins per plane is square
    if (np.sqrt(num_bins_per_plane) % 1 == 0.0):
        # allocate bins in the right chunks
        
        
    else:
        # Realistically the right way to do this is to make a list of divisors
        # and choose the "middle" two (of an orderd list of divisors), then have
        # the r,z parts map correctly to that. 
        raise NotImplementedError("haven't decided how to do this yet")



    if debug:
        #print(f"""Process {rank} on Plane {plane} (plane-ordered num) {per_plane_rank} with {num_bins_per_plane} bins per plane""")

    
    
