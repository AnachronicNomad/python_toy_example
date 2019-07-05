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
    
    plane = np.int(np.floor(rank % num_planes))
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
        root = np.sqrt(num_bins_per_plane)
        dr = (r[1] - r[0]) / root
        dz = (z[1] - z[0]) / root
        r_point = dr * (per_plane_rank % root) + r[0]
        z_point = dz * np.floor(per_plane_rank / root) + z[0]
    else:
        # Realistically the right way to do this is to make a list of divisors
        # and choose the "middle" two (of an orderd list of divisors), then have
        # the r,z parts map correctly to that. 
        raise NotImplementedError("haven't decided how to do this yet, num_planes_per_bin non-square")

    bin_ = Bin(r_point, z_point, v, plane)

    if debug:
        print("Process", rank, "on plane", plane, "plane_rank", plane_rank, "with r=", r_point, "z=", z_point)

    """Now generate the particles per each bin"""
    for i in range(0, np.random.randint(1, max_particles_per_bin+1)):
        w = np.random.weibull(a=np.random.randint(0,7))
        v = np.random.uniform(low=v[0], high=v[1])

        bin_.add_particle(Particle(weight=w, velocity=v))

    # make sure the constraint vector is updated
    bin_.update_constraints()

    
    
