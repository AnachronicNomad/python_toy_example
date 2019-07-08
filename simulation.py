# file: simulation.py 
import numpy as np
from mpi4py import MPI
from data_structures.sim_structures import *

def determine_spatial_coords_from_global_rank(comm_world, num_procs, rank,
                                              num_planes,
                                              r=(0.5, 1.5), z=(0.0), v=(0.0,1.0)):
    """
        If there are more processes than planes, split up v-space
        Assuming a single r,z coordinate for each phi_plane
    """
    r_point = (r[1] - r[0] / 2)
    z_point = (z[1] - z[0] / 2)
    v_range = (v[0], v[1])
    phi_plane = 0

    if num_procs < num_planes:
        raise ValueError("There must be at least NUM_PROCS to represent the phi-planes")
    elif (num_procs / num_planes) % 1.0 != 0.0:
        # if the number of processes aren't a multiple of the number of planes, bad
        raise ValueError("NUM_PROCS must be multiple of NUM_PLANES")
    elif num_planes == num_procs:
        phi_plane = rank
    else:
        # we have multiples of phi_planes and want to split up v_space
        multiple = num_procs / num_planes

        """ phi_plane here is given by block ordering. """
        phi_plane = np.int(np.floor(rank / multiple))
        plane_order = rank % multiple

        v_step = (v[1] - v[0]) / multiple
        v_range = ((v_step*plane_order) + v[0], v[0] + ((v_step*(plane_order+1))))

    return (r_point, z_point, v_range, phi_plane)


def determine_group_comm_from_global_rank(comm_world, num_procs, rank, bin):
    """
        For now, only velocity is split up
    """
    # Find scheme to assign a unique "color" to same (r,z) group
    new_comm = comm_world.Split(color=1, key=rank)

    return new_comm


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

    """Now generate the particles per bin, update constraint vec and matrix"""
    pos = determine_spatial_coords_from_global_rank(comm_world, num_procs, rank, 
                                                    num_planes, r, z, v)
    bin_ = Bin(pos[0], pos[1], pos[2], pos[3])
    for i in range(0, np.random.randint(1, max_particles_per_bin+1)):
        w = np.random.weibull(a=np.random.randint(0,7))
        v = np.random.uniform(low=bin_.v[0], high=bin_.v[1])

        bin_.add_particle(Particle(weight=w, velocity=v))

    # make sure the constraint vector is updated
    bin_.update_constraints()
    bin_.build_constraint_mat()

    print("Process", rank, "has coords plane=", bin_.plane_id, "r=", bin_.r, "z=", bin_.z,
            "v=", bin_.v)
    

    """Create group and communicator for axisymmetry"""
    axi_comm = determine_group_comm_from_global_rank(comm_world, num_procs, rank, bin_)
    
    
    return None
    
