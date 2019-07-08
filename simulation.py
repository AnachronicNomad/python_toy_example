# file: simulation.py 
import numpy as np
from mpi4py import MPI
from data_structures.sim_structures import *

def determine_spatial_coords_from_global_rank(comm_world, num_procs, rank,
                                              num_planes,
                                              r=(0.5, 1.5), z=(0.0,1.0), v=(0.0,1.0)):
    """
        If there are more processes than planes, split up v-space
        Assuming a single r,z coordinate for each phi_plane

        Note: Strongly consider having each (φ,r,z) spatial bin use MPI_spawn 
        to spawn a collection of velocity-space bins and manage these that way. 

        This would create a group/communicator that can aggregate the values from v-space
        bins to a process tied to a spatial location in (φ,r,z); and then these "vertex" 
        processes can solve the constraint problem axisymmetrically in their relevant 
        groups.
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

    return (phi_plane, r_point, z_point, v_range)


def determine_group_comm_from_global_rank(comm_world, num_procs, rank, bin):
    """
        For now, only velocity is split up. Since this is currently a reduced
        example with controlled behavior, we're just going to put them all 
        on the same communicator. 

        WARNING: In order to reduce heavy/repeated I/O on the MPI instance, 
        it's preferable to have each individual process build an MPI_Group of 
        process IDs that it correlates with, and then create the communicator
        from that group using 
    """
    # Find scheme to assign a unique "color" to identical (r,z) group
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
    global_rank = comm_world.Get_rank() # what bin 

    # Assign this "bins" spatial coordinates
    pos = determine_spatial_coords_from_global_rank(comm_world, num_procs, 
                                                    global_rank, num_planes, 
                                                    r, z, v)


    """Now generate the particles per bin, update constraint vec and matrix"""
    bin_ = Bin(pos[0], pos[1], pos[2], pos[3])
    for i in range(0, np.random.randint(1, max_particles_per_bin+1)):
        w = np.float64(np.random.weibull(a=np.random.randint(0,7)))
        v = np.float64(np.random.uniform(low=bin_.v[0], high=bin_.v[1]))

        bin_.add_particle(Particle(weight=w, velocity=v))

    # make sure the constraint vector is updated
    bin_.update_constraints()
    bin_.build_constraint_mat()

    #print("Process", global_rank, "has coords plane=", bin_.plane_id, "r=", bin_.r, "z=", bin_.z,
            #"v=", bin_.v)
    

    """Create group/communicator for axisymmetry"""
    axi_comm = determine_group_comm_from_global_rank(comm_world, num_procs, 
                                                     global_rank, bin_)

    # good practice to use local rank when dealing with group/comm specific processes
    local_rank = axi_comm.Get_rank()

    ## Gather the rank and number of particles in the bins for the root process.
    #sendbuf = np.array([ local_rank, len(bin_.particles) ], dtype='int64')
    #recvbuf_num_particles = np.empty([axi_comm.Get_size(), 2], dtype='int64')

    ## WARNING: This will force a synchronization across all processes. 
    ## There is an implicit MPI barrier in an Allgather.

    #axi_comm.Allgather(sendbuf, recvbuf_num_particles)
    #print("Process ", global_rank, "reports\n", recvbuf_num_particles)

    #num_particles = 0
    #for row in recvbuf_num_particles:
    #   num_particles += row[1]

    #print("There's a total number of ", num_particles, "particles")

    ## Gather all of the constraint matrices to the leading process for the 
    ## communicator, and have it rebuilt the full constraint matrix
    ## to send back out to the rest of the communicator. 

    # gather the constraint matrices on leading process for communicator
    data = bin_.constraint_mat
    data = axi_comm.gather(data, root=0)
    tmp = np.empty([3, 1])
    if local_rank == 0:
        for array in data:
            tmp = np.hstack([tmp, array])
        #print(tmp)

    # broadcast back to everybody in communicator
    if local_rank == 0:
        full_matrix = np.array(tmp)
    else:
        full_matrix = None

    full_matrix = axi_comm.bcast(full_matrix, root=0)

    #print("Process", global_rank, "reports recv", full_matrix.shape)


    return None
    
