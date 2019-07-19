# file: simulation.py 
import numpy as np
from mpi4py import MPI
from data_structures.sim_structures import *


def sim(num_nodes=8, vp_max=1.0, mu_max=1.0, npart_max=5, 
        nnodes_plane=4, vpsize=1, musize=1):

    #
    # This is set-up to emulate what's going on in the real
    #

    # MPI Init not needed
    sml_comm = MPI.COMM_WORLD
    num_procs = sml_comm.Get_size()
    procno = sml_comm.Get_rank()

    #
    # build grid on main process, send to all processes.
    #
    if procno == 0:
        grid = []
        for i in range(0,num_nodes):
            print("Node", i, "has particles: ")
            particles = []
            for npart in range(1, np.int(np.random.randint(low=1, high=npart_max))+1):
                w = np.float64(np.random.weibull(a=np.random.randint(0,7))) + 1e-5
                vp = np.float64(np.random.uniform(low=1e-5, high=vp_max))
                mu = np.float64(np.random.uniform(low=1e-5, high=mu_max))

                print("weight", w, "v_para", vp, "vperp", mu)

                particles.append(Particle(weight=w, vp=vp, mu=mu))

            grid.append([i, # node_id
                         particles]) 
    else:
        grid = None


    #print("Proc ", procno, "hitting barrier")
    sml_comm.Barrier()
    grid = sml_comm.bcast(grid, root=0)
    #print("Proc", procno, "has grid", grid)
    #sml_comm.Barrier()

    #
    # Assign nodes to processes
    # 
    nnodes_proc = np.int(num_nodes / num_procs)

    inode1 = np.int(procno * nnodes_proc)
    inode2 = np.int(inode1 + nnodes_proc - 1)

    #print("Proc", procno, "has nodes", inode1, "-", inode2)

    #
    # Create bins and do constraint matrix allotment
    #
    bins = np.empty(shape=(nnodes_proc, vpsize, musize), dtype=object)
    #print(bins.shape)

    for bin_vp in range(0, vpsize):
        for bin_mu in range(0, musize):
            for nnode in range(inode1, inode2+1):
                bin_ = Bin()

                for ptl in grid[nnode][1]:
                    bin_.add_particle(ptl)

                #print(bin_)

                bin_.update_constraints()
                bin_.build_constraint_mat()

                node_idx = np.int(nnode % nnodes_proc)
                bins[node_idx, bin_vp, bin_mu] = bin_

    #print("Proc", procno, "has the following constraint matrices:")
    #for bin_vp in range(0, vpsize):
    #    for bin_mu in range(0, musize):
    #        for nnode in range(inode1, inode2+1):
                #print("v para", bin_vp, "mu", bin_mu, "node", nnode)
    #            node_idx = np.int(nnode % nnodes_proc)
                #print(bins[node_idx, bin_vp, bin_mu].constraint_mat)


    #
    # Share the constraint matrices
    #
    # The problem description is that: 
    #
    # 1. Each bin in (v_para, v_perp/mu) has some number of particles.
    # 2. Each bin has a constraint matrix.
    #       - rows == number of constraints
    #       - cols == number of particles in bin
    # 3. Each bin also has constraints (sum for all particles in bin)
    #
    # Question - Is it ONLY the constraint matrix that is being shared, once, axisymmetrically?
    # 
    # 4. Each bin related to each node, needs all the constraint matrices
    #    from all bins in the other nodes in the toroidal direction. 
    #    
    #    ^^^THIS STATEMENT IS MISLEADING. ELABORATION:
    #       
    #       Only those bins, on each node (which share an R,Z coordinate), 
    #       which share a range of v_para and mu 
    #       need to share their constraint matrices. 
    #

    #
    # First, define a new communicator for an axisymmetric grouping 
    #
    # Don't know if this should be done in the loop.
    # It's possible that the process is assigned multiple nodes
    # which are dispersed across multiple planes. 
    #

    if (inode2 - inode1) > nnodes_plane:
        raise ValueError("Process handling nodes across planes, edge case.")
    
    axi_color = np.int(0)
    axi_comm = sml_comm.Split(color=axi_color, key=procno)

    for bin_vp in range(0, vpsize):
        for bin_mu in range(0, musize):
            for nnode in range(inode1, inode2+1):
                # Sum all constraint vectors
                axi_comm.Barrier()
                axi_constraint_vec = np.zeros(shape=(5,),)
                cvec = bins[nnode,bin_vp,bin_mu]

                axi_constraint_vec = axi_comm.Allreduce()
                
                #
                
                #



    # Need to know beforehand how many particles are in each bin that's being
    # gathered in order to allocate the right size data structure for each 
    # processing which will build the axisymmetric constraint matrix

    #for bin_vp in range(0, vpsize):
    #    for bin_mu in range(0, musize):
    #        for nnode in range(inode1, inode2+1):
    #            data = sml_comm.gather(bins[nnode, bin_vp, bin_mu], 
    #                root=(procno % nnodes_plane))


    # Figure out how to gather it all on the root process
    data = bin_.constraint_mat
    data = sml_comm.gather(data, root=0)
    tmp = np.empty([5, 1])
    if procno == 0:
        for array in data:
            tmp = np.hstack([tmp, array])
        #print(tmp)



    #
    # Solve optimization
    #

#    soln = np.linalg.lstsq(full_matrix, bin_.constraints, rcond=None)

    #
    # Distribute weight 
    #


    # MPI Finalize not needed
    return 0


#def old():


 #   bin_.add_particle(Particle(weight=w, velocity=v))

#    # make sure the constraint vector is updated
#    bin_.update_constraints()
#    bin_.build_constraint_mat()
    ## with the appropriate size. 

#    data = bin_.constraint_mat
#    data = axi_comm.gather(data, root=0)
#    tmp = np.empty([3, 1])
#    if local_rank == 0:
#        for array in data:
#            tmp = np.hstack([tmp, array])
        #print(tmp)

    # broadcast back to everybody in communicator
#    if local_rank == 0:
#        full_matrix = np.array(tmp)
#    else:
#        full_matrix = None

#    full_matrix = axi_comm.bcast(full_matrix, root=0)
