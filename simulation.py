# file: simulation.py 
import numpy as np
from mpi4py import MPI
from data_structures.sim_structures import *


def sim(num_nodes=20, vp_max=1.0, mu_max=1.0, npart_max=3, nnodes_plane=4):
    sml_comm = MPI.COMM_WORLD
    num_procs = sml_comm.Get_size()
    procno = sml_comm.Get_rank()

    #
    # build grid on main process, send to all processes.
    #
    if procno == 0:
        grid = []
        for i in range(0,num_nodes):
            particles = []
            for npart in range(1, np.int(np.random.randint(npart_max))+1):
                w = np.float64(np.random.weibull(a=np.random.randint(0,7)))
                vp = np.float64(np.random.uniform(low=1e-5, high=vp_max))
                mu = np.float64(np.random.uniform(low=1e-5, high=mu_max))

                particles.append(Particle(weight=w, vp=vp, mu=mu))

            grid.append([i, # node_id
                         particles]) 
    else:
        grid = None


    print("Proc ", procno, "hitting barrier")
    sml_comm.Barrier()

    grid = sml_comm.bcast(grid, root=0)

    print("Proc", procno, "has grid", grid)

    sml_comm.Barier()

    start_node = 0
    end_node = 1

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

    #print("Process", global_rank, "reports recv", full_matrix.shape)

    ## Now that everybody in the communicator has the axisymmetric constraint
    ## matrix, time to solve for the new weights with it's own constraint
    ## vector
#    soln = np.linalg.lstsq(full_matrix, bin_.constraints, rcond=None)

#    print("Process", global_rank, "now has solution", soln)


#    return None
    
