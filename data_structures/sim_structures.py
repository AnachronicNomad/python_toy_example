# file: sim_structures.py
import numpy as np

"""
Note: These currently exist as dataclasses w/o the attribute to serve
as respresentations prior to working out the idea further
"""

class Particle():
    def __init__(self, weight, vp, mu):
        self.weight = np.float32(weight)
        self.vp = np.float64(vp)
        self.mu = np.float64(mu)


class Bin():
    def __init__(self):
        self.npart = np.int32(0)
        self.ptl_inx = []
        self.total_weight = np.int32()

        self.constraints = np.zeros((5,), dtype=np.float32) 
        self.constraint_mat = None 


    def add_particle(self, particle):
        self.particles.append(particle)


    def update_particles(self, weights):
        particles = self.particles
        for i in len(weights):
            particles[i].weight = weight[i]

    
    def update_constraints(self):
        self.constraints[0] = sum(p.weight for p in self.particles)
        self.constraints[1] = sum(w*vp for w,vp in zip([p.weight for p in self.particles], [p.vp for p in self.particles]))
        self.constraints[2] = sum(w*(vp**2) for w,vp in zip([p.weight for p in self.particles], [p.vp for p in self.particles])) 
        self.constraints[3] = sum(w*mu for w,mu in zip([p.weight for p in self.particles], [p.mu for p in self.particles]))
        self.constraints[4] = sum(w*(mu**2) for w,mu in zip([p.weight for p in self.particles], [p.mu for p in self.particles])) 

    def build_constraint_mat(self):
        tmp = [[1.0, p.vp, p.vp**2, p.mu, p.mu**2] for p in self.particles]
        self.constraint_mat = np.stack(tmp, axis=-1)



class Plane():
    def __init__(self, id, bins=[]):
        self.id = id
        self.bins = bins
