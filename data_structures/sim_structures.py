# file: sim_structures.py
import numpy as np

"""
Note: These currently exist as dataclasses w/o the attribute to serve
as respresentations prior to working out the idea further
"""

class Particle():
    def __init__(self, weight, velocity):
        self.weight = weight # we'll make this an integer from 1-100
        self.velocity = velocity # this should be np.float64


class Bin():
    def __init__(self, r, z, v, plane_id):
        self.r = r # the spatial coordinate in terms of radial distance from
                   # center of torus
        self.z = z
        self.v = (v[0], v[1]) # lower, upper bounds of velocity space
                              # represented by this bin
        self.plane_id = plane_id # φindex of the plane this bin is tied to
        
        self.particles = []   # Populate particles per bin during the sim

        self.constraints = np.array([0,0,0]) # build this with values after particles
        self.constraint_mat = None # placeholder, build after particles added


    def add_particle(self, particle):
        self.particles.append(particle)


    def update_particles(self, weights):
        particles = self.particles
        for i in len(weights):
            particles[i].weight = weight[i]

    
    def update_constraints(self):
        self.constraints[0] = sum(p.weight for p in self.particles)
        self.constraints[1] = sum(w*v for w,v in zip([p.weight for p in self.particles], [p.velocity for p in self.particles]))
        self.constraints[2] = sum(w*(v**2) for w,v in zip([p.weight for p in self.particles], [p.velocity for p in self.particles])) 


    def build_constraint_mat(self):
        tmp = [[1, p.velocity, p.velocity**2] for p in self.particles]
        self.mat = np.stack(tmp, axis=-1)



class Plane():
    def __init__(self, id, bins=[]):
        self.id = id
        self.bins = bins
