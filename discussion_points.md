# Discussion Points

## Communicator Implementation 

> Note: Scroll down to bolded section that talks about the 3 options for the real-real. 

> Note: I only found out about "Parallel Programming in MPI and OpenMP" by Victor Eijkhout from TACC recently, it was the text for Boulder's open topics in high performance scientific computing . I am still in the process of reading it, but am using the HTML version as a reference text in addition to online examples. 

The general form for ranks, groups, and communicators is to 
1. Construct a group reference from ranks in group for COMM_WORLD that is a subset of the Group for COMM_WORLD 
2. Create a communicator from that group
3. Use that communicator for all further processing for that logical group. (You have to provide a communicator handle for any collective operation in MPI, ) 

The reasons for this are:
* while thread safe, operations on communicators can become blocking for all 
  processes associated with that communicator, and this is mitigated by using 
  Groups. None of the collective operations are interrupt-safe, at least in the MPICH implementation.  
* MPICH implementation says MPI_Comm_split does an ALL_GATHER on the entire 
  parent communicator to find which process IDs have the specified color. 
* ALL_GATHER has an implicit MPI_Barrier and synchronization prior to and 
  directly after communication to fill the buffer and determine the correct 
  communicator for each color provided to it. 
* This means that not only is there intense I/O for splitting communicators 
  (already drainingly slow in cluster applications) in a scenario with 
  hundreds of thousands of processes, but the entirety of COMM_WORLD is forced 
  to synchronize for each different "color" for a new communicator during that 
  'setup' phase.
* Group calculation can occur locally on each thread, and communicator 
  creation (using the constructor that takes an MPI_Group of ranks) can occur 
  asynchronously. 
* The groups can be persisted outside of the time-evolution of the simulation, 
  therefore axisymmetric group & related communicator setup & tear-down only 
  have a one-time cost. While this is true of communicators as well, it may 
  require changes in sml_module to pass references to every process of its 
  axisymmetric communicator, 
* Usage of groups means severe reduction in the amount of I/O that occurs, in 
  addition to having asynchronous startup. 

__So far there are 3 options that I've come up with so far for figuring out how to take a processes' global rank and using that to define an axisymmetric communicator.__

Defining an additional axisymmetric communicator may be necessary in MPI to reduce complexity and maintain sanity of the code, depending on how the constraint matrices need to be shared to individual processes. 
It may be possible to devise a scheme wherein all processes are able to take their simulation communicator rank and calculate the rank of the process it needs to send its information to, but this means it is possible that each process may have up to N^2 blocking send/recvs during the communication step, where N is the number of φ-planes. 
Additionally, it is organizationally difficult to achieve this, making the code error-prone during development, and may take additional resources to ensure effectiveness. 

It is recommended that the communicator creation occurs prior to resampling (or potentially outside resampling altogether) during simulation start-up, after the mesh grid has been defined, since the number of bins representing spatial (φ,r,z) bins does not change across iterations of the time-evolution. Then, have each bin maintain a reference to its axisymmetric group outside of the resampling procedure.  
This allows for the communicator handle to be passed to the process at arbitrary times during execution. The following assume that the process only knows its rank during execution. 

1. MPI_Comm_split 

If the process global rank (or the simulation MPI communicator) is isomorphic to a single 5d phase space bin, then it's possible to assign a "color" (communicator designation or ID) to each process. 
Another way to think of this is a "color" assignment based on (r,z) coordinates for each process. However, this color must be represented as an integer for MPI, and I haven't thought of a good way to map unique (r,z) coordinates to a single integer. 
One way to accomplish it is through division and modulus of those coordinates (after division and modulus on the rank to assign them (r,z) coords in the first place, if these values are not already known by the process). 
Another way is to do something like 
```
	color = r + (rank % num_z_vals) * z
```
or, since this is analogous to finding a length of a unique value z=(x,y) in a cartesian graph, ```color = r^2 + z^2```, assuming z=0 is actually underneath the volume in question (r,z > 0). 
Assuming the z-plane bisects the volume (z positive or negative), this could be altered to ```color = r^2 + m * (z)^2```. Where m=1 if z>0, m=2 if z<0.  

From there, each process synchronizes on the call to MPI_Comm_split, and fills its reference to a new communicator that is supplied by the MPI instance.

What is currently unclear is if multiple "colors" supplied as an argument also causes additional synchronizations. (I don't suspect this is the case.) However, the amount of time MPI will spend building the communicator increases with the number of created axisymmetric communicators. 
Following that, the usual scatter, gather, and other operations are then available on this communicator to share the constraint matrices. 

__Pro__
* Comparatively quick/easy to implement. 
* Easy to quickly drop in to verify that the rest of the program is working correctly with the axisymmetric resampling and constraint-based quadratic optimization.
* Easy, small implementation allows this to be dropped into resamp_mod and 
  later replaced or modified over time with a more efficient or fault-tolerant implementation (not relying on process rank for figuring out the axisymmetric communicator, not putting the communication structure inside the resampling code, etc.)

__Con__
* Can force entire process structure to synchronize at least once [0]. If care is not taken and depending on where this is called, it could occur during _every_ iteration of the time-evolution, or in a worst case scenario, during each resampling. 
* Variable velocity space grid, variable spatial mesh topology can cause a headache for using the process rank to calculate what (r,z) value the process is associated with. 
* This implementation may not be sustainable or easily maintainable with changes to the rest of XGC over time. 
* If not handled carefully and placed in an appropriate location, may cause serialization to occur. 
* May require intensive I/O between all the processes (since they're distributed across compute nodes). 
* Requires additional knowledge from sml_module regarding total number of φ-planes, total number of vertices per φ-plane, size of v-space grid. May require additional calls to MPI instance to find out if the number of MPI objects is approaching the limit. 

__Status Quo__
* Still requires some number of axisymmetric communicators 
* Memory load per process still large to hold the axisymmetric constraint matrix. 
* Still going to have issues with how to share the constraint matrices given that they may have variable column dimensions. (Addressed in later section)

__Sources__

- [0] https://www.mpich.org/static/docs/latest/www3/MPI_Comm_split.html

- [1] http://pages.tacc.utexas.edu/~eijkhout/pcse/html/mpi-comm.html#Splittingacommunicator


2. Relying on MPI_Group, inter- and intra-communications

Each process may independently determine its axisymmetric Group, or the process ranks it shares its constraint matrix info with. 
Then, have each process create/lookup the communicator for that Group. 

>I am unsure if each process holds a unique communicator or a reference to a communicator held by the MPI implementation; or if MPI implementations account for not creating a new communicator for an identical group (to avoid unnecessary copies).
>
> I believe that MPI does _not_ create per-process copies of the communicator when created from a group. It wouldn't make sense that per-process copies would exist and you could still do an Allgather on them. 

There are two mutually exclusive assumptions that can be made, which affect if a Group is created and how the communicator is generated. 

___

(2a) There are no pre-existing groups which must be relied on. 

Using the process rank and information about number of φ-planes, spatial grid topology, velocity phase space numbering, etc., determine the array of ranks which this process must communicate its constraint matrix to. 


(2b) There is an existing group structure which should be relied on, additional groups are not preferred. (MPI object limit, memory constraints)

The way to accomplish this is to define an inter-communicator between two groups, where an existing group in this case may represent a φ-plane.

___

Then, the communicator may share the constraint matrices (implementation discussed in later section). 


__Pro__
* Non-blocking communication w.r.t. the simulation COMM_WORLD. Will not (_necessarily_) force synchronization across COMM_WORLD, while creating the communicator. The communicator may or may not be duplicated per process, but either way is nonblocking. 
* Reduces some I/O overhead compared to MPI_Comm_split 

__Con__
* As simulation resolution increases, the initial memory load per process of keeping track of all the processes for the axisymmetric communicator may be very large (In the case of . 
* For (2b), only blocking process-to-process communication is allowed. However, this only blocks the two processes sharing information with each other. 

__Status Quo__
* For (2a), have to use a collection of functions like MPI_Group_incl and exclusion function to modify the group. 

__Sources__

- http://pages.tacc.utexas.edu/~eijkhout/pcse/html/mpi-comm.html#Inter-communicators

- https://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-1.1/node111.htm#Node111

- https://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-1.1/node114.htm


3. MPI Virtual Topology (distributed Graph or Cartesian Grid)

The process keeps 

__Pro__

__Con__

__Status Quo__

__Sources__

- https://www.mcs.anl.gov/~balaji/pubs/2011/eurompi/eurompi11.mpi-topology.pdf

- http://pages.tacc.utexas.edu/~eijkhout/pcse/html/mpi-topo.html

- https://www3.nd.edu/~zxu2/acms60212-40212/Lec-08.pdf


## Algorithm for sharing the constraint matrix, then performing the constraint optimization. 

### Sharing the constraint matrix

1. Use collective communication methods. 

* __[a]__ Allgather to all processes in group so they may solve the constraint matrix individually

* __[b]__ Gather all matrices to a single process per axisymmetric group 

2. Use MPI data types, packing, and defined structs. 

3. Use shared memory regions

One way to share it may be by creating a shared memory location (referred to as an MPI_Window) using MPI_Win_create. 

https://www.mpich.org/static/docs/v3.2.x/www3/MPI_Win_create.html

4. MPI has a facility for caching data through attachment to a communicator. 

https://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-1.1/node118.htm#Node118

__Sources__

## Further Recommendations/Ideas

1. Enforce communications structure in a separate module called from sml_module, or within sml_module itself. 

2. Make each MPI process tied only to a spatial mesh grid coordinate, (φ,r,z). Have each process call MPI_Spawn to generate worker processes representing the velocity phase plane or grid and use the created worker-process communicator to handle resampling. 

## Questions about original implementation

1. XGCa already exists. Where should I be looking inside that codebase for how it handles axisymmetric communications? What is it already communicating axisymmetrically?

2. I found the grid_type, which already has a list of all (R,Z) coordinates for the system.  



