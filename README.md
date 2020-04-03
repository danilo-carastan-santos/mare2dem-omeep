# mare2dem-omeep
MARE2DEM source files customized for the OMEEP project. The original source code can be found at https://mare2dem.ucsd.edu/?page_id=108

To compile, run the command line make CLUSTER=omeep
INTEL_PATH=<intel_path> in the repository root dorectory, where <intel_path> is the path to the "intel" directory (where the compiler are installed, it is often /opt/intel).

From the original source code download page (https://mare2dem.ucsd.edu/?page_id=108) there are some data set examples. To perform a test execution, go to the inversion_CSEM directory (from the examples root directory) and run the command mpirun -n <nb_processes> <path>/MARE2DEM Demo.0.resistivity, where <path> is the path that leads to the MARE2DEM binary, and <nb_processes> is the number of MPI processes. At least two MPI processes are needed to run the application.  
 
The makefile compilation commands have the scorep instrumentation already set. More information about scorep can be found at the link http://www.vi-hps.org/projects/score-p/
