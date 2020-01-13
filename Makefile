# 
# Makefile for the MARE2DEM Parallel MPI 2.5D EM Modeling code
# 
# Kerry Key
# Scripps Institution of Oceanography
# kkey@ucsd.edu
# 
# Feb   2013    Added ScaLAPACK library
# 
# March 2012    Minor tweaks for code release to SEMC
#
# April 2010    Added support for Triton compute cluster 
# 
# March 2010	After reading the first four chapters of the O'Reilly GNU Make book...
#               Redesigned the Makefile so that it will also make Metis and SuperLU
#               by overriding the default variables in their respective Makefiles. 
#               In other words, you only need to set your system's specific settings 
#               here, rather than in the Metis and SuperLU makefiles as well.  
# 
# 
# 
#  Usage:   make CLUSTER=<ClusterName>      where <ClusterName> is from the list below and 
#                                           sets the corresponding options for your system.
#                                             
# 
#---------------------------------------------------------------------
# Specify your  MPI Compilers and arguments:
#----------------------------------------------------------------------

CLUSTER_LC := $(shell echo $(CLUSTER) | tr A-Z a-z)
# 
# macpro with Intel Compilers:
# 
ifeq "$(CLUSTER_LC)" "omeep"
   FC_NO_INSTR      = /opt/intel/compilers_and_libraries/linux/mpi/intel64/bin/mpif90
   FC      = $(PREP) $(FC_NO_INSTR)
   FFLAGS  = -O2 -debug all -m64  -fc=/opt/intel/compilers_and_libraries/linux/bin/intel64/ifort# optimized code  
#   FFLAGS  = -O2 -m64  -fpp # optimized code
#   FFLAGS =  -fpp  -m64 -warn all -check  -traceback -fpe0   -stand f03 -fstack-security-check # for debugging code
   CC_NO_INSTR      = /opt/intel/compilers_and_libraries/linux/mpi/intel64/bin/mpicc
#   CC      = $(PREP) $(CC_NO_INSTR)
   CC      = $(CC_NO_INSTR)
   CFLAGS  = -O2 -m64 
   # You only need these commands if you are compiling the Metis and SuperLU libs in MARE2DEM/Source/SuperLU and /Metis:
   ARCH = /opt/intel/compilers_and_libraries/linux/bin/intel64/xiar  # use this with the intel icc compiler and optimization -O2 or faster
   ARCHFLAGS = ruv
   RANLIB = ranlib
   BLASDEF = -DUSE_VENDOR_BLAS  
   
# Use this for the threaded (openmp) version of lapack and blas in the Intel Math Kernel Library:
    MKLPATH=$(MKLROOT)/lib
    MKLINCLUDE=$(MKLROOT)/include
    BLASLIB =  -L$(MKLPATH) -I$(MKLINCLUDE)  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
 
    SUPERLU_CDEFS = -DAdd_    #  -DAdd_ or -DNoChange, this is for the stupid underscore difference between C and Fortran, 
                                 # use nm <objectfile.o> to see how the symbol endings in the c and fortran files 
endif

# 
#  Dual:
# 
ifeq "$(CLUSTER_LC)" "dual"
   FC=/vend/intel/impi/4.0.1.007/bin64/mpif90
   FFLAGS = -O2 -fpp -fc=/vend/intel/composerxe-2011.2.137/bin/intel64/ifort
   CC      = /vend/intel/impi/4.0.1.007/bin64/mpicc
   CFLAGS  = -O2 
   # You only need these commands if you are compiling the Metis and SuperLU libs in MARE2DEM/Source/SuperLU and /Metis:
   ARCH = /vend/intel/composerxe-2011.2.137/bin/intel64/xiar  # use this with the intel icc compiler and optimization -O2 or faster
   ARCHFLAGS = ruv
   RANLIB = ranlib
   BLASDEF = -DUSE_VENDOR_BLAS  
   MKLPATH=$(MKLROOT)/lib/intel64
   MKLINCLUDE=$(MKLROOT)/include
   BLASLIB =  -L$(MKLPATH) -I$(MKLINCLUDE)  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
   
   SUPERLU_CDEFS = -DAdd_    #  -DAdd_ or -DNoChange, this is for the stupid underscore difference between C and Fortran, 
                                 # use nm <objectfile.o> to see how the symbol endings in the c and fortran files 
endif

# 
#  Lonestar TACC Dell Linux cluster:
# 
ifeq "$(CLUSTER_LC)" "lonestar"
   FC      = mpif90  
   FFLAGS  =  -O2 -fpp -xT -m64
   CC      = mpicc  
   CFLAGS  =  -O2 -xT -m64
   # You only need these commands if you are compiling the Metis and SuperLU libs in MARE2DEM/Source/SuperLU and /Metis:
   ARCH = xiar # use this with the intel icc compiler and optimization -O2 or faster
   ARCHFLAGS = ruv
   RANLIB = ranlib
   BLASDEF = -DUSE_VENDOR_BLAS  
   BLASLIB = -Wl,-rpath,$(TACC_MKL_LIB) -I$(TACC_MKL_INC) -L$(TACC_MKL_LIB) \
	 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread


   SUPERLU_CDEFS =  -DAdd_   #  -DAdd_ or -DNoChange, this is for the stupid underscore difference between C and Fortran, 
                                 # use nm <objectfile.o> to see how the symbol endings in the c and fortran files 
endif


#
# Triton Compute Cluster at SDSC:
#
# Using intel and openmpi. NOTE that you need to modify your shell to load modules intel and openmpi_mx.
#
ifeq "$(CLUSTER_LC)" "triton"
   FC	   = mpif90
   FFLAGS  = -O2 -fpp
   CC	   = mpicc
   CFLAGS  = -O2
   # You only need these commands if you are compiling the Metis and SuperLU libs in MARE2DEM/Source/SuperLU and /Metis:
   ARCH = xiar
   ARCHFLAGS = ruv
   RANLIB = ranlib
   BLASDEF = -DUSE_VENDOR_BLAS
   MKLPATH=${MKL_ROOT}/lib/em64t/
   BLASLIB =  -L$(MKLPATH)    -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
    
   SUPERLU_CDEFS = -DAdd_    #  -DAdd_ or -DNoChange, this is for the stupid underscore difference between C and Fortran
                                 # use nm <objectfile.o> to see how the symbol endings in the c and fortran files
endif

#
# TSCC: Triton Shared Compute Cluster (New Spring 2013). Replaces "triton"
#
# Using intel and openmpi.  
#
ifeq "$(CLUSTER_LC)" "tscc" 
   FC	   = mpif90
   FFLAGS  = -O2 -fpp  -m64  
   CC	   = mpicc
   CFLAGS  = -O2
   # You only need these commands if you are compiling the Metis and SuperLU libs in MARE2DEM/Source/SuperLU and /Metis:
   ARCH = xiar
   ARCHFLAGS = ruv
   RANLIB = ranlib
   BLASDEF = -DUSE_VENDOR_BLAS
   MKLPATH=${MKL_ROOT}/lib/intel64/  
   BLASLIB =  -L$(MKLPATH)    -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
    
   SUPERLU_CDEFS = -DAdd_    #  -DAdd_ or -DNoChange, this is for the stupid underscore difference between C and Fortran
                                 # use nm <objectfile.o> to see how the symbol endings in the c and fortran files
endif

 
#---------------------------------
# clusters beneath this line have NOT been updated for the March 2012 code updates. The settings below may need to be 
# tweaked slightly.


# Intel error checking: use these for debugging source code:
#  FFLAGS =  -m64 -stand f03 -warn all -fstack-security-check -check all     

# gfortran error checking:
# FFLAGS = -m64  -Wall -g -fbounds-check -Wuninitialized -O -ftrapv -fimplicit-none -fno-automatic



 
#-------------------------------------------------------------------------------
# Libraries required by MARE2DEM:
#-------------------------------------------------------------------------------
 
# 
# SuperLU library:
# 
SUPERLU_dir     = ./libraries/SuperLU_4.0
LIBSUPERLU      = $(SUPERLU_dir)/libsuperlu_4.0.a
SUPERLU_HEADER = -I$(SUPERLU_dir)/SRC

# 
# Metis library:
# 
ARMETIS ='$(ARCH) $(ARCHFLAGS)'  # Metis AR needs the flags attached

METIS_dir  =  ./libraries/metis-4.0
LIBMETIS   = $(METIS_dir)/libmetis.a

# 
# SuiteSparse library for UMFPACK:
#  
Suitesparse =./libraries/SuiteSparse

UMFdir  = $(Suitesparse)/UMFPACK

LibUMF  =   $(UMFdir)/Lib/libumfpack.a
LibAMD  =   $(Suitesparse)/AMD/Lib/libamd.a 
LibCAMD =   $(Suitesparse)/CAMD/Lib/libcamd.a 
LibCHLM =   $(Suitesparse)/CHOLMOD/Lib/libcholmod.a
LibCLMD =   $(Suitesparse)/COLAMD/Lib/libcolamd.a
LibCCMD =   $(Suitesparse)/CCOLAMD/Lib/libccolamd.a

#LibMETS =   $(METIS_dir)/libmetis.a
 
LIBUMF  = $(LibUMF) $(LibAMD) $(LibCHLM)  $(LibCLMD) $(LibCAMD) $(LibCCMD)   #$(LIBMETIS)
#LIBUMF  = $(LibUMF) $(LibAMD) $(LibCAMD) $(LibCLMD) $(LibCCMD) $(LibCHLM)  #$(LIBMETIS)
ARUMF   = $(ARCH) $(ARCHFLAGS)
IncUMF  = -I$(Suitesparse)/UMFPACK/Include  \
          -I$(Suitesparse)/UFconfig \
          -I$(Suitesparse)/AMD/Include  \
          -I$(Suitesparse)/CCOLAMD/Include\
          -I$(Suitesparse)/CAMD/Include\
          -I$(Suitesparse)/CHOLMOD/Include \
          -I$(Suitesparse)/COLAMD/Include   
        
# 
# ScaLAPACK Library"
#
SCALAPACK_dir =./libraries/scalapack-2.0.2
LIBSCALAPACK  = $(SCALAPACK_dir)/libscalapack.a


#-------------------------------------------------------------------------------
# Compilation Commands
# 
# You shouldn't need to muck with anything beneath here
# 
#-------------------------------------------------------------------------------

#TARGETS: clean MARE2DEM  TestForSlivers  
 

ifndef CLUSTER
    MARE2DEM: checkCLUSTERarg
endif
  
ifeq "$(CLUSTER)" "lonestar"
    MARE2DEM: buildMARE2DEM
else ifeq "$(CLUSTER)" "omeep"
    MARE2DEM: buildMARE2DEM
else ifeq "$(CLUSTER)" "triton"
    MARE2DEM: buildMARE2DEM  
else ifeq "$(CLUSTER)" "tscc"
    MARE2DEM: buildMARE2DEM      
else ifeq "$(CLUSTER)" "dual"
    MARE2DEM: buildMARE2DEM
else
    MARE2DEM:  checkCLUSTERarg    
endif

 

# 
# Error message if user doesn't specify cluster variable on input, i.e.,: make CLUSTER=lonestar 
#
checkCLUSTERarg: 
	@printf "\n\n\n !!!!!!!!! Error making MARE2DEM !!!!!!!!! \n\n";
	@printf "\n    CLUSTER variable is incorrect or undefined    \n\n";
	@printf "    Usage: make CLUSTER=<myclustername>   \n\n";
	@printf "    Currently supported clusters: \n\n"; 
	@printf "        $ make CLUSTER=tscc    \n";
	@printf "        $ make CLUSTER=omeep  \n";
	@printf "        $ make CLUSTER=lonestar    \n";
	@printf "        $ make CLUSTER=triton    \n";
	@printf "        $ make CLUSTER=dual      \n";
	@printf "         \n\n\n\n\n";	
# kwk debug: still need to check that CLUSTER variable is an acceptable value...

all: $(TARGETS)

# 
#  Cleaning functions:
# 
clean:  clean_mare2dem clean_metis clean_superlu clean_umfpack clean_scalapack
		
clean_mare2dem:
	@printf "#\n#\n# Cleaning MARE2DEM \n#\n#\n"; 
	rm -f *.o *.mod  MARE2DEM;

clean_metis:
	@printf "#\n#\n# Cleaning Metis Graph Partitioning Library: \n#\n#\n"; \
	cd $(METIS_dir); make clean; rm -f *.a;

clean_superlu:
	@printf "#\n#\n# Cleaning SuperLU Sparse Linear Solver Library: \n#\n#\n"; \
	cd $(SUPERLU_dir); pwd; make clean; 
	
clean_umfpack:
	@printf "#\n#\n# Cleaning UMFPack Sparse Linear Solver Library: \n#\n#\n"; \
	cd $(UMFdir); pwd; make purge  	
	
clean_scalapack:
	@printf "#\n#\n# Cleaning ScaLAPACK Library: \n#\n#\n"; \
	cd $(SCALAPACK_dir); pwd; make clean;  	
		
# 
# SuperLU Library build:
# 
# Note that these are only executed if libsuperlu.a and libmetis.a can't be found:
# For example, this is not called if $(LIBSUPERLU) points to
# your cluster's own superlu library that is outside MARE2DEM/Source
# 
# 
$(LIBSUPERLU): 
	@printf "#\n#\n# Making SuperLU Sparse Linear Solver Library: \n#\n#\n"; \
	cd $(SUPERLU_dir); \
	make superlulib CC=$(CC_NO_INSTR) CFLAGS='$(CFLAGS)' \
	FORTRAN=$(FC_NO_INSTR) FFLAGS='$(FFLAGS)' CDEFS=$(SUPERLU_CDEFS) BLASDEF=$(BLASDEF) \
	ARCH=$(ARCH) ARCHFLAGS=$(ARCHFLAGS) RANLIB=$(RANLIB); cd $(CURDIR);

# 
#  Metis Library build:
# 
$(LIBMETIS):
	@printf "#\n#\n# Making Metis Graph Partitioning Library: \n#\n#\n"; \
	cd $(METIS_dir); \
	make default CC=$(CC_NO_INSTR) OPTFLAGS='$(CFLAGS)' AR=$(ARMETIS)  RANLIB=$(RANLIB); cd $(CURDIR);

# 
# Suitesparse UMFPack Build:
# 
$(LIBUMF):  
	@printf "#\n#\n# Making UMFPack Sparse Linear Solver Library: \n#\n#\n"; \
	cd $(UMFdir);\
	make CC=$(CC_NO_INSTR) CFLAGS='$(CFLAGS)' \
	F77=$(FC_NO_INSTR) F77FLAGS='$(FFLAGS)'  BLAS='$(BLASLIB)' \
	AR='$(ARUMF)'  RANLIB=$(RANLIB) METIS=$(PWD)/$(LIBMETIS) METIS_PATH=$(PWD)/$(METIS_dir);\
	cd $(CURDIR);

#
# ScaLAPACK library build:
#
$(LIBSCALAPACK):
	@printf "#\n#\n# Making ScaLAPACK Library: \n#\n#\n"; \
	cd $(SCALAPACK_dir); \
	make lib CC=$(CC_NO_INSTR) CCFLAGS='$(CFLAGS)' \
	FC=$(FC_NO_INSTR) FCFLAGS='$(FFLAGS)'  CDEFS=$(SUPERLU_CDEFS) \
	ARCH=$(ARCH) ARCHFLAGS=$(ARCHFLAGS) RANLIB=$(RANLIB); cd $(CURDIR);


# 
# MARE2DEM build:
# 

MARE2DEM_Files	=  kdtree2.o fem2D_utilities.o  binaryTree.o  call_triangle.o sort.o \
			  c_fortran_zgssv.o  superlu_zsolver.o \
			  umf4_f77zwrapper.o umfpack_zsolver.o \
			  string_helpers.o triangle.o mt1D.o EMconstants.o em2dkx.o   Occam.o\
		      c_fortran_triangle.o FilterModules.o   \
		      mare2dem_common.o  spline_kx_module.o mare2dem_worker.o mare2dem_io.o \
		      mare2dem_mpi.o   EM2D.o RunMARE2DEM.o 

buildMARE2DEM:	$(LIBSUPERLU) $(LIBMETIS) $(LIBUMF) $(LIBSCALAPACK) $(MARE2DEM_Files) 
	        $(FC) $(FFLAGS) $(MARE2DEM_Files) $(LIBSUPERLU) $(LIBUMF) $(LIBMETIS) $(LIBSCALAPACK) \
	        $(BLASLIB) -o MARE2DEM

# 
# Test build:
# 
Test_Files	= binaryTree.o call_triangle.o triangle.o  c_fortran_triangle.o  \
		      TestForSlivers.o  
		
TestForSlivers:	 $(Test_Files) 
	$(FC) $(FFLAGS) $(Test_Files)   -o TestForSlivers
	
# 
# Triangle build:
# 
# Don't optimize Triangle on the G5 since unexpected calamities will happen:
# its way faster than the 2D EM solves and doesn't need optimization anyway
# 
TRILIBDEFS = -DTRILIBRARY   # Triangle compile flag, don't change this!

triangle.o:  triangle.c  triangle.h
	$(CC)  $(TRILIBDEFS) -c -o $(BIN)triangle.o triangle.c
		
# 
# Pattern matching rules for Fortran and C files (these replace old style .SUFFIXES rules):
# 
# General Fortran compile:
%.o: %.f90 
	$(FC) $(FFLAGS)    -c -o $@ $^
	
# General C compile:
%.o : %.c
	$(CC) $(CFLAGS) $(CDEFS) $(SUPERLU_HEADER) $(IncUMF) -c -o $@ $< 


