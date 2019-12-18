############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          SLmake.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   February 15, 2000
#
#  Modified:        October 13, 2011
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################
#
#  C preprocessor definitions:  set CDEFS to one of the following:
#
#     -DNoChange (fortran subprogram names are lower case without any suffix)
#     -DUpCase   (fortran subprogram names are upper case without any suffix)
#     -DAdd_     (fortran subprogram names are lower case with "_" appended)

CDEFS         = -DAdd_

#
#  The fortran and C compilers, loaders, and their flags
#

#FC            = mpif90
#CC            = mpicc 
NOOPT         = -O0
#FCFLAGS       = -O3
#CCFLAGS       = -O3


    FC      = /usr/local/bin/mpif90  
   FCFLAGS  =  -O2 -m64  -fpp -i8   #-O2 -m64 -warn all -check bounds
   CC      = /usr/local/bin/mpicc  
   CCFLAGS  = -O2 -m64 

FCLOADER      = $(FC)
CCLOADER      = $(CC)
FCLOADFLAGS   = $(FCFLAGS)
CCLOADFLAGS   = $(CCFLAGS)

#
#  The archiver and the flag(s) to use when building archive (library)
#  Also the ranlib routine.  If your system has no ranlib, set RANLIB = echo
#


#ARCH          = ar
#ARCHFLAGS     = cr
#RANLIB        = ranlib
   ARCH = /opt/intel/Compiler/11.1/084/bin/intel64/xiar  # use this with the intel icc compiler and optimization -O2 or faster
   ARCHFLAGS = ruv
   RANLIB = ranlib
   
#
#  The name of the ScaLAPACK library to be created
#

SCALAPACKLIB  = libscalapack.a

#
#  BLAS, LAPACK (and possibly other) libraries needed for linking test programs
#

#BLASLIB       = -lblas
#LAPACKLIB     = -llapack
#LIBS          = $(LAPACKLIB) $(BLASLIB)
    
    MKLPATH=$(MKLROOT)/lib/em64t
    MKLINCLUDE=$(MKLROOT)/include
    BLASLIB = -L$(MKLPATH) -I$(MKLINCLUDE)  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core 
    LIBS    = $(BLASLIB) 
 

