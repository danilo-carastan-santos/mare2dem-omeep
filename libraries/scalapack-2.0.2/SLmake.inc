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

FC            = mpif90
CC            = mpicc 
NOOPT         = -O0
FCFLAGS       = -O3
CCFLAGS       = -O3


 
NOOPT = -O0
#FCFLAGS = -O3
#CCFLAGS = -O3
FCLOADER = $(FC)
CCLOADER = $(CC)
FCLOADFLAGS = $(FCFLAGS)
CCLOADFLAGS = $(CCFLAGS)
ARCH = ar
ARCHFLAGS = cr
RANLIB = ranlib
SCALAPACKLIB = libscalapack.a
LIBS=-lmkl_intel_ilp64 -lmkl_sequential -lmkl_core 
