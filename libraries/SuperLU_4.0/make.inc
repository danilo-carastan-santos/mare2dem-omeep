############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   October 2, 1995
#
#  Modified:	    February 4, 1997  Version 1.0
#		    November 15, 1997 Version 1.1
#		    September 1, 1999 Version 2.0
#
############################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _osx  #_linux

#
#  The name of the libraries to be created/linked to
#
SuperLUroot	= ../   #$(HOME)/Codes/SuperLU_4.0
SUPERLULIB   	= ../libsuperlu_4.0.a  #$(SuperLUroot)/lib/libsuperlu_4.0.a

#BLASLIB   	= $(SuperLUroot)/lib/libblas.a
BLASDEF 	= -DUSE_VENDOR_BLAS

## ATLAS BLAS causes single-precision to fail
#BLASLIB   	= -L/usr/lib/atlas -lblas

## This BLAS causes single-precision to fail the test in SuperLU
#BLASLIB 	= -L/usr/lib -lblas

# BLASLIB 	= -L$(HOME)/lib/GotoBLAS -lgoto

# OS X developers tools Blas:
#BLASLIB 	=  -framework vecLib

# MKL BLAS:
#MKLPATH=$MKLROOT/lib/em64t
#MKLINCLUDE=$MKLROOT/include

#BLASLIB =-L$MKLPATH -I$MKLINCLUDE -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
 BLASLIB = 
 
TMGLIB       	= libtmglib.a
LIBS		= $(SUPERLULIB) $(BLASLIB)  ../../metis-4.0/libmetis.a

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = ar
ARCHFLAGS    = cr
RANLIB       = ranlib

CC           = icc
CFLAGS       = -DPRNTlevel=0 -O2 -m64 # -O3 #-g
FORTRAN	     = ifort #f90 #g77
FFLAGS       = -O2 -m64
LOADER       = $(CC) #$(FORTRAN)
LOADOPTS     =

#
#  C preprocessor defs for compilation for the Fortran interface
#  (-DNoChange, -DAdd_, -DAdd__, or -DUpCase)
#
CDEFS        = -DAdd_ #-DNoChange #-DAdd_
#
# The directory in which Matlab is installed
#
MATLAB	     = /usr/sww/matlab

