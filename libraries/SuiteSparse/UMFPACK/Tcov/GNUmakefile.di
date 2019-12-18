all: go

include UFconfig/UFconfig.mk

go: run
	- ( cd UMFPACK/Source ; ./ucov.di )
	- ( cd AMD/Source     ; ./acov.di )

run: prog
	- ./ut > ut.out
	- tail ut.out
	#- $(RM) ut.out

prog:
	( cd UMFPACK ; make library )
	( cd AMD ; make library )
	$(CC) -DDINT $(CFLAGS) $(UMFPACK_CONFIG) -IUMFPACK/Source -IUMFPACK/Include -IAMD/Source -IAMD/Include -IUFconfig -o ut ut.c UMFPACK/Source/libumfpack.a AMD/Source/libamd.a CHOLMOD/Lib/libcholmod.a CAMD/Lib/libcamd.a COLAMD/Lib/libcolamd.a metis-4.0/libmetis.a CCOLAMD/Lib/libccolamd.a $(LIB)

utcov:
	- ( cd UMFPACK/Source ; ./ucov.di )
	- ( cd AMD/Source     ; ./acov.di )


purge:
	( cd UMFPACK ; make purge )
	( cd AMD ; make purge )
