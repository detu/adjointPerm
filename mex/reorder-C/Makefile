MATLABVERSION = 76
MEX           = mex -largeArrayDims -g


PEDANTIC=-Wall -pedantic -Wextra -Wformat-nonliteral             \
	-Wcast-align -Wpointer-arith -Wbad-function-cast         \
	-Wmissing-prototypes -Wstrict-prototypes                 \
	-Wmissing-declarations -Winline -Wundef -Wnested-externs \
	-Wcast-qual -Wshadow -Wconversion -Wwrite-strings        \
	-Wno-conversion -Wchar-subscripts -Wredundant-decls


GCC        = gcc -fPIC #-std=c99
LIBINCGCC  = -L/work/lib -L/usr/local/src/SuiteSparse-3.2/lib \
	     -L/home/bska/lib/lib/suitesparse/3.1.0/lib

ICC        = /opt/intel/cce/10.1.015/bin/icc -fPIC
LIBINCICC  = -L/opt/intel/cce/10.1.015/lib/


CC     = $(GCC)
LIBINC = $(LIBINCGCC)

INCJUMBO    = -I/work/opt/matlab76/extern/include/ \
	      -I/work/include
INC    = $(INCJUMBO)                                   \
	-I/usr/local/matlab76/extern/include/          \
	-I/usr/local/include                           \
	-I/usr/local/src/SuiteSparse-3.2/include       \
	-I/home/bska/lib/lib/suitesparse/3.1.0/include \
	-I/usr/local/share/matlab/R2008a/extern/include

#INC    = -I/usr/local/matlab14/extern/include


LIB    = $(LIBINC)  -lm -lblas -llapack -lblas  -lumfpack -lamd
CFLAGS = -g $(INC) -Wall  -DMEXFILE=1 -DMATLABVERSION=$(MATLABVERSION) #$(PEDANTIC)
SHELL  = /bin/bash

all:    tags depend implicitupwind

# Write date to today.h
today.h:;	@(echo char today[512]=\"Compiled `date +%a\ %b\ %d\ %Y\ %r`\"\;>today.h)

# Force recompilation each time to include today.h
implicitupwind.o:	today.h
	$(CC) $(CFLAGS)  -c -o implicitupwind.o implicitupwind.c

OBJ3 =  implicitupwind.o  utils.o mexutils.o system.o  strongcomp.o \
	sparse.o nlsolvers.o lu.o

implicitupwind:	today.h $(OBJ3)
	$(MEX) $(OBJ3) -o $@ $(LIBINC) $(LIB) >/dev/null

install:
	#cp twophaseUpwReorder.m ../../transport
	cd ../../transport/ &&	ln -s ../mex/reorder-C/implicitupwind.mexa64
.PHONY: clean depend all today.h

clean:; rm -f *~ $(OBJ3) today.h implicitupwind.mex*  \
	TAGS; makedepend

tags : $(OBJ:.o=.c)  $(wildcard *.h)
	etags *.c *.h

depend :; makedepend $(INC) -f Makefile *.c 2>/dev/null

# DO NOT DELETE THIS LINE -- make depend depends on it.

implicitupwind.o: /usr/include/string.h /usr/include/features.h
implicitupwind.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
implicitupwind.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
implicitupwind.o: /usr/include/xlocale.h
implicitupwind.o: /usr/local/share/matlab/R2008a/extern/include/mex.h
implicitupwind.o: /usr/local/share/matlab/R2008a/extern/include/matrix.h
implicitupwind.o: /usr/local/share/matlab/R2008a/extern/include/tmwtypes.h
implicitupwind.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
implicitupwind.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
implicitupwind.o: /usr/include/bits/posix2_lim.h
implicitupwind.o: /usr/include/bits/xopen_lim.h /usr/include/bits/stdio_lim.h
implicitupwind.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
implicitupwind.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
implicitupwind.o: /usr/include/bits/endian.h /usr/include/sys/types.h
implicitupwind.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
implicitupwind.o: /usr/include/time.h /usr/include/sys/select.h
implicitupwind.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
implicitupwind.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
implicitupwind.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
implicitupwind.o: /usr/include/assert.h /usr/include/stdio.h
implicitupwind.o: /usr/include/libio.h /usr/include/_G_config.h
implicitupwind.o: /usr/include/wchar.h /usr/include/bits/wchar.h
implicitupwind.o: /usr/include/gconv.h /usr/include/bits/sys_errlist.h
implicitupwind.o: /usr/include/math.h /usr/include/bits/huge_val.h
implicitupwind.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
implicitupwind.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
implicitupwind.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
implicitupwind.o: sparse.h strongcomp.h system.h utils.h mexutils.h
implicitupwind.o: /usr/include/signal.h /usr/include/bits/signum.h
implicitupwind.o: /usr/include/bits/siginfo.h /usr/include/bits/sigaction.h
implicitupwind.o: /usr/include/bits/sigcontext.h /usr/include/bits/sigstack.h
implicitupwind.o: /usr/include/sys/ucontext.h /usr/include/bits/sigthread.h
lu.o: /usr/include/stdlib.h /usr/include/features.h /usr/include/sys/cdefs.h
lu.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
lu.o: /usr/include/gnu/stubs-64.h /usr/include/bits/waitflags.h
lu.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
lu.o: /usr/include/bits/endian.h /usr/include/xlocale.h
lu.o: /usr/include/sys/types.h /usr/include/bits/types.h
lu.o: /usr/include/bits/typesizes.h /usr/include/time.h
lu.o: /usr/include/sys/select.h /usr/include/bits/select.h
lu.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
lu.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
lu.o: /usr/include/alloca.h /usr/include/stdio.h /usr/include/libio.h
lu.o: /usr/include/_G_config.h /usr/include/wchar.h /usr/include/bits/wchar.h
lu.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
lu.o: /usr/include/bits/sys_errlist.h lu.h
mexutils.o: /usr/local/share/matlab/R2008a/extern/include/mex.h
mexutils.o: /usr/local/share/matlab/R2008a/extern/include/matrix.h
mexutils.o: /usr/local/share/matlab/R2008a/extern/include/tmwtypes.h
mexutils.o: /usr/include/limits.h /usr/include/features.h
mexutils.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
mexutils.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
mexutils.o: /usr/include/bits/posix1_lim.h /usr/include/bits/local_lim.h
mexutils.o: /usr/include/linux/limits.h /usr/include/bits/posix2_lim.h
mexutils.o: /usr/include/bits/xopen_lim.h /usr/include/bits/stdio_lim.h
mexutils.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
mexutils.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
mexutils.o: /usr/include/bits/endian.h /usr/include/xlocale.h
mexutils.o: /usr/include/sys/types.h /usr/include/bits/types.h
mexutils.o: /usr/include/bits/typesizes.h /usr/include/time.h
mexutils.o: /usr/include/sys/select.h /usr/include/bits/select.h
mexutils.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
mexutils.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
mexutils.o: /usr/include/alloca.h /usr/include/assert.h /usr/include/stdio.h
mexutils.o: /usr/include/libio.h /usr/include/_G_config.h
mexutils.o: /usr/include/wchar.h /usr/include/bits/wchar.h
mexutils.o: /usr/include/gconv.h /usr/include/bits/sys_errlist.h
mexutils.o: /usr/include/string.h sparse.h mexutils.h /usr/include/signal.h
mexutils.o: /usr/include/bits/signum.h /usr/include/bits/siginfo.h
mexutils.o: /usr/include/bits/sigaction.h /usr/include/bits/sigcontext.h
mexutils.o: /usr/include/bits/sigstack.h /usr/include/sys/ucontext.h
mexutils.o: /usr/include/bits/sigthread.h
nlsolvers.o: /usr/include/stdlib.h /usr/include/features.h
nlsolvers.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
nlsolvers.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
nlsolvers.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
nlsolvers.o: /usr/include/endian.h /usr/include/bits/endian.h
nlsolvers.o: /usr/include/xlocale.h /usr/include/sys/types.h
nlsolvers.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
nlsolvers.o: /usr/include/time.h /usr/include/sys/select.h
nlsolvers.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
nlsolvers.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
nlsolvers.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
nlsolvers.o: /usr/include/stdio.h /usr/include/libio.h
nlsolvers.o: /usr/include/_G_config.h /usr/include/wchar.h
nlsolvers.o: /usr/include/bits/wchar.h /usr/include/gconv.h
nlsolvers.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
nlsolvers.o: /usr/include/string.h /usr/include/assert.h /usr/include/math.h
nlsolvers.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
nlsolvers.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
nlsolvers.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
nlsolvers.o: /usr/include/bits/mathcalls.h nlsolvers.h
sparse.o: /usr/include/stdlib.h /usr/include/features.h
sparse.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
sparse.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
sparse.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
sparse.o: /usr/include/endian.h /usr/include/bits/endian.h
sparse.o: /usr/include/xlocale.h /usr/include/sys/types.h
sparse.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
sparse.o: /usr/include/time.h /usr/include/sys/select.h
sparse.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
sparse.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
sparse.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
sparse.o: /usr/include/stdio.h /usr/include/libio.h /usr/include/_G_config.h
sparse.o: /usr/include/wchar.h /usr/include/bits/wchar.h /usr/include/gconv.h
sparse.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
sparse.o: /usr/include/string.h /usr/include/assert.h /usr/include/math.h
sparse.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
sparse.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
sparse.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
sparse.o: /usr/include/bits/mathcalls.h
sparse.o: /usr/local/include/SuiteSparse/umfpack.h
sparse.o: /usr/local/include/SuiteSparse/UFconfig.h /usr/include/limits.h
sparse.o: /usr/include/bits/posix1_lim.h /usr/include/bits/local_lim.h
sparse.o: /usr/include/linux/limits.h /usr/include/bits/posix2_lim.h
sparse.o: /usr/include/bits/xopen_lim.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_symbolic.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_numeric.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_solve.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_free_symbolic.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_free_numeric.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_defaults.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_qsymbolic.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_wsolve.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_triplet_to_col.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_col_to_triplet.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_transpose.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_scale.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_get_lunz.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_get_numeric.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_get_symbolic.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_save_numeric.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_load_numeric.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_save_symbolic.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_load_symbolic.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_get_determinant.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_report_status.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_report_info.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_report_control.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_report_matrix.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_report_triplet.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_report_vector.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_report_symbolic.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_report_numeric.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_report_perm.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_timer.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_tictoc.h
sparse.o: /usr/local/include/SuiteSparse/amd.h
sparse.o: /usr/local/include/SuiteSparse/umfpack_global.h lu.h sparse.h
strongcomp.o: /usr/include/stdlib.h /usr/include/features.h
strongcomp.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
strongcomp.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
strongcomp.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
strongcomp.o: /usr/include/endian.h /usr/include/bits/endian.h
strongcomp.o: /usr/include/xlocale.h /usr/include/sys/types.h
strongcomp.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
strongcomp.o: /usr/include/time.h /usr/include/sys/select.h
strongcomp.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
strongcomp.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
strongcomp.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
strongcomp.o: /usr/include/stdio.h /usr/include/libio.h
strongcomp.o: /usr/include/_G_config.h /usr/include/wchar.h
strongcomp.o: /usr/include/bits/wchar.h /usr/include/gconv.h
strongcomp.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
strongcomp.o: /usr/include/string.h /usr/include/assert.h sparse.h
strongcomp.o: strongcomp.h system.h
system.o: /usr/include/stdlib.h /usr/include/features.h
system.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
system.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
system.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
system.o: /usr/include/endian.h /usr/include/bits/endian.h
system.o: /usr/include/xlocale.h /usr/include/sys/types.h
system.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
system.o: /usr/include/time.h /usr/include/sys/select.h
system.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
system.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
system.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
system.o: /usr/include/stdio.h /usr/include/libio.h /usr/include/_G_config.h
system.o: /usr/include/wchar.h /usr/include/bits/wchar.h /usr/include/gconv.h
system.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
system.o: /usr/include/assert.h /usr/include/string.h /usr/include/math.h
system.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
system.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
system.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
system.o: /usr/include/bits/mathcalls.h sparse.h nlsolvers.h system.h
utils.o: /usr/include/stdio.h /usr/include/features.h
utils.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
utils.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
utils.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
utils.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
utils.o: /usr/include/bits/wchar.h /usr/include/gconv.h
utils.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
utils.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
utils.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
utils.o: /usr/include/bits/endian.h /usr/include/xlocale.h
utils.o: /usr/include/sys/types.h /usr/include/time.h
utils.o: /usr/include/sys/select.h /usr/include/bits/select.h
utils.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
utils.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
utils.o: /usr/include/alloca.h /usr/include/sys/time.h /usr/include/string.h
utils.o: /usr/include/assert.h
utils.o: /usr/local/share/matlab/R2008a/extern/include/mex.h
utils.o: /usr/local/share/matlab/R2008a/extern/include/matrix.h
utils.o: /usr/local/share/matlab/R2008a/extern/include/tmwtypes.h
utils.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
utils.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
utils.o: /usr/include/bits/posix2_lim.h /usr/include/bits/xopen_lim.h utils.h
