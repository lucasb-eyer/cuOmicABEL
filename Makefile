include ./make.inc

SRCDIR = ./src
DRIVER = ./HP-GWAS

CFLAGS+=-g -Wall -I $(SRCDIR)/  # -D__WORDSIZE=64
LDLIBS += -lm

# Use these if you want GPU support built in too. Requires a valid CUDA install.
CFLAGS += -DWITH_GPU -I$(CUDA_ROOT)/include
LDLIBS += -L$(CUDA_ROOT)/lib64 -lcublas

SRCS = $(SRCDIR)/CLAK_GWAS.c $(SRCDIR)/fgls_chol.c $(SRCDIR)/fgls_chol_gpu.c $(SRCDIR)/fgls_eigen.c $(SRCDIR)/wrappers.c $(SRCDIR)/timing.c $(SRCDIR)/statistics.c $(SRCDIR)/REML.c $(SRCDIR)/optimization.c $(SRCDIR)/ooc_BLAS.c $(SRCDIR)/double_buffering.c $(SRCDIR)/utils.c $(SRCDIR)/GWAS.c $(SRCDIR)/databel.c
OBJS = $(SRCS:.c=.o)

.PHONY: all clean

all: $(DRIVER)

$(DRIVER): $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) $(LDLIBS) -o $@

clean:
	$(RM) $(OBJS)
	$(RM) $(DRIVER) $(EXE_GPU);
	$(RM) $(SRCDIR)/*mod*
	$(RM) $(SRCDIR)/*opari_GPU*


src/CLAK_GWAS.o: src/CLAK_GWAS.c src/wrappers.h src/utils.h src/GWAS.h \
 src/databel.h src/timing.h src/REML.h src/fgls_chol.h src/fgls_eigen.h \
 src/double_buffering.h
src/GWAS.o: src/GWAS.c src/utils.h src/GWAS.h src/databel.h src/wrappers.h
src/REML.o: src/REML.c src/options.h src/blas.h src/lapack.h src/ooc_BLAS.h \
 src/wrappers.h src/utils.h src/GWAS.h src/databel.h src/statistics.h \
 src/optimization.h src/REML.h
src/databel.o: src/databel.c src/databel.h src/wrappers.h
src/double_buffering.o: src/double_buffering.c src/GWAS.h src/databel.h \
 src/wrappers.h src/double_buffering.h
src/fgls_chol.o: src/fgls_chol.c src/blas.h src/lapack.h src/options.h \
 src/GWAS.h src/databel.h src/wrappers.h src/timing.h \
 src/double_buffering.h src/utils.h src/fgls_chol.h
src/fgls_chol_gpu.o: src/fgls_chol_gpu.c src/blas.h src/lapack.h src/options.h \
 src/GWAS.h src/databel.h src/wrappers.h src/timing.h \
 src/double_buffering.h src/utils.h src/fgls_chol_gpu.h
src/fgls_eigen.o: src/fgls_eigen.c src/blas.h src/lapack.h src/options.h \
 src/GWAS.h src/databel.h src/wrappers.h src/timing.h \
 src/double_buffering.h src/ooc_BLAS.h src/utils.h src/fgls_eigen.h
src/ooc_BLAS.o: src/ooc_BLAS.c src/blas.h src/lapack.h src/options.h \
 src/GWAS.h src/databel.h src/wrappers.h src/utils.h \
 src/double_buffering.h src/ooc_BLAS.h
src/optimization.o: src/optimization.c src/optimization.h
src/statistics.o: src/statistics.c src/statistics.h
src/timing.o: src/timing.c src/timing.h
src/utils.o: src/utils.c \
 src/GWAS.h src/databel.h src/wrappers.h src/utils.h
src/wrappers.o: src/wrappers.c src/GWAS.h src/databel.h src/wrappers.h
