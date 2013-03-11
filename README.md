cuOmicABEL
==========

Adding CUDA multi-GPU support to OmicABEL for computing GWAS quickly.

Installing
==========
You need the following dependencies:

### GotoBLAS2
OmicABEL heavily relies on [GotoBLAS2](http://www.tacc.utexas.edu/tacc-projects/gotoblas2).
Use the following to install it:

```bash
$ cd ~
$ wget http://www.tacc.utexas.edu/documents/13601/b58aeb8c-9d8d-4ec2-b5f1-5a5843b4d47b -O GotoBLAS2-1.13.tar.gz
$ tar -xzvf GotoBLAS2-1.13.tar.gz
$ cd GotoBLAS2
$ make TARGET=NEHALEM
```

You can drop the `TARGET=NEHALEM` option, but you will then get compile errors
on modern Intel-processors because it doesn't recognize them. `make clean` in
case that happens.

### LAPACK
Same here; install it using the following commands:

```bash
$ cd ~
$ wget http://www.netlib.org/lapack/lapack-3.4.2.tgz
$ tar -xzvf lapack-3.4.2.tgz
$ mv lapack-3.4.2 lapack
$ cd lapack
$ cp make.inc.example make.inc
```

Edit that file so that at the bottom, it shows the following:

```bash
BLASLIB      = ~/GotoBLAS2/libgoto2.a
```

Then, compile it by running `make`. At some point, the tests will likely fail
due to a linker error. We could fix this by linking the tests with pthreads,
but we do not really care about them anyways; as long as there is a
`liblapack.a` file, you're fine. Still need to create a differently named link
to it for OmicABEL: `ln -s liblapack.a lapack_LINUX.a`.

### cuOmicABEL

In case you installed GotoBLAS2 and LAPACK in other directories, you'll need to
edit the `make.inc` file accordingly.

Besides that, cuOmicABEL should be compiled simply by running `make`.
