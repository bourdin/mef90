# Installing mef90 and vDef
vDef is build on top of mef90, a parallel unstructured finite element library and relies heavily on a patched version of petsc-3.3.
Building vDef requires:
  * C and fortran compilers
  * MPI 
  * git, mercurial, python3, cmake are required in order to build PETSc
  
Additionally, some utilities requires the following python modules:
  * numpy
  * matplotlib
  * JSON and argparse

In all that follows, it is assumed that the environment variable MEF90_DIR points to the root of the mef90 installation:
```bash
    [bourdin@bbserv ~]$ echo $MEF90_DIR
    /home/bourdin/Development/mef90
```
The actual content of the $MEF90_DIR folder may be somewhat different from the one shown here.

## Building PETSc:
  * As of 2023-01-19, vDef uses teh development branch of petsc ("main" branch).
  * More instructions are provided at https://petsc.org/release/install
  * Clone petsc:
    ```bash
       [bourdin@bbserv HPC]$ git clone --single-branch --branch main https://gitlab.com/petsc/petsc.git
    ```
  * Set the PETSC_DIR environment variable to point to the extracted folder
    ```bash
       [bourdin@bbserv petsc]$ echo $PETSC_DIR
       /opt/HPC/petsc
    ```
   * Set the PETSC_ARCH environment variable to a meaningful value. This value will be used by mef90 in order to allow out-of-tree build. THis way, several installation of PETSc and vDef (with various optimisation or debugging informations, for instance) can co-exist on a system.

    
### Configure petsc. 
  vDef requires the following external packages: `exodusii`, `netcdf, `pnetcdf`, `hdf5`, `metis`, `parmetis`, and `zlib`.
  Fortran bindings for the `exodusII` libraries are needed (`--with-exodusii-fortran-bindings` option)
  The ml and hypre preconditioners are not mandatory but can drastically improve solver performances in some problems.

  As part of its setup, petsc will download and compile dependencies. On a system without internet access, one can get a list of all packages that need download then compile petsc. This is a 2 steps process:
  1. Run the configure script with `--with-packages-download-dir=<directory>` option. This will return a list of packages required and their location.
  2. Download the packages from a machine with internet access and place them in the directory specified in step 1 of the build system.
  3. Re-run the command from step 1.

  For instance, 
  ```bash
SiMini:petsc-main (main)$ ./configure --download-exodusii=1 --download-hdf5=1 --download-metis=1 --download-netcdf=1 --download-parmetis=1 --download-pnetcdf=1 -download-zlib --with-packages-download-dir=downloads
Download the following packages to /opt/HPC/petsc/downloads 

zlib ['http://www.zlib.net/zlib-1.2.11.tar.gz', 'http://ftp.mcs.anl.gov/pub/petsc/externalpackages/zlib-1.2.11.tar.gz']
hdf5 ['https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.2/src/hdf5-1.12.2.tar.bz2', 'http://ftp.mcs.anl.gov/pub/petsc/externalpackages/hdf5-1.12.2.tar.bz2']
netcdf ['https://github.com/Unidata/netcdf-c/archive/v4.9.0.tar.gz', 'http://ftp.mcs.anl.gov/pub/petsc/externalpackages/netcdf-4.9.0.tar.gz']
pnetcdf ['git://https://github.com/parallel-netcdf/pnetcdf', 'https://parallel-netcdf.github.io/Release/pnetcdf-1.12.3.tar.gz', 'http://ftp.mcs.anl.gov/pub/petsc/externalpackages/pnetcdf-1.12.3.tar.gz']
hypre ['git://https://github.com/hypre-space/hypre', 'https://github.com/hypre-space/hypre/archive/aff83e81fce70563edc614edfef76fa998ba96fd.tar.gz']
metis ['git://https://bitbucket.org/petsc/pkg-metis.git', 'https://bitbucket.org/petsc/pkg-metis/get/v5.1.0-p10.tar.gz']
parmetis ['git://https://bitbucket.org/petsc/pkg-parmetis.git', 'https://bitbucket.org/petsc/pkg-parmetis/get/v4.0.3-p8.tar.gz']
mumps ['https://graal.ens-lyon.fr/MUMPS/MUMPS_5.5.1.tar.gz', 'http://ftp.mcs.anl.gov/pub/petsc/externalpackages/MUMPS_5.5.1.tar.gz']
scalapack ['git://https://github.com/Reference-ScaLAPACK/scalapack', 'https://github.com/Reference-ScaLAPACK/scalapack/archive/5bad7487f496c811192334640ce4d3fc5f88144b.tar.gz']
triangle ['git://https://bitbucket.org/petsc/pkg-triangle', 'https://bitbucket.org/petsc/pkg-triangle/get/v1.3-p2.tar.gz']
exodusii ['git://https://github.com/gsjaardema/seacas.git', 'https://github.com/gsjaardema/seacas/archive/v2022-08-01.tar.gz']
ctetgen ['git://https://bitbucket.org/petsc/ctetgen', 'https://bitbucket.org/petsc/ctetgen/get/ctetgen-0.10.tar.gz']
```
  Then run the script again


#### Linux system
   * On a RHEL linux system with the GNU compiler suite, the following configuration is a good starting point for a build with optimization
     ```bash
       ./configure                         \
    --COPTFLAGS='-O3 -mcpu=native'   \
    --CXXOPTFLAGS='-O3 -mcpu=native' \
    --FOPTFLAGS='-O3 -mcpu=native'   \
    --CFLAGS='-Wimplicit-function-declaration -Wunused'\
    --FFLAGS='-ffree-line-length-none -fallow-argument-mismatch -Wunused'        \
    --download-exodusii=1             \
    --download-fblaslapack=1          \
    --download-hdf5=1                 \
    --download-metis=1                \
    --download-netcdf=1               \
    --download-parmetis=1             \
    --download-pnetcdf=1              \
    --download-zlib=1                 \
    --with-debugging=0                \
    --with-exodusii-fortran-bindings  \
    --with-mpi-dir=$MPI_HOME          \
    --with-shared-libraries=1         \
    --with-x11=0
     ```

#### On a macOS system with gcc 12.2
A version of MPI (both `mpich` and `open-mpi` installed from `homebrew` are tested) must be installed and their binaries, headers, and libraries must be in the standard search paths (or use the `--with-mpi-dir` optiona):
  
Once these are installed, configure petsc (with debugging) with
```bash
./configure \
    --COPTFLAGS='-O3 -mcpu=native'   \
    --CXXOPTFLAGS='-O3 -mcpu=native' \
    --FOPTFLAGS='-O3 -mcpu=native'   \
    --CFLAGS='-Wimplicit-function-declaration -Wunused'\
    --FFLAGS='-ffree-line-length-none -fallow-argument-mismatch -Wunused'        \
    --download-exodusii=1             \
    --download-hdf5=1                 \
    --download-netcdf=1               \
    --download-parmetis=1             \
    --download-pnetcdf=1              \
    --download-zlib=1                 \
    --with-debugging=0                \
    --with-exodusii-fortran-bindings  \
    --with-shared-libraries=1         \
    --with-x11=0
```
or substitute your favorite compiler optimizations and disable debugging for an optimized build

In case of problems with X11, try `--with-x=0`


### Build petsc
  * Follow the on-screen instruction to compile petsc from there (`make PETSC_DIR=.... PETSC_ARCH=...`)
  * (Recommended) Add `$PETSC_DIR/bin` and `$PETSC_DIR/$PETSC_ARCH/bin` to `$PATH` and `$PETSC_DIR/$PETSC_ARCH/lib` to `$PYTHONPATH`:
  ```bash
     [bourdin@bbserv ~]$ export PATH=$PETSC_DIR/bin:$PETSC_DIR/$PETSC_ARCH/bin:$PATH
     [bourdin@bbserv ~]$ export PYTHONPATH=$PYTHONPATH:$PETSC_DIR/$PETSC_ARCH/lib
  ```

### Recommended: snlp (required for plasticity)

Set $SNLP_DIR to the location where snlp will be installed. Remark that SNLP relies on PETSc for its makefile system, so using multiple builds of PETSc will require using multiple builds of SNLP. Then clone, build, and install snlp
 ```bash
 [bourdin@bbserv ~]$ git clone https://github.com/bourdin/snlp.git
 [bourdin@bbserv ~]$ make
 [bourdin@bbserv ~]$ make install
 ```

## Building vDef
From there, it should be as simple as 
   ```bash
      [bourdin@bbserv ~]$ cd $MEF90_DIR; make
   ```
Note that the default setting is to link with shared libraries, and set their path using rpath (so that `$LD_LIBRARY_PATH` or equivalent does not have to be set).
This means that `$PETSC_DIR/$PETSC_ARCH/lib` needs to be readable from the compute nodes. If PETSc libraries are moved, use chrpath to change the search path after building vDef

## Testing:
  run `make test` in `$MEF90_DIR/HeatXfer`, `$MEF90_DIR/ThermoElasticity`, and `$MEF90_DIR/vDef`
  Differences in number of iterations, or round-off error are acceptable
  
  Note that make test will try to run mpi jobs directly. It may be necessary to run make tests in an interactive MPI job session.
  The MPI job launcher can be changed by setting the MPIEXEC environment variable, and the number of processors by setting NP

