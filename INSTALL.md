# Installing mef90 and vDef
vDef is build on top of mef90, a parallel unstructured finite element library and relies heavily on a patched version of petsc-3.3.
Building vDef requires:
  * C, C++, and fortran compilers
  * MPI with support for C++ and fortran90 (standard)
  * mercurial, python are required in order to build PETSc
  
Additionally, some utilities requires the following python modules:
  * numpy
  * matplotlib
  * JSON and argparse

In all that follows, it is assumed that the environment variable MEF90_DIR points to the root of the mef90 installation:
```bash
    [bourdin@head ~]$ echo $MEF90_DIR
    /home/bourdin/Development/mef90/mef90-sieve
    [bourdin@head ~]$ ls $MEF90_DIR
    bin       Makefile.include  mef90version.h  patches  tags        ThermoElasticity
    HeatXfer  m_DefMech         m_HeatXfer      python   TestMeshes  vDef
    Makefile  MEF90             objs            sanson   Tests
```
The actual content of the $MEF90_DIR folder may be somewhat different from the one shown here.

## Building PETSc-3.3:
  * More instructions are provided at http://www.mcs.anl.gov/petsc/petsc-as/documentation/installation.html
  * vDef is *NOT* compatible with petsc-3.4 and later. Port is under way.
  * Clone a heavily patched version of petsc-3.3:
    ```bash
       [bourdin@head petsc-3.3-p7]$ git clone --single-branch --branch mef90-3.3 https://gitlab.com/blaisebourdin/petsc.git petsc-3.3-mef90
    ```
  * Set the PETSC_DIR environment variable to point to the extracted folder
    ```bash
       [bourdin@head petsc-3.3-mef90]$ echo $PETSC_DIR
       /share/apps/petsc-3.3-mef90
    ```
   * set the PETSC_ARCH environment variable to a meaningful value. this value will be used by mef90 in order to allow out of tree build. 
     *IMPORTANT*: when using the gcc compiler family, PETSC_ARCH *MUST* contain the string "gcc". mef90 tests for this string in order to set the proper compiler options for free-form fortran source code.
     ```bash
       [bourdin@head petsc-3.3-mef90]$ echo $PETSC_ARCH
       Linux-gcc4.4-mef90-O
     ```

    
### Configure petsc-3.3. 
  mef90 requires the following external packages: `netcdf` and `exodusii`, `metis` and `parmetis`, `chaco`, `boost`, `triangle`. 
  multi-threaded build require a recent version of cmake.
  the ml and hypre preconditioners are not mandatory but can drastically improve solver performances in some problems.
  petsc must be configured with fortran datatypes, sieve, and C++ as the C language.

  As part of its setup, petsc will download and compile dependencies. On a system without internet access, one can get a list of all packages that need download then compile petsc. This is a 2 steps process:
  1. Run the configure script with `--with-packages-download-dir=<directory>` option. This will return a list of packages required and their location (for some reason, the location of the requested MPI will not show up)
  2. Download the packages from a machine with internet access and place them in the directory specified in step 1 of the build system.
  3. Re-run the command from step 1.

  For instance, 
  ```bash
    [MacBookGray:petsc-3.3-mef90]$ ./configure --download-chaco=1 --download-exodusii=1 --download-hypre=1 --download-ml=1 --download-netcdf=1 --download-sowing=1 --download-yaml=1 --with-packages-download-dir=tarballs/ 
     ===============================================================================
                  Configuring PETSc to compile on your system                       
     ===============================================================================
     Download the following packages to /opt/HPC/petsc-3.3-mef90/tarballs 
     
     netcdf ['https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-c-4.7.4.tar.gz']
     chaco ['http://ftp.mcs.anl.gov/pub/petsc/externalpackages/Chaco-2.2.tar.gz']
     exodusii ['https://github.com/gsjaardema/seacas/archive/0ac6baa46262ab17a7b3c0a07b8832ce8d088324.tar.gz']
     yaml ['http://pyyaml.org/download/libyaml/yaml-0.1.4.tar.gz']
     hypre ['http://ftp.mcs.anl.gov/pub/petsc/externalpackages/hypre-2.8.0b.tar.gz']
     sowing ['http://ftp.mcs.anl.gov/pub/petsc/externalpackages/sowing-1.1.26-p1.tar.gz']
     ml ['http://ftp.mcs.anl.gov/pub/petsc/externalpackages/ml-6.2-p3.tar.gz']
  ``` 
  Then run the script again


#### Linux system

   * On a RHEL linux system with the GNU compiler suite, the following configuration is a good starting point for a build with optimization
     ```bash
       ./configure                         \
         COPTFLAGS='-O3 -march=native'     \
         CXXOPTFLAGS='-O3 -march=native'   \
         FOPTFLAGS='-O3 -march=native'     \
         --download-boost=1                \
         --download-chaco=1                \
         --download-exodusii=1             \
         --download-metis=1                \
         --download-netcdf=1               \
         --download-parmetis=1             \
         --download-sowing=1               \
         --download-triangle=1             \
         --download-yaml=1                 \
         --with-clanguage=C++              \
         --with-cmake=cmake                \
         --with-debugging=0                \
         --with-fortran-datatypes          \
         --with-mpi-dir=$MPI_HOME          \
         --with-shared-libraries=1         \
         --with-sieve
     ```
  
   * With the intel compilers, and MKL, I use
     ```bash
      ./configure                           \
          --COPTFLAGS=-O3                   \
          --CXXOPTFLAGS=-O3                 \
          --FOPTFLAGS=-O3                   \
          --LDFLAGS=-lstdc++                \
          --download-boost=1                \
          --download-chaco=1                \
          --download-exodusii=1             \
          --download-metis=1                \
          --download-netcdf=1               \
          --download-parmetis=1             \
          --download-sowing=1               \
          --download-sowing-cc=/opt/rh/devtoolset-8/root/usr/bin/gcc     \
          --download-sowing-cxx=/opt/rh/devtoolset-8/root/usr/bin/g++    \
          --download-sowing-cpp=/opt/rh/devtoolset-8/root/usr/bin/cpp    \
          --download-sowing-cxxcpp=/opt/rh/devtoolset-8/root/usr/bin/cpp \
          --download-triangle=1             \
          --download-yaml=1                 \
          --with-blas-lapack-dir=$MKLROOT   \
          --with-clanguage=C++              \
          --with-cmake=cmake                \
          --with-debugging=0                \
          --with-fortran-datatypes          \
          --with-mpi-dir=$MPI_HOME          \
          --with-pic                        \
          --with-shared-libraries=1         \
          --with-sieve                      \
          --with-vendor-compilers=intel
     ```
The version of sowing included in this patched version of petsc requires a recent version of gcc and will not compile with the intel compiler, which is why one needs to explicitly specify the path of a recent gcc.

#### On a macOS system with gcc 10.2 (tested on intel-based macOS catalina and big sur) 
The version of boost included with petsc-3.3 does not compile with the system clang compilers, hence using a version of MPI compiled with the system C++ compiler (such as the one from homebrew) is not possible. Instead, it is recommended to let pets compile its own version of MPICH using gcc.

The `autoconf automake cmake gcc make` packages need to be installed (with homebrew for instance):
  
Once these are installed, configure petsc (with debugging) with
```bash
   ./configure                                   \
       --CC=/usr/local/bin/gcc-10                \
       --CFLAGS='-fgnu89-inline -Wno-deprecated' \
       --CXX=/usr/local/bin/g++-10               \
       --CXXFLAGS='-Wno-deprecated'              \
       --FC=/usr/local/bin/gfortran-10           \
       --FFLAGS='-fallow-argument-mismatch'      \
       --download-chaco=1                        \
       --download-exodusii=1                     \
       --download-hypre=1                        \
       --download-ml=1                           \
       --download-mpich=1                        \
       --download-mpich-shared=1                 \
       --download-netcdf=1                       \
       --download-sowing=1                       \
       --download-sowing-cc='gcc-10 -std=c89'    \
       --download-yaml=1                         \
       --with-boost=1                            \
       --with-boost-dir=/usr/local               \
       --with-c2html=0                           \
       --with-clanguage=C++                      \
       --with-clib-autodetect=0                  \
       --with-cmake=cmake                        \
       --with-debugging=1                        \
       --with-fortranlib-autodetect=0            \
       --with-fortran-datatypes=1                \
       --with-shared-libraries=1                 \
       --with-sieve                              \
       --with-x11=1                              
```
or substitute your favorite compiler optimizations and disable debugging for an optimized build

In case of problems with X11, try `--with-x=0`


### Build petsc-3.3
  Follow the on-screen instruction to compile petsc from there (`make PETSC_DIR=.... PETSC_ARCH=...`)

### Optional: exodus python bindings
  Many of the mesh translator or python scripts included in mef90 require the python binding for exodus are installed in `$PETSC_DIR/$PETSC_ARCH/lib`
  It is therefore recommended to add this folder to the `PYTHONPATH`:
  ```bash
     [bourdin@head ~]$ export PYTHONPATH=$PYTHONPATH:$PETSC_DIR/$PETSC_ARCH/lib
  ```

### Recommended: snlp (required for plasticity and soon for some models of unilateral contact)

Set $SNLP_DIR to the location where snlp will be installed. Remark that SNLP relies on PETSc for its makefile system, so using mumtiple builds of PETSc will require using multiple builds of SNLP. THen clone, build, and install snlp
 ```bash
 [bourdin@head ~]$ git clone https://github.com/bourdin/snlp.git
 [bourdin@head ~]$ make
 [bourdin@head ~]$ make install
 ```

## Building vDef
From there, it should be as simple as 
   ```bash
      [bourdin@head ~]$ cd $MEF90_DIR; make
   ```
Note that the default setting is to link with shared libraries, and set their path using rpath (so that `$LD_LIBRARY_PATH` or equivalent does not have to be set).
This means that `$PETSC_DIR/$PETSC_ARCH/lib` needs to be readable from the compute nodes. If PETSc libraries are moved, use chrpath to change the search path after building vDef

For some reason, on macOS only, the rpath for the exodus library is not properly set, so that the environment variable `DYLD_FALLBACK_LIBRARY_PATH` must be set to `$PETSC_DIR/$PETSC_ARCH/lib`, possibly in the shell configuration files.

## Testing:
  run `make test` in `$MEF90_DIR/HeatXfer`, `$MEF90_DIR/ThermoElasticity`, and `$MEF90_DIR/vDef`
  Differences in number of iterations, or round-off error are acceptable
  
  Note that make test will try to run mpi jobs directly. It may be necessary to run make tests in an interactive MPI job session.
  The MPI job launcher can be changed by setting the MPIEXEC environment variable, and the number of processors by setting NP

