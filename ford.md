---
project: mef90 / vDef
summary: Fortran + PETSc reference implementation of the variational approach to fracture
author: Blaise Bourdin
project_website: https://github.com/bourdin/mef90
src_dir: MEF90
         m_DefMech
         m_HeatXfer
         vDef
         HeatXfer
output_dir: doc/html
fpp_extensions: F90
preprocess: true
preprocessor: /opt/homebrew/bin/gcc-16 -traditional-cpp -E -x f95-cpp-input -D__GFORTRAN__
mathjax_config: doc/mathjax-config.js
css: doc/user.css
page_dir: doc/pages
md_extensions: doi_link
docmark: >
predocmark: !!
macro: MEF90_DIM=3
       MEF90_ELEMENTTYPE=MEF90Element3DVect
       MEF90_ELEMENTTYPE_VECT
exclude: **/m_MEF90_MassMatrixInterface.F90
         **/m_MEF90_NormsInterface.F90
         **/vDefBT.F90
         **/vDefHF.F90
         **/vDefP.F90
         **/vDefUpa.F90
include: /opt/HPC/src/petsc/main/include
         /opt/HPC/src/petsc/main/main-tahoe-gcc15.2-arm64-g/include
         /opt/HPC/src/petsc/main/main-tahoe-gcc15.2-arm64-g/externalpackages/git.exodusii/packages/seacas/libraries/exodus_for/include
         MEF90
display: public
         protected
         private
sort: alpha
source: false
graph: false
search: true
---

Documentation for the mef90 / vDef codebase, a Fortran + PETSc reference
implementation of the variational approach to fracture, generated with
[FORD](https://forddocs.readthedocs.io) from the sources at
[github.com/bourdin/mef90](https://github.com/bourdin/mef90).

Start from **Modules** for the structure of the library, **Derived Types** for
the context and options objects, or **Procedures** to look up a single routine.

In the source listings, identifiers naming a documented module, type, or
procedure link to its page, and calls into PETSc link to the corresponding
manual page on petsc.org.

Several files are compiled more than once by the Makefiles with different
macros, so their pages document one representative instantiation
(`MEF90_DIM=3`, `MEF90_ELEMENTTYPE=MEF90Element3DVect`); references to the other
instantiations appear unlinked.

For how these pages are written, built, and published, see
[Maintaining these pages](page/index.html).
