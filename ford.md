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

Documentation for the mef90 / vDef codebase, generated with
[FORD](https://forddocs.readthedocs.io).

Note on preprocessor templates: several files are compiled multiple times by the
Makefiles with different macros (`MEF90_DIM=2/3`, `MEF90_ELEMENTTYPE=...`). FORD
preprocesses each file once, so these pages document a single representative
instantiation: `MEF90_DIM=3` with `MEF90_ELEMENTTYPE=MEF90Element3DVect`.
References to the other instantiations (e.g. `*_MEF90Element2DScal`) appear
unlinked.

To rebuild these pages:

```
python3 -m venv .venv-ford && .venv-ford/bin/pip install ford
.venv-ford/bin/ford ford.md
```
