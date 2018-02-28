export PETSC_DIR=/opt/HPC/petsc-3.3-p7
export PETSC_ARCH=gcc-mef90-O
export MEF90_DIR=/opt/HPC/mef90
export SNLP_DIR=/opt/HPC/snlp-gcc-mef90-O
export MPI_DIR=/usr/lib64/mpich-3.2
export PATH=${MPI_DIR}/bin:${PETSC_DIR}/bin:${PETSC_DIR}/${PETSC_ARCH}/bin:${MEF90_DIR}/bin:${MEF90_DIR}/bin/${PETSC_ARCH}:${PATH}
export PYTHONPATH=${MEF90_DIR}/python:${PYTHONPATH}

