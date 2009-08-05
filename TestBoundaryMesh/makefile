ALL: TestExodus TestMeshExodus
include ${PETSC_DIR}/conf/base

TestExodus: TestExodus.o  chkopts
	-${CLINKER} -o TestExodus TestExodus.o ${PETSC_LIB}

TestMeshExodus: TestMeshExodus.o  chkopts
	-${CLINKER} -o TestMeshExodus TestMeshExodus.o ${PETSC_LIB}
