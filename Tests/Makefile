all: TestBoundaryElement3D TestBoundaryElement2D

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables
include ${MEF90_DIR}/MEF90/Makefile.include



${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}
	
TestBoundaryElement3D:      ${PETSC_ARCH} TestBoundaryElement3D.o ${PETSC_ARCH}/TestBoundaryElement3D.o
	-${FLINKER} -o ${PETSC_ARCH}/TestBoundaryElement3D ${PETSC_ARCH}/TestBoundaryElement3D.o ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestBoundaryElement2D:      ${PETSC_ARCH} TestBoundaryElement2D.o ${PETSC_ARCH}/TestBoundaryElement2D.o
	-${FLINKER} -o ${PETSC_ARCH}/TestBoundaryElement2D ${PETSC_ARCH}/TestBoundaryElement2D.o ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestMaterials:  ${PETSC_ARCH} TestMaterials.o ${PETSC_ARCH}/TestMaterials.o
	-${FLINKER} -o ${PETSC_ARCH}/TestMaterials ${PETSC_ARCH}/TestMaterials.o ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}



