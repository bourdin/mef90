MEF90_DIM=2
all: SimplePoissonNG${MEF90_DIM}D

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables
include $(MEF90_DIR)/MEF90/Makefile.include

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}
	
SimplePoissonNG${MEF90_DIM}D: ${PETSC_ARCH} SimplePoissonNG.o   ${PETSC_ARCH}/${MEF90_DIM}SimplePoissonNG.o \
                             m_Poisson.o ${PETSC_ARCH}/${MEF90_DIM}m_Poisson.o 
	${FLINKER} -o ${PETSC_ARCH}/SimplePoissonNG${MEF90_DIM}D ${PETSC_ARCH}/${MEF90_DIM}SimplePoissonNG.o  ${PETSC_ARCH}/${MEF90_DIM}m_Poisson.o  ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

SimplePoissonNG.o: ${PETSC_ARCH} m_Poisson.o ${PETSC_ARCH}/${MEF90_DIM}m_Poisson.o

m_Poisson.o: ${PETSC_ARCH}

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}


