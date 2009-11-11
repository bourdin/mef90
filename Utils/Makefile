MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 

all: Partitioner EXOGetBB

include ${PETSC_DIR}/conf/base


FFLAGS = ${FC_MODULE_FLAG}${PETSC_ARCH} -I. -I${MEF90DIR} --free-line-length-none
#FFLAGS = -fmod=${PETSC_ARCH} -I. -I${MEF90DIR} -ffree-line-length-huge

.F90.o: 
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${FFLAGS} ${FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/$@ $<


${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

Partitioner:   ${PETSC_ARCH} Partitioner.o ${PETSC_ARCH}/Partitioner.o
	-${FLINKER} -o ${PETSC_ARCH}/Partitioner ${PETSC_ARCH}/Partitioner.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

EXOGetBB:   ${PETSC_ARCH} EXOGetBB.o ${PETSC_ARCH}/EXOGetBB.o
	-${FLINKER} -o ${PETSC_ARCH}/EXOGetBB ${PETSC_ARCH}/EXOGetBB.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}
