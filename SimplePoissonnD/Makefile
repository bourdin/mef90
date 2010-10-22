MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 

all: SimplePoisson${DIM}

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables
FFLAGS = -DPB_${DIM} -module ${PETSC_ARCH} -I. -I${MEF90DIR}


.F90.o: 
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${FFLAGS} ${FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/${DIM}$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}
	
SimplePoisson${DIM}: ${PETSC_ARCH} SimplePoisson.o   ${PETSC_ARCH/${DIM}SimplePoisson.o \
                             m_SimplePoisson.o ${PETSC_ARCH}/${DIM}m_SimplePoisson.o 
	${FLINKER} -o ${PETSC_ARCH}/SimplePoisson${DIM} ${PETSC_ARCH}/${DIM}SimplePoisson.o  ${PETSC_ARCH}/${DIM}m_SimplePoisson.o  ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
	
SimplePoisson.o:   ${PETSC_ARCH} m_SimplePoisson.o ${PETSC_ARCH}/${DIM}m_SimplePoisson.o

m_SimplePoisson.o: ${PETSC_ARCH}

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}

