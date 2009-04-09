PROJ      = TestRuptStruct
MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 
M_RUPTSTRUCTDIR  = ../m_RuptStruct/${PETSC_ARCH}
M_RUPTSTRUCTOBJS = ${M_RUPTSTRUCTDIR}/m_RuptStruct.o

all: ${PROJ}

include ${PETSC_DIR}/conf/base

FFLAGS = -I. -I${MEF90DIR} -I${M_RUPTSTRUCTDIR}
.F90.o: 
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${FFLAGS} ${FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

#${PROJ}.o: ${PETSC_ARCH}

${PROJ}: ${PETSC_ARCH} ${PROJ}.o ${PETSC_ARCH}/${PROJ}.o
	${FLINKER} -o ${PETSC_ARCH}/${PROJ} ${PETSC_ARCH}/${PROJ}.o ${MEF90OBJS} ${M_RUPTSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
	
