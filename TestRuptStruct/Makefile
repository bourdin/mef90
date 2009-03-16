PROJ      = TestRuptStruct
MEF90DIR  = ../MEF90
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 
M_RUPTSTRUCTDIR  = ../m_RuptStruct
M_RUPTSTRUCTOBJS = $(M_RUPTSTRUCTDIR)/m_RuptStruct.o

all: ${PROJ}

include ${PETSC_DIR}/conf/base

FFLAGS = -DPB_${DIM} -I. -I${MEF90DIR} -I${M_RUPTSTRUCTDIR}

${PROJ}: ${PROJ}.o
	${FLINKER} -o ${PROJ} ${PROJ}.o ${MEF90OBJS} ${M_RUPTSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
	
