PROJ      = SimplePoisson
MEF90DIR  = ../MEF90
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 

all: ${PROJ}${DIM}

include ${PETSC_DIR}/conf/base

OBJS = ${PROJ}${DIM}.o m_${PROJ}${DIM}.o
FFLAGS = -DPB_${DIM} -I. -I${MEF90DIR} -ftrace=full -ffree-line-length-huge

${PROJ}${DIM}: ${PROJ}${DIM}.o m_${PROJ}${DIM}.o
	${FLINKER} -o ${PROJ}${DIM} ${OBJS} ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
	
${PROJ}${DIM}.o: m_${PROJ}${DIM}.o
m_${PROJ}${DIM}.o: ${MEF90OBJS}

# Dependences
#${PROJ}${DIM}: ${PROJ}${DIM}.o m_${PROJ}${DIM}.o
#${PROJ}${DIM}.o: m_${PROJ}${DIM}.o
