PROJ      = TestSieve
MEF90DIR  = ../MEF90
MEF90OBJS = $(MEF90DIR)/m_Constantes.o    $(MEF90DIR)/m_AlgebLin.o    \
            $(MEF90DIR)/m_MEF_Types.o     $(MEF90DIR)/m_MEF_Elements.o\
            $(MEF90DIR)/m_MEF90.o         $(MEF90DIR)/m_Utils.o       \
            $(MEF90DIR)/m_MEF_EXO.o       $(MEF90DIR)/m_MEF_MPI.o     \
            $(MEF90DIR)/m_MEF_Sieve.o 
all: ${PROJ}

include ${PETSC_DIR}/conf/base

OBJS = ${PROJ}.o m_${PROJ}.o

FFLAGS = -I. -I${MEF90DIR} -I${M_RUPT_STRUCT_DIR} ${EXO_CFLAGS} -ftrace=full -ffree-line-length-huge

${PROJ}: ${PROJ}.o
	-${FLINKER} -o ${PROJ} ${OBJS} ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${EXO_LFLAGS}
	
TestNodeSet: TestNodeSet.o
	-${FLINKER} -o TestNodeSet TestNodeSet.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${EXO_LFLAGS}

# Dependences
{PROJ}: ${PROJ}.o m_${PROJ}.o
${PROJ}.o: m_${PROJ}.o
