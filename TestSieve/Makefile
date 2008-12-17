PROJ      = TestSieve
MEF90DIR  = ../MEF90
MEF90OBJS = $(MEF90DIR)/m_Constantes.o    $(MEF90DIR)/m_AlgebLin.o    \
            $(MEF90DIR)/m_MEF_Types.o     $(MEF90DIR)/m_MEF_InitElem.o\
            $(MEF90DIR)/m_MEF90.o         $(MEF90DIR)/m_Utils.o       \
            $(MEF90DIR)/m_MEF_MPI.o
            
#EXO_CFLAGS = -I$(HPCDIST_DIR)/include
#EXO_LFLAGS = -L$(HPCDIST_DIR)/cbind/src -L$(HPCDIST_DIR)/forbind/src -lexoIIv2for -lexoIIv2c -lnetcdf -lmetis -lparmetis

all: ${PROJ}

include ${PETSC_DIR}/conf/base

OBJS = ${PROJ}.o m_${PROJ}.o

FFLAGS = -I. -I${MEF90DIR} -I${M_RUPT_STRUCT_DIR} ${EXO_CFLAGS}

${PROJ}: ${PROJ}.o
	-${FLINKER} -o ${PROJ} ${OBJS} ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${EXO_LFLAGS}

# Dependences
{PROJ}: ${PROJ}.o m_${PROJ}.o
${PROJ}.o: m_${PROJ}.o
