PROJ      = TestSieve
MEF90DIR  = ../MEF90
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 

all: ${PROJ} TestLocal  TestAssembly

include ${PETSC_DIR}/conf/base

OBJS = ${PROJ}.o m_${PROJ}.o

FFLAGS = -I. -I${MEF90DIR} -I${M_RUPT_STRUCT_DIR} ${EXO_CFLAGS} -ftrace=full -ffree-line-length-huge

${PROJ}: ${PROJ}.o
	-${FLINKER} -o ${PROJ} ${OBJS} ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
	
TestNodeSet: TestNodeSet.o
	-${FLINKER} -o TestNodeSet TestNodeSet.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestLocal: TestLocal.o
	-${FLINKER} -o TestLocal TestLocal.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	
TestAssembly: TestAssembly.o
	-${FLINKER} -o TestAssembly TestAssembly.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestSectionInt: TestSectionInt.o
	-${FLINKER} -o TestSectionInt TestSectionInt.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

# Dependences
{PROJ}: ${PROJ}.o m_${PROJ}.o
${PROJ}.o: m_${PROJ}.o
