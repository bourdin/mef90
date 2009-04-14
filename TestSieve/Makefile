MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 

all: TestSieve TestSieve3D TestNodeSet TestLocal TestAssembly TestSectionInt

include ${PETSC_DIR}/conf/base


#FFLAGS = -module ${PETSC_ARCH} -I. -I${MEF90DIR} -ftrace=full -ffree-line-length-huge
FFLAGS = -fmod=${PETSC_ARCH} -I. -I${MEF90DIR} -ftrace=full -ffree-line-length-huge

.F90.o: 
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${FFLAGS} ${FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/$@ $<


${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

m_TestSieve.o:    ${PETSC_ARCH}
TestSieve.o:      ${PETSC_ARCH} m_TestSieve.o ${PETSC_ARCH}/m_TestSieve.o
TestSieve:        ${PETSC_ARCH} m_TestSieve.o ${PETSC_ARCH}/m_TestSieve.o \
                TestSieve.o   ${PETSC_ARCH}/TestSieve.o 
	-${FLINKER} -o ${PETSC_ARCH}/TestSieve  ${PETSC_ARCH}/TestSieve.o ${PETSC_ARCH}/m_TestSieve.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

m_TestSieve3D.o:  ${PETSC_ARCH}
TestSieve3D.o:    ${PETSC_ARCH} m_TestSieve3D.o ${PETSC_ARCH}/m_TestSieve3D.o
TestSieve3D:      ${PETSC_ARCH} m_TestSieve3D.o ${PETSC_ARCH}/m_TestSieve3D.o \
                  TestSieve3D.o   ${PETSC_ARCH}/TestSieve3D.o 
	-${FLINKER} -o ${PETSC_ARCH}/TestSieve3D  ${PETSC_ARCH}/TestSieve3D.o ${PETSC_ARCH}/m_TestSieve3D.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

TestNodeSet:      ${PETSC_ARCH} TestNodeSet.o ${PETSC_ARCH}/TestNodeSet.o
	-${FLINKER} -o ${PETSC_ARCH}/TestNodeSet ${PETSC_ARCH}/TestNodeSet.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestLocal:        ${PETSC_ARCH} TestLocal.o ${PETSC_ARCH}/TestLocal.o
	-${FLINKER} -o ${PETSC_ARCH}/TestLocal ${PETSC_ARCH}/TestLocal.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestAssembly:     ${PETSC_ARCH} TestAssembly.o ${PETSC_ARCH}/TestAssembly.o
	-${FLINKER} -o ${PETSC_ARCH}/TestAssembly ${PETSC_ARCH}/TestAssembly.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestSectionInt:   ${PETSC_ARCH} TestSectionInt.o ${PETSC_ARCH}/TestSectionInt.o
	-${FLINKER} -o ${PETSC_ARCH}/TestSectionInt ${PETSC_ARCH}/TestSectionInt.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestCoordinates:   ${PETSC_ARCH} TestCoordinates.o ${PETSC_ARCH}/TestCoordinates.o
	-${FLINKER} -o ${PETSC_ARCH}/TestCoordinates ${PETSC_ARCH}/TestCoordinates.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}
