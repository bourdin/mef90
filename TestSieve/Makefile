MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 

all: TestSieve TestSieve3D TestNodeSet TestLocal TestAssembly TestSectionInt

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables

ifeq (${FC},gfortran)
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} -I. -I../MEF90/${PETSC_ARCH} ${FC_FLAGS} --free-line-length-none  -ffree-form
else
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG} ${PETSC_ARCH} -I. -I../MEF90/${PETSC_ARCH} ${FC_FLAGS}
endif

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/$@ $<

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

TestDuplicate:     ${PETSC_ARCH} TestDuplicate.o ${PETSC_ARCH}/TestDuplicate.o
	-${FLINKER} -o ${PETSC_ARCH}/TestDuplicate   ${PETSC_ARCH}/TestDuplicate.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestSectionInt:   ${PETSC_ARCH} TestSectionInt.o ${PETSC_ARCH}/TestSectionInt.o
	-${FLINKER} -o ${PETSC_ARCH}/TestSectionInt ${PETSC_ARCH}/TestSectionInt.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestMeshDistribute:   ${PETSC_ARCH} TestMeshDistribute.o ${PETSC_ARCH}/TestMeshDistribute.o
	-${FLINKER} -o ${PETSC_ARCH}/TestMeshDistribute ${PETSC_ARCH}/TestMeshDistribute.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestCoordinates:   ${PETSC_ARCH} TestCoordinates.o ${PETSC_ARCH}/TestCoordinates.o
	-${FLINKER} -o ${PETSC_ARCH}/TestCoordinates ${PETSC_ARCH}/TestCoordinates.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestScatter:   ${PETSC_ARCH} TestScatter.o ${PETSC_ARCH}/TestScatter.o
	-${FLINKER} -o ${PETSC_ARCH}/TestScatter ${PETSC_ARCH}/TestScatter.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestSectionFromVec:   ${PETSC_ARCH} TestSectionFromVec.o ${PETSC_ARCH}/TestSectionFromVec.o
	-${FLINKER} -o ${PETSC_ARCH}/TestSectionFromVec ${PETSC_ARCH}/TestSectionFromVec.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestUpdate:   ${PETSC_ARCH} TestUpdate.o ${PETSC_ARCH}/TestUpdate.o
	-${FLINKER} -o ${PETSC_ARCH}/TestUpdate ${PETSC_ARCH}/TestUpdate.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

Partitioner:   ${PETSC_ARCH} Partitioner.o ${PETSC_ARCH}/Partitioner.o
	-${FLINKER} -o ${PETSC_ARCH}/Partitioner ${PETSC_ARCH}/Partitioner.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestUpdateRestrict:   ${PETSC_ARCH} TestUpdateRestrict.o ${PETSC_ARCH}/TestUpdateRestrict.o
	-${FLINKER} -o ${PETSC_ARCH}/TestUpdateRestrict ${PETSC_ARCH}/TestUpdateRestrict.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestFibration:   ${PETSC_ARCH} TestFibration.o ${PETSC_ARCH}/TestFibration.o
	-${FLINKER} -o ${PETSC_ARCH}/TestFibration ${PETSC_ARCH}/TestFibration.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestFibrationVertex:   ${PETSC_ARCH} TestFibrationVertex.o ${PETSC_ARCH}/TestFibrationVertex.o
	-${FLINKER} -o ${PETSC_ARCH}/TestFibrationVertex ${PETSC_ARCH}/TestFibrationVertex.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestSectionToArray:   ${PETSC_ARCH} TestSectionToArray.o ${PETSC_ARCH}/TestSectionToArray.o
	-${FLINKER} -o ${PETSC_ARCH}/TestSectionToArray ${PETSC_ARCH}/TestSectionToArray.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestSaveMesh:   ${PETSC_ARCH} TestSaveMesh.o ${PETSC_ARCH}/TestSaveMesh.o
	-${FLINKER} -o ${PETSC_ARCH}/TestSaveMesh ${PETSC_ARCH}/TestSaveMesh.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

FindMPIProblem:   ${PETSC_ARCH} FindMPIProblem.o ${PETSC_ARCH}/FindMPIProblem.o
	-${FLINKER} -o ${PETSC_ARCH}/FindMPIProblem ${PETSC_ARCH}/FindMPIProblem.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestH5:   ${PETSC_ARCH} TestH5.o ${PETSC_ARCH}/TestH5.o
	-${FLINKER} -o ${PETSC_ARCH}/TestH5 ${PETSC_ARCH}/TestH5.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestMeshLoad:   ${PETSC_ARCH} TestMeshLoad.o ${PETSC_ARCH}/TestMeshLoad.o
	-${FLINKER} -o ${PETSC_ARCH}/TestMeshLoad ${PETSC_ARCH}/TestMeshLoad.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	
TestLabel:   ${PETSC_ARCH} TestLabel.o ${PETSC_ARCH}/TestLabel.o
		-${FLINKER} -o ${PETSC_ARCH}/TestLabel ${PETSC_ARCH}/TestLabel.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestCorruption:   ${PETSC_ARCH} TestCorruption.o m_TestCorruption.o ${PETSC_ARCH}/TestCorruption.o ${PETSC_ARCH}/m_TestCorruption.o
	-${FLINKER} -o ${PETSC_ARCH}/TestCorruption ${PETSC_ARCH}/TestCorruption.o ${PETSC_ARCH}/m_TestCorruption.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
m_TestCorruption.o: ${PETSC_ARCH}

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}
