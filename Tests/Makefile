MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o        \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o      \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o           \
            $(MEF90DIR)/m_MEF_Sieve.o      $(MEF90DIR)/m_MEF_BoundaryElements.o

all: TestBoundaryElement3D TestBoundaryElement2D

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
        MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} -I./${PETSC_ARCH} -I../MEF90/${PETSC_ARCH} ${FC_FLAGS} -ffree-line-length-none -ffixed-line-length-none  -ffree-form
else
        MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} -I./${PETSC_ARCH} -I../MEF90/${PETSC_ARCH} ${FC_FLAGS}
endif

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

TestBoundaryElement3D:      ${PETSC_ARCH} TestBoundaryElement3D.o ${PETSC_ARCH}/TestBoundaryElement3D.o
	-${FLINKER} -o ${PETSC_ARCH}/TestBoundaryElement3D ${PETSC_ARCH}/TestBoundaryElement3D.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestBoundaryElement2D:      ${PETSC_ARCH} TestBoundaryElement2D.o ${PETSC_ARCH}/TestBoundaryElement2D.o
	-${FLINKER} -o ${PETSC_ARCH}/TestBoundaryElement2D ${PETSC_ARCH}/TestBoundaryElement2D.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
