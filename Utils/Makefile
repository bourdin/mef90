MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 

all: Partitioner EXOGetBB

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables

ifneq (,$(findstring xlf,$(FC)))
	FC_MODULE_OUTPUT_FLAG=-qmoddir=
endif

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS=-g ${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} -I./${PETSC_ARCH} -I../MEF90/${PETSC_ARCH} -I${M_VARFRACSTRUCTDIR} ${FC_FLAGS} -ffree-line-length-none -ffixed-line-length-none  -ffree-form
else
	MEF90_FC_FLAGS=-g ${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} -I./${PETSC_ARCH} -I../MEF90/${PETSC_ARCH} -I${M_VARFRACSTRUCTDIR} ${FC_FLAGS}
endif

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

Partitioner:   ${PETSC_ARCH} Partitioner.o ${PETSC_ARCH}/Partitioner.o
	-${FLINKER} -o ${PETSC_ARCH}/Partitioner ${PETSC_ARCH}/Partitioner.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

EXOGetBB:   ${PETSC_ARCH} EXOGetBB.o ${PETSC_ARCH}/EXOGetBB.o
	-${FLINKER} -o ${PETSC_ARCH}/EXOGetBB ${PETSC_ARCH}/EXOGetBB.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}
