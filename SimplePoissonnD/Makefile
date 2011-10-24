MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 

all: SimplePoisson${DIM} TransientHeat${DIM}

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS=-DPB_${DIM} ${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} -I./${PETSC_ARCH} -I../MEF90/${PETSC_ARCH} ${FC_FLAGS} -ffree-line-length-none -ffixed-line-length-none  -ffree-form
else
	MEF90_FC_FLAGS=-DPB_${DIM} ${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} -I./${PETSC_ARCH} -I../MEF90/${PETSC_ARCH} ${FC_FLAGS}
endif


.F90.o: 
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/${DIM}$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}
	
SimplePoisson${DIM}: ${PETSC_ARCH} SimplePoisson.o   ${PETSC_ARCH}/${DIM}SimplePoisson.o \
                             m_Poisson.o ${PETSC_ARCH/${DIM}m_Poisson.o 
	echo ${FLINKER} -o ${PETSC_ARCH}/SimplePoisson${DIM} ${PETSC_ARCH}/${DIM}SimplePoisson.o  ${PETSC_ARCH}/${DIM}m_Poisson.o  ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
	${FLINKER} -o ${PETSC_ARCH}/SimplePoisson${DIM} ${PETSC_ARCH}/${DIM}SimplePoisson.o  ${PETSC_ARCH}/${DIM}m_Poisson.o  ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

TransientHeat${DIM}: ${PETSC_ARCH} TransientHeat.o   ${PETSC_ARCH}/${DIM}TransientHeat.o \
                             m_TransientHeat.o ${PETSC_ARCH}/${DIM}TransientHeat.o \
                             m_Poisson.o ${PETSC_ARCH}/${DIM}m_Poisson.o
	${FLINKER} -o ${PETSC_ARCH}/TransientHeat${DIM}  ${PETSC_ARCH}/${DIM}m_Poisson.o  ${PETSC_ARCH}/${DIM}TransientHeat.o  ${PETSC_ARCH}/${DIM}m_TransientHeat.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

SimplePoisson.o:   ${PETSC_ARCH} m_Poisson.o ${PETSC_ARCH}/${DIM}m_Poisson.o

m_Poisson.o: ${PETSC_ARCH}

TransientHeat.o:   ${PETSC_ARCH} m_Poisson.o ${PETSC_ARCH}/${DIM}m_Poisson.o \
                             m_TransientHeat.o ${PETSC_ARCH}/${DIM}m_TransientHeat.o

m_TransientHeat.o: ${PETSC_ARCH} m_Poisson.o ${PETSC_ARCH}/${DIM}m_Poisson.o

TSPoisson${DIM}: ${PETSC_ARCH} TSPoisson.o   ${PETSC_ARCH}/${DIM}TSPoisson.o m_TSPoisson.o ${PETSC_ARCH}/${DIM}m_TSPoisson.o 
	${FLINKER} -o ${PETSC_ARCH}/TSPoisson${DIM} ${PETSC_ARCH}/${DIM}TSPoisson.o  ${PETSC_ARCH}/${DIM}m_TSPoisson.o  ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TSPoisson.o:   ${PETSC_ARCH} m_TSPoisson.o ${PETSC_ARCH}/${DIM}m_TSPoisson.o

m_TSPoisson.o: ${PETSC_ARCH}


clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}

