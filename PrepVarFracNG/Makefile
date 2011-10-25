MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o      $(MEF90DIR)/m_MEF_IO.o 
M_VARFRACSTRUCTDIR  = ../m_VarFrac_Struct/${PETSC_ARCH}
M_VARFRACSTRUCTOBJS = ${M_VARFRACSTRUCTDIR}/m_VarFrac_Struct.o

all: PrepVarFracNG

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

PrepVarFracNG: ${PETSC_ARCH} PrepVarFracNG.o ${PETSC_ARCH}/PrepVarFracNG.o m_PrepVarFracNG.o ${PETSC_ARCH}/m_PrepVarFracNG.o
	${FLINKER} -o ${PETSC_ARCH}/PrepVarFracNG ${PETSC_ARCH}/PrepVarFracNG.o ${PETSC_ARCH}/m_PrepVarFracNG.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

PrepVarFracNG.o: ${PETSC_ARCH} m_PrepVarFracNG.o ${PETSC_ARCH}/m_PrepVarFracNG.o

m_PrepVarFracNG.o: ${PETSC_ARCH}

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rm -Rf ${PETSC_ARCH}/

debug::
	@echo MEF90_FC_FLAGS is ${MEF90_FC_FLAGS}
	@echo FC is ${FC}
	@echo FC_FLAGS is ${FC_FLAGS}
	@echo PETSC_MAKE_STOP_ON_ERROR is ${PETSC_MAKE_STOP_ON_ERROR}
	@echo FCPPFLAGS is ${FCPPFLAGS}
	@echo FC_MODULE_FLAG is ${FC_MODULE_FLAG}
	@echo FC_MODULE_OUTPUT_FLAG is ${FC_MODULE_OUTPUT_FLAG}
