BINDIR = ${MEF90_DIR}/bin/${PETSC_ARCH}
OBJDIR = ${MEF90_DIR}/obj/${PETSC_ARCH}

MEF90DIR  = ${OBJDIR}/MEF90/
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o      $(MEF90DIR)/m_MEF_IO.o 

M_VARFRACSTRUCTDIR = ${OBJDIR}/m_VarFrac_Struct
M_VARFRACSTRUCTOBJS = ${M_VARFRACSTRUCTDIR}/m_VarFrac_Struct.o

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


ifneq (,$(findstring xlf,$(FC)))
	FC_MODULE_OUTPUT_FLAG=-qmoddir=
endif

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${OBJDIR}/${PROJ} ${FC_FLAGS}  -ffixed-line-length-none  -ffree-line-length-none -ffree-form -I${MEF90DIR}
else
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${OBJDIR}/${PROJ} ${FC_FLAGS} -I${MEF90DIR}
endif


