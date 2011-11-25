MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o      $(MEF90DIR)/m_MEF_IO.o 

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


ifneq (,$(findstring xlf,$(FC)))
	FC_MODULE_OUTPUT_FLAG=-qmoddir=
endif

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} ${FC_FLAGS} -ffixed-line-length-none  -ffree-line-length-none -ffree-form -I../MEF90/${PETSC_ARCH}
else
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} ${FC_FLAGS} -I../MEF90/${PETSC_ARCH}
endif


