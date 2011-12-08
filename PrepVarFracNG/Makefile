BIN = PrepVarFracNG
#Order is important
SRCS = m_PrepVarFracNG.F90 PrepVarFracNG.F90 

MEF90_INCS += -I${OBJDIR}

#For using heat module computation should be launched with DIM and WITH_HEAT
ifdef WITH_HEAT
	ifneq (,$(findstring gcc,$(PETSC_ARCH)))
		MEF90_FLAGS +=  ${FC_D_FLAG}WITH_HEAT 
	else
		MEF90_FLAGS += -DWITH_HEAT 
	endif
endif 

include ../Makefile.inc
