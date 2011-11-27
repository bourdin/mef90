BIN =ELast 
#Order is important
SRCS = m_Elast.F90\
	   Elast.F90

TAO = True

MEF90_INCS += -I${OBJDIR} ${TAO_INCLUDE} 

ifneq (,$(findstring xlf,$(FC)))
	FC_D_FLAG=-WF,-D
else
  	FC_D_FLAG=-D
endif

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FLAGS += ${FC_D_FLAG}PB_${DIM} ${FC_D_FLAG}WITH_TAO 
else
	MEF90_FLAGS += -DPB_${DIM} -DWITH_TAO 
endif

include ../Makefile.inc  
