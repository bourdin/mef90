BIN = VarFracQS VarFracQSRecomputeEnergy 
#Order is important
SRCS = \
	   m_VarFracQS_Types.F90  \
	   m_VarFracQS_Post.F90  \
	   m_VarFracQS_V.F90 \
	   m_VarFracQS_U.F90  \
	   m_VarFracQS.F90  \

#Compilation takes DIM option 
#For using heat module computation should be launched with WITH_HEAT
ifdef WITH_HEAT
	ifneq (,$(findstring gcc,$(PETSC_ARCH)))
		MEF90_FLAGS +=  ${FC_D_FLAG}WITH_HEAT 
	else
		MEF90_FLAGS += -DWITH_HEAT 
	endif
endif 

WITH_HEAT = True
TAO = True 
MEF90_INCS += -I${OBJDIR} ${TAO_INCLUDE} 

include ../Makefile.inc  
