BIN = VarFracQS VarFracQSRecomputeEnergy 
#Order is important
SRCS = \
	   m_VarFracQS_Types.F90  \
	   m_VarFracQS_Post.F90  \
	   m_VarFracQS_V.F90 \
	   m_VarFracQS_U.F90  \
	   m_VarFracQS.F90  \
	   VarFracQS.F90  \
#	   VarFracQSRecomputeEnergy.F90

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
