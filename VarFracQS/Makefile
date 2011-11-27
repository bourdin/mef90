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

#Compilation takes DIM option 

TAO = True 
MEF90_INCS += -I${OBJDIR} ${TAO_INCLUDE} 

include ../Makefile.inc  
