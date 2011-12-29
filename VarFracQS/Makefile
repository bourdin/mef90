BIN = VarFracQS VarFracQSRecomputeEnergy 
#Order is important
SRCS = \
	   m_VarFracQS_Types.F90  \
	   m_VarFracQS_Post.F90  \
	   m_VarFracQS_V.F90 \
	   m_VarFracQS_U.F90  \
	   m_VarFracQS.F90  \

#Compilation takes DIM option 
############### Test if Options are defined properly 
#ugly hack 
ifndef DIM
True:
	@echo The DIM option is mandatory for VarFracQS compilation 
	@echo ------------------ Stop Compilation ---------------- 
endif

#For using heat module computation should be launched with OPTS=Heat
ifdef OPTS
ifneq (${OPTS}, Heat)
True:
	@echo The OPTS variable if defined can only take the value Heat 
	@echo ------------------ Stop Compilation ---------------- 
else 
SRCS += m_VarFracQS_T.F90  
endif
endif 


ifeq (${OPTS}, Heat)
	ifneq (,$(findstring gcc,$(PETSC_ARCH)))
		MEF90_FLAGS +=  ${FC_D_FLAG}HEAT 
	else
		MEF90_FLAGS += -DHEAT 
	endif
endif 

TAO = True 
MEF90_INCS += -I${OBJDIR} ${TAO_INCLUDE} 



include ../Makefile.inc  
