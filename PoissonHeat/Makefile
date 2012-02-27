BIN = TransientHeat 
#Order is important
SRCS = \
	   m_MEF_Integrate.F90  \
	   m_TransientHeat.F90 \

#Compilation takes DIM option 
############### Test if Options are defined properly 
ifndef DIM
True:
	@echo The DIM option is mandatory for SimplePoisson && TransientHeat compilation 
	@echo --------------------------- Stop Compilation ----------------------------- 
endif

MEF90_INCS += -I${OBJDIR} 

include ../Makefile.inc  
