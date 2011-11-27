BIN = SimplePoisson TransientHeat 
#Order is important
SRCS = \
	   m_Poisson.F90  \
	   m_MEF_Integrate.F90  \
	   m_TransientHeat.F90 \

#Compilation takes DIM option 

MEF90_INCS += -I${OBJDIR} 

include ../Makefile.inc  
