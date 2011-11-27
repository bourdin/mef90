BIN = SimplePoisson 
#Order is important
SRCS = \
	   m_Poisson.F90  \
	   m_MEF_Integrate.F90  \
	   m_TransientHeat.F90 \
	   SimplePoisson.F90   \
#	   TransientHeat.F90  

MEF90_INCS += -I${OBJDIR} 
#FFLAGS is used in FLINKER def in petsc/rules
#+= si pas defini avant ne rajoute rien
FFLAGS = ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
FFLAGS += -L${LIBDIR} -lMEF90 
FFLAGS += -lVar_Struct 

ifneq (,$(findstring xlf,$(FC)))
	FC_D_FLAG=-WF,-D
else
  	FC_D_FLAG=-D
endif

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FLAGS += ${FC_D_FLAG}PB_${DIM} 
else
	MEF90_FLAGS += -DPB_${DIM} 
endif

#order of includes is important
include ../Makefile.inc  
