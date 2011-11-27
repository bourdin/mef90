BIN = PrepVarFracNG
#Order is important
SRCS = m_PrepVarFracNG.F90 PrepVarFracNG.F90 

MEF90_INCS += -I${OBJDIR}
#FFLAGS is used in FLINKER def in petsc/rules
#+= si pas defini avant ne rajoute rien
FFLAGS = ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
FFLAGS += -L${LIBDIR} -lMEF90 
FFLAGS+= -lVar_Struct 


include ../Makefile.inc
