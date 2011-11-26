BIN = PrepVarFrac
#Order is important
SRCS = m_PrepVarFrac.F90 PrepVarFrac.F90 

MEF90_INCS += -I../MEF90 -I../m_VarFrac_Struct 
#FFLAGS is used in FLINKER def in petsc/rules
#+= si pas defini avant ne rajoute rien
FFLAGS = ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
FFLAGS += -L../MEF90 -lMEF90 
FFLAGS+= -L../m_VarFrac_Struct -lVar_Struct 

#librairies : conventions minsuule et pas d underscore 

include ../Makefile.inc
