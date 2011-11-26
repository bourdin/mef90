BIN = PrepVarFrac
#Order is important
SRCS =	\
	m_MEF_Parameters.F90      \
	m_MEF_Utils.F90           \
	m_MEF_LinAlg.F90          \
	m_MEF_Types.F90	      \
	m_MEF_MPI.F90             \
	m_MEF_Elements.F90        \
	m_MEF_BoundaryElements.F90\
	m_MEF_EXO.F90             \
	m_MEF_Sieve.F90           \
	m_MEF_IO.F90             \
	m_MEF90.F90
SRCS +=	m_VarFrac_Struct.F90
SRCS += m_PrepVarFrac.F90 PrepVarFrac.F90 

VPATH = ../MEF90 ../m_VarFrac_Struct 
#MEF90_INCS += -I../MEF90 -I../m_VarFrac_Struct 
#FFLAGS is used in FLINKER def in petsc/rules
#+= si pas defini avant ne rajoute rien
FFLAGS = ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
#FFLAGS += -L../MEF90 -lMEF90 
#FFLAGS+= -L../m_VarFrac_Struct -lVar_Struct 

#librairies : conventions minsuule et pas d underscore 

include ../Makefile.inc
