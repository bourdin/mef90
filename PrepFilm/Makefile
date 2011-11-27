BIN = PrepFilm
#Order is important
SRCS = m_PrepFilm.F90 PrepFilm.F90 

MEF90_INCS += -I${OBJDIR}
#FFLAGS is used in FLINKER def in petsc/rules
FFLAGS = ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
FFLAGS += -L${LIBDIR} -lMEF90 
FFLAGS+= -lVar_Struct 

include ../Makefile.inc

