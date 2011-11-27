BIN =ELast 
#Order is important
SRCS = m_Elast.F90 

MEF90_INCS += -I${OBJDIR} ${TAO_INCLUDE} 
#FFLAGS is used in FLINKER def in petsc/rules
#+= si pas defini avant ne rajoute rien
FFLAGS = ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
FFLAGS += -L${LIBDIR} -lMEF90 
FFLAGS += -lVar_Struct 
FFLAGS += ${TAO_LIB} ${TAO_FORTRAN_LIB} 

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FLAGS += ${FC_D_FLAG}PB_${DIM} ${FC_D_FLAG}WITH_TAO 
else
	MEF90_FLAGS += -DPB_${DIM} -DWITH_TAO 
endif

#order of includes is important
include ../Makefile.inc  
include ${TAO_DIR}/bmake/tao_common_variables
