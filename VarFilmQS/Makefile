BIN = VarFilmQS 
#Order is important
SRCS = m_VarFilmQS_Types.F90  \
	   m_VarFilmQS_Post.F90  \
	   m_VarFilmQS_V.F90 \
	   m_VarFilmQS_W.F90  \
	   m_VarFilmQS_U.F90  \
	   m_VarFilmQS.F90  \
	   VarFilmQS.F90  


MEF90_INCS += -I${OBJDIR} ${TAO_INCLUDE} 
#FFLAGS is used in FLINKER def in petsc/rules
#+= si pas defini avant ne rajoute rien
FFLAGS = ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
FFLAGS += -L${LIBDIR} -lMEF90 
FFLAGS += -lVar_Struct 
FFLAGS += ${TAO_LIB} ${TAO_FORTRAN_LIB} 

#ifneq (,$(findstring xlf,$(FC)))
#	FC_D_FLAG=-WF,-D
#else
#  	FC_D_FLAG=-D
#endif

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FLAGS += ${FC_D_FLAG}PB_${DIM} ${FC_D_FLAG}WITH_TAO 
else
	MEF90_FLAGS += -DPB_${DIM} -DWITH_TAO 
endif

#order of includes is important
include ../Makefile.inc  
include ${TAO_DIR}/bmake/tao_common_variables
