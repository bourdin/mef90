BIN = VarFracQS
#Order is important
SRCS = \
	   m_VarFracQS_Types.F90  \
	   m_VarFracQS_Post.F90  \
	   m_VarFracQS_V.F90 \
	   m_VarFracQS_U.F90  \
	   m_VarFracQS.F90  \
	   VarFracQS.F90  \
	   VarFracQSRecomputeEnergy.F90

#include ${TAO_DIR}/bmake/tao_common
#

MEF90_INCS += -I../MEF90 -I../m_VarFrac_Struct 
#FFLAGS is used in FLINKER def in petsc/rules
#+= si pas defini avant ne rajoute rien
FFLAGS = ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
FFLAGS += -L../MEF90 -lMEF90 
FFLAGS+= -L../m_VarFrac_Struct -lVar_Struct 

ifneq (,$(findstring xlf,$(FC)))
	FC_D_FLAG=-WF,-D
else
  	FC_D_FLAG=-D
endif

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS += ${FC_D_FLAG}PB_${DIM} ${FC_D_FLAG}WITH_TAO 
else
	MEF90_FC_FLAGS += -DPB_${DIM} -DWITH_TAO 
endif


include ../Makefile.inc  
