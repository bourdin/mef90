BIN = VarFracQS 
#Order is important
SRCS = \
	   m_VarFracQS_Types.F90  \
	   m_VarFracQS_Post.F90  \
	   m_VarFracQS_V.F90 \
	   m_VarFracQS_U.F90  \
	   m_VarFracQS.F90  \
	   VarFracQS.F90  #\
#	   VarFracQSRecomputeEnergy.F90


#MEF90_INCS += -I../MEF90 -I../m_VarFrac_Struct -I/opt/HPC/tao-lite -I/opt/HPC/tao-lite/include 
MEF90_INCS += -I${OBJDIR} ${TAO_INCLUDE} 
#FFLAGS is used in FLINKER def in petsc/rules
#+= si pas defini avant ne rajoute rien
FFLAGS = ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
FFLAGS += -L${LIBDIR} -lMEF90 
FFLAGS += -lVar_Struct 
FFLAGS += ${TAO_LIB} ${TAO_FORTRAN_LIB} 
#FFLAGS += -L/opt/HPC/tao-lite/lib/${PETSC_ARCH} -ltaofortran -ltaopetsc -ltao

ifneq (,$(findstring xlf,$(FC)))
	FC_D_FLAG=-WF,-D
else
  	FC_D_FLAG=-D
endif

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FLAGS += ${FC_D_FLAG}PB_${DIM} ${FC_D_FLAG}WITH_TAO 
else
	MEF90_FLAGS += -DPB_${DIM} -DWITH_TAO 
endif

#order of includes is important
include ../Makefile.inc  
include ${TAO_DIR}/bmake/tao_common_variables
