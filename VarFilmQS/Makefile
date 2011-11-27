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
TAO =True

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
