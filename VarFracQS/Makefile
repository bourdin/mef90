PROJ      = VarFracQS
OBJDIR_PROJ = ${OBJDIR}/${PROJ}

all: VarFracQS${DIM} VarFracQSRecomputeEnergy${DIM}

include ../Makefile.inc  
include ${TAO_DIR}/bmake/tao_common

ifneq (,$(findstring xlf,$(FC)))
	FC_D_FLAG=-WF,-D
else
  	FC_D_FLAG=-D
endif

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS=${FC_D_FLAG}PB_${DIM} ${FC_D_FLAG}WITH_TAO ${FC_MODULE_OUTPUT_FLAG}${OBJDIR_PROJ} -I./${OBJDIR_PROJ} -I${MEF90DIR} -I${M_VARFRACSTRUCTDIR} ${FC_FLAGS} -ffree-line-length-none -ffixed-line-length-none  -ffree-form

else
	MEF90_FC_FLAGS=-DPB_${DIM} -DWITH_TAO ${FC_MODULE_OUTPUT_FLAG}${OBJDIR_PROJ} -I./${OBJDIR_PROJ} -I${MEF90DIR} -I${M_VARFRACSTRUCTDIR} ${FC_FLAGS}
endif

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${OBJDIR_PROJ}/${DIM}$@ $<
dir:
	@mkdir -p ${OBJDIR_PROJ}
	@mkdir -p ${BINDIR}

VarFracQS${DIM}: dir VarFracQS.o   ${OBJDIR_PROJ}/${DIM}VarFracQS.o \
                                      m_VarFracQS_Types.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Types.o \
                                      m_VarFracQS_U.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_U.o \
                                      m_VarFracQS_V.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_V.o \
                                      m_VarFracQS.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS.o \
                                      m_VarFracQS_Post.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Post.o
	${FLINKER} -o ${BINDIR}/VarFracQS${DIM} ${OBJDIR_PROJ}/${DIM}VarFracQS.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Types.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_U.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_V.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS.o  ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Post.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

VarFracQS.o:  dir m_VarFracQS.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS.o \
                                      m_VarFracQS_Types.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Types.o \
                                      m_VarFracQS_U.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_U.o \
                                      m_VarFracQS_V.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_V.o \
                                      m_VarFracQS_Post.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Post.o

VarFracQSRecomputeEnergy${DIM}: dir VarFracQSRecomputeEnergy.o   ${OBJDIR_PROJ}/${DIM}VarFracQSRecomputeEnergy.o \
                                          m_VarFracQS_Types.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Types.o \
                                          m_VarFracQS_U.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_U.o \
                                          m_VarFracQS_V.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_V.o \
                                          m_VarFracQS.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS.o \
                                          m_VarFracQS_Post.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Post.o
	${FLINKER} -o ${BINDIR}/VarFracQSRecomputeEnergy${DIM} ${OBJDIR_PROJ}/${DIM}VarFracQSRecomputeEnergy.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Types.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_U.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_V.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS.o  ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Post.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

VarFracQSRecomputeEnergy.o:   dir m_VarFracQS.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS.o \
                                              m_VarFracQS_Types.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Types.o \
                                              m_VarFracQS_U.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_U.o \
                                              m_VarFracQS_V.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_V.o \
                                              m_VarFracQS_Post.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Post.o

m_VarFracQS.o: dir m_VarFracQS_Types.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Types.o \
                                      m_VarFracQS_U.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_U.o \
                                      m_VarFracQS_V.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_V.o \
                                      m_VarFracQS_Post.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Post.o
m_VarFracQS_U.o:  dir m_VarFracQS_Types.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Types.o \
                                m_VarFracQS_Post.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Post.o
m_VarFracQS_V.o:  dir m_VarFracQS_Types.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Types.o \
                                m_VarFracQS_Post.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Post.o
m_VarFracQS_Post.o:  dir m_VarFracQS_Types.o ${OBJDIR_PROJ}/${DIM}m_VarFracQS_Types.o
m_VarFracQS_Types.o: dir

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}

debug::
	@echo MEF90_FC_FLAGS is ${MEF90_FC_FLAGS}
	@echo FC is ${FC}
	@echo FC_FLAGS is ${FC_FLAGS}
	@echo PETSC_MAKE_STOP_ON_ERROR is ${PETSC_MAKE_STOP_ON_ERROR}
	@echo FCPPFLAGS is ${FCPPFLAGS}
	@echo FC_MODULE_FLAG is ${FC_MODULE_FLAG}
	@echo FC_MODULE_OUTPUT_FLAG is ${FC_MODULE_OUTPUT_FLAG}
