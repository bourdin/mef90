PROJ = PrepVarFracNG
OBJDIR_PROJ = ${OBJDIR}/${PROJ}

M_HEATDIR  = ${OBJDIR}/Heat
M_HEATOBJS =  ${M_HEATDIR}/${DIM}m_Poisson.o ${M_HEATDIR}/${DIM}m_TransientHeat.o ${M_HEATDIR}/${DIM}m_MEF_Integrate.o

all: PrepVarFracNG

include ../Makefile.inc  
include ../Makefile.mk  

ifneq (,$(findstring xlf,$(FC)))
	FC_MODULE_OUTPUT_FLAG=-qmoddir=
endif

ifdef WITH_HEAT
ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS=-DPB_${DIM} -D${WITH_HEAT} -g ${FC_MODULE_OUTPUT_FLAG}${OBJDIR_PROJ} -I./${OBJDIR_PROJ} -I${MEF90DIR} -I${M_VARFRACSTRUCTDIR} -I${M_HEATDIR} ${FC_FLAGS} -ffree-line-length-none -ffixed-line-length-none  -ffree-form
else
	MEF90_FC_FLAGS=-DPB_${DIM} -D${WITH_HEAT} -g ${FC_MODULE_OUTPUT_FLAG}${OBJDIR_PROJ} -I./${OBJDIR_PROJ} -I${MEF90DIR} -I${M_VARFRACSTRUCTDIR} -I${M_HEATDIR} ${FC_FLAGS}
endif
else
ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS= -g ${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} -I./${OBJDIR_PROJ} -I./${MEF90DIR} -I${M_VARFRACSTRUCTDIR} -I${M_HEATDIR} ${FC_FLAGS} -ffree-line-length-none -ffixed-line-length-none  -ffree-form
else
	MEF90_FC_FLAGS= -g ${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} -I./${OBJDIR_PROJ} -I/${MEF90DIR} -I${M_VARFRACSTRUCTDIR} -I${M_HEATDIR} ${FC_FLAGS}
endif
endif


.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${OBJDIR_PROJ}/$@ $<

ifdef WITH_HEAT
PrepVarFracNG: dir PrepVarFracNG.o ${OBJDIR_PROJ}/PrepVarFracNG.o m_PrepVarFracNG.o ${OBJDIR_PROJ}/m_PrepVarFracNG.o
	${FLINKER} -o ${BINDIR}/PrepVarFracNG ${OBJDIR_PROJ}/PrepVarFracNG.o ${OBJDIR_PROJ}/m_PrepVarFracNG.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${M_HEATOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
else
PrepVarFracNG: dir PrepVarFracNG.o ${OBJDIR_PROJ}/PrepVarFracNG.o m_PrepVarFracNG.o ${OBJDIR_PROJ}/m_PrepVarFracNG.o
	${FLINKER} -o ${BINDIR}/PrepVarFracNG ${OBJDIR_PROJ}/PrepVarFracNG.o ${OBJDIR_PROJ}/m_PrepVarFracNG.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
endif

PrepVarFracNG.o: dir m_PrepVarFracNG.o ${OBJDIR_PROJ}/m_PrepVarFracNG.o

m_PrepVarFracNG.o: dir

debug::
	@echo MEF90_FC_FLAGS is ${MEF90_FC_FLAGS}
	@echo FC is ${FC}
	@echo FC_FLAGS is ${FC_FLAGS}
	@echo PETSC_MAKE_STOP_ON_ERROR is ${PETSC_MAKE_STOP_ON_ERROR}
	@echo FCPPFLAGS is ${FCPPFLAGS}
	@echo FC_MODULE_FLAG is ${FC_MODULE_FLAG}
	@echo FC_MODULE_OUTPUT_FLAG is ${FC_MODULE_OUTPUT_FLAG}
	@echo ${M_HEATDIR}
	@echo ${M_HEATOBJS}
	@echo ${WITH_HEAT}
