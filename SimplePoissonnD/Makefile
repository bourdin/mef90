PROJ = Heat
OBJDIR_PROJ = ${OBJDIR}/${PROJ}

all: SimplePoisson${DIM} TransientHeat${DIM}

include ../Makefile.inc

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS=-DPB_${DIM} ${FC_MODULE_OUTPUT_FLAG}${OBJDIR_PROJ} -I./${OBJDIR_PROJ} -I${MEF90DIR} ${FC_FLAGS} -ffree-line-length-none -ffixed-line-length-none  -ffree-form
else
	MEF90_FC_FLAGS=-DPB_${DIM} ${FC_MODULE_OUTPUT_FLAG}${OBJDIR_PROJ} -I./${OBJDIR_PROJ} -I${MEF90DIR} ${FC_FLAGS} 
endif


.F90.o: 
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${OBJDIR_PROJ}/${DIM}$@ $<

dir:
	@mkdir -p ${OBJDIR_PROJ}
	@mkdir -p ${BINDIR}
	
SimplePoisson${DIM}: dir SimplePoisson.o  ${OBJDIR_PROJ}/${DIM}SimplePoisson.o m_Poisson.o ${OBJDIR_PROJ}/${DIM}m_Poisson.o 
	${FLINKER} -o ${BINDIR}/SimplePoisson${DIM} ${OBJDIR_PROJ}/${DIM}SimplePoisson.o  ${OBJDIR_PROJ}/${DIM}m_Poisson.o  ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

TransientHeat${DIM}: dir TransientHeat.o   ${OBJDIR_PROJ}/${DIM}TransientHeat.o \
                             m_TransientHeat.o ${OBJDIR_PROJ}/${DIM}TransientHeat.o \
                             m_MEF_Integrate.o ${OBJDIR_PROJ}/${DIM}m_MEF_Integrate.o
	${FLINKER} -o ${BINDIR}/TransientHeat${DIM}  ${OBJDIR_PROJ}/${DIM}m_Poisson.o  ${OBJDIR_PROJ}/${DIM}TransientHeat.o  ${OBJDIR_PROJ}/${DIM}m_TransientHeat.o ${OBJDIR_PROJ}/${DIM}m_MEF_Integrate.o ${MEF90OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

SimplePoisson.o: dir m_Poisson.o ${OBJDIR_PROJ}/${DIM}m_Poisson.o

m_Poisson.o: dir

m_MEF_Integrate.o: dir 

TransientHeat.o:   dir m_Poisson.o ${OBJDIR_PROJ}/${DIM}m_Poisson.o \
                             m_TransientHeat.o ${OBJDIR_PROJ}/${DIM}m_TransientHeat.o

m_TransientHeat.o: dir m_Poisson.o ${OBJDIR_PROJ}/${DIM}m_Poisson.o \
                             m_MEF_Integrate.o ${OBJDIR_PROJ}/${DIM}m_MEF_Integrate.o


clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}


