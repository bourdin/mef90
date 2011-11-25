PROJ = PrepVarFrac
OBJDIR_PROJ = ${OBJDIR}/${PROJ}

all: PrepVarFrac PrepVarFracPolyCrystal

include ../Makefile.inc  


.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} -I${M_VARFRACSTRUCTDIR} ${FCPPFLAGS} -o ${OBJDIR_PROJ}/$@ $<

dir:
	@mkdir -p ${OBJDIR_PROJ}
	@mkdir -p ${BINDIR}

PrepVarFrac: dir PrepVarFrac.o ${OBJDIR_PROJ}/PrepVarFrac.o m_PrepVarFrac.o ${OBJDIR_PROJ}/m_PrepVarFrac.o
	${FLINKER} -o ${BINDIR}/PrepVarFrac ${OBJDIR_PROJ}/PrepVarFrac.o ${OBJDIR_PROJ}/m_PrepVarFrac.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

PrepVarFracPolyCrystal: dir PrepVarFracPolyCrystal.o ${OBJDIR_PROJ}/PrepVarFracPolyCrystal.o m_PrepVarFrac.o ${OBJDIR_PROJ}/m_PrepVarFrac.o
	${FLINKER} -o ${BINDIR}/PrepVarFracPolyCrystal ${OBJDIR_PROJ}/PrepVarFracPolyCrystal.o ${OBJDIR_PROJ}/m_PrepVarFrac.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

PrepVarFrac.o: dir m_PrepVarFrac.o ${OBJDIR_PROJ}/m_PrepVarFrac.o
PrepVarFracPolyCrystal.o: dir m_PrepVarFrac.o ${OBJDIR_PROJ}/m_PrepVarFrac.o

m_PrepVarFrac.o: dir  

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rm -Rf ${PETSC_ARCH}/

debug::
	@echo MEF90_FC_FLAGS is ${MEF90_FC_FLAGS}
	@echo FC is ${FC}
	@echo FC_FLAGS is ${FC_FLAGS}
	@echo PETSC_MAKE_STOP_ON_ERROR is ${PETSC_MAKE_STOP_ON_ERROR}
	@echo FCPPFLAGS is ${FCPPFLAGS}
	@echo FC_MODULE_FLAG is ${FC_MODULE_FLAG}
	@echo FC_MODULE_OUTPUT_FLAG is ${FC_MODULE_OUTPUT_FLAG}
