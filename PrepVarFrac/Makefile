M_VARFRACSTRUCTDIR  = ../m_VarFrac_Struct/${PETSC_ARCH}
M_VARFRACSTRUCTOBJS = ${M_VARFRACSTRUCTDIR}/m_VarFrac_Struct.o

all: PrepVarFrac PrepVarFracPolyCrystal

include ../Makefile.inc  


.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} -I${M_VARFRACSTRUCTDIR} ${FCPPFLAGS} -o ${PETSC_ARCH}/$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

PrepVarFrac: ${PETSC_ARCH} PrepVarFrac.o ${PETSC_ARCH}/PrepVarFrac.o m_PrepVarFrac.o ${PETSC_ARCH}/m_PrepVarFrac.o
	${FLINKER} -o ${PETSC_ARCH}/PrepVarFrac ${PETSC_ARCH}/PrepVarFrac.o ${PETSC_ARCH}/m_PrepVarFrac.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

PrepVarFracPolyCrystal: ${PETSC_ARCH} PrepVarFracPolyCrystal.o ${PETSC_ARCH}/PrepVarFracPolyCrystal.o m_PrepVarFrac.o ${PETSC_ARCH}/m_PrepVarFrac.o
	${FLINKER} -o ${PETSC_ARCH}/PrepVarFracPolyCrystal ${PETSC_ARCH}/PrepVarFracPolyCrystal.o ${PETSC_ARCH}/m_PrepVarFrac.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

PrepVarFrac.o: ${PETSC_ARCH} m_PrepVarFrac.o ${PETSC_ARCH}/m_PrepVarFrac.o
PrepVarFracPolyCrystal.o: ${PETSC_ARCH} m_PrepVarFrac.o ${PETSC_ARCH}/m_PrepVarFrac.o

m_PrepVarFrac.o: ${PETSC_ARCH}

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
