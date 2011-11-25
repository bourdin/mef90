PROJ = PrepVarFrac
OBJDIR_PROJ = ${OBJDIR}/${PROJ}

all: PrepVarFrac PrepVarFracPolyCrystal

include ../Makefile.inc  
include ../Makefile.mk  


.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} -I${M_VARFRACSTRUCTDIR} ${FCPPFLAGS} -o ${OBJDIR_PROJ}/$@ $<

PrepVarFrac: dir PrepVarFrac.o ${OBJDIR_PROJ}/PrepVarFrac.o m_PrepVarFrac.o ${OBJDIR_PROJ}/m_PrepVarFrac.o
	${FLINKER} -o ${BINDIR}/PrepVarFrac ${OBJDIR_PROJ}/PrepVarFrac.o ${OBJDIR_PROJ}/m_PrepVarFrac.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

PrepVarFracPolyCrystal: dir PrepVarFracPolyCrystal.o ${OBJDIR_PROJ}/PrepVarFracPolyCrystal.o m_PrepVarFrac.o ${OBJDIR_PROJ}/m_PrepVarFrac.o
	${FLINKER} -o ${BINDIR}/PrepVarFracPolyCrystal ${OBJDIR_PROJ}/PrepVarFracPolyCrystal.o ${OBJDIR_PROJ}/m_PrepVarFrac.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

PrepVarFrac.o: dir m_PrepVarFrac.o ${OBJDIR_PROJ}/m_PrepVarFrac.o
PrepVarFracPolyCrystal.o: dir m_PrepVarFrac.o ${OBJDIR_PROJ}/m_PrepVarFrac.o

m_PrepVarFrac.o: dir  

