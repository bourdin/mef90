PROJ = PrepFilm
OBJDIR_PROJ = ${OBJDIR}/${PROJ}

M_FILMSTRUCTDIR  = ${OBJDIR}/m_Film_Struct
M_FILMSTRUCTOBJS = ${M_FILMSTRUCTDIR}/m_Film_Struct.o

all: PrepFilm

include ../Makefile.inc

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${OBJDIR_PROJ} -I./${OBJDIR_PROJ} -I${MEF90DIR} -I${PETSC_ARCH} -I${M_FILMSTRUCTDIR} ${FC_FLAGS}  -ffree-form -ffree-line-length-none -ffixed-line-length-none
else
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${OBJDIR_PROJ} -I./${OBJDIR_PROJ} -I/${MEF90DIR} -I${M_FILMSTRUCTDIR} ${FC_FLAGS}
endif

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${OBJDIR_PROJ}/$@ $<

dir:
	@mkdir -p ${OBJDIR_PROJ}
	@mkdir -p ${BINDIR}

PrepFilm: dir PrepFilm.o ${OBJDIR_PROJ}/PrepFilm.o m_PrepFilm.o ${OBJDIR_PROJ}/m_PrepFilm.o
	${FLINKER} -o ${BINDIR}/PrepFilm ${OBJDIR_PROJ}/PrepFilm.o ${OBJDIR_PROJ}/m_PrepFilm.o ${MEF90OBJS} ${M_FILMSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

PrepFilm.o: dir m_PrepFilm.o ${OBJDIR_PROJ}/m_PrepFilm.o

m_PrepFilm.o: dir

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rm -Rf ${PETSC_ARCH}/
