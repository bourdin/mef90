PROJ      = m_Film_Struct
OBJDIR_PROJ = ${OBJDIR}/${PROJ}
all: ${PROJ}.o

include ../Makefile.inc 
include ../Makefile.mk 


.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${OBJDIR_PROJ}/$@ $<

dir:
	@mkdir -p ${OBJDIR_PROJ}

$(PROJ).o: dir


