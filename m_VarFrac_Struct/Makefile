PROJ      = m_VarFrac_Struct

all: ${PROJ}.o 

include ../Makefile.inc  
include ../Makefile.mk 

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${M_VARFRACSTRUCTDIR}/$@ $<

dir:
	@mkdir -p ${M_VARFRACSTRUCTDIR} 

$(PROJ).o: dir 

