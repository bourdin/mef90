PROJ      = m_Film_Struct
OBJDIR_PROJ = ${OBJDIR}/${PROJ}
all: ${PROJ}.o
include ../Makefile.inc 


.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${OBJDIR_PROJ}/$@ $<

dir:
	@mkdir -p ${OBJDIR_PROJ}

$(PROJ).o: dir

clean::
	@rm -rf ${PETSC_ARCH}/
	@rmdir ${PETSC_ARCH}

debug::
	@echo MEF90_FC_FLAGS is ${MEF90_FC_FLAGS}
	@echo FC is ${FC}
	@echo FC_FLAGS is ${FC_FLAGS}
	@echo PETSC_MAKE_STOP_ON_ERROR is ${PETSC_MAKE_STOP_ON_ERROR}
	@echo FCPPFLAGS is ${FCPPFLAGS}
	@echo FC_MODULE_FLAG is ${FC_MODULE_FLAG}
	@echo FC_MODULE_OUTPUT_FLAG is ${FC_MODULE_OUTPUT_FLAG}

