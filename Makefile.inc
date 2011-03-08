###
### a few declarations common to all mef90 projects
###

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} ${FC_FLAGS} -ffixed-line-length-none  -ffree-line-length-none -ffree-form -I${PETSC_ARCH}
else
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG}${PETSC_ARCH} ${FC_FLAGS}
endif

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/$@ $<

debug::
	@echo MEF90_FC_FLAGS is ${MEF90_FC_FLAGS}
	@echo COMPILER is ${COMPILER}
	@echo FC is ${FC}
	@echo FC_FLAGS is ${FC_FLAGS}
	@echo MEF90_FC_FLAGS is ${MEF90_FC_FLAGS}
	@echo PETSC_MAKE_STOP_ON_ERROR is ${PETSC_MAKE_STOP_ON_ERROR}
	@echo FCPPFLAGS is ${FCPPFLAGS}
	@echo FC_MODULE_FLAG is ${FC_MODULE_FLAG}
	@echo FC_MODULE_OUTPUT_FLAG is ${FC_MODULE_OUTPUT_FLAG}
	@echo COMPILER is ${COMPILER}

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

