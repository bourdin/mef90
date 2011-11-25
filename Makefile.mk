dir:
	mkdir -p ${OBJDIR}
	mkdir -p ${OBJDIR_PROJ}
	mkdir -p ${MEF90DIR}

clean::
	@rm -Rf ${OBJDIR_PROJ}/*o ${OBJDIR_PROJ}/*mod
	@rmdir ${OBJDIR_PROJ}

debug::
	@echo MEF90_FC_FLAGS is ${MEF90_FC_FLAGS}
	@echo FC is ${FC}
	@echo FC_FLAGS is ${FC_FLAGS}
	@echo PETSC_MAKE_STOP_ON_ERROR is ${PETSC_MAKE_STOP_ON_ERROR}
	@echo FCPPFLAGS is ${FCPPFLAGS}
	@echo FC_MODULE_FLAG is ${FC_MODULE_FLAG}
	@echo FC_MODULE_OUTPUT_FLAG is ${FC_MODULE_OUTPUT_FLAG}


