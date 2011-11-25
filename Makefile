SUBDIRS = MEF90 m_VarFrac_Struct PrepVarFrac 


all: 
	for d in $(SUBDIRS); do $(MAKE) -C $$d DIM=${DIM}; done


clean::
	@rm -Rf ${OBJDIR}/*o ${OBJDIR}/*mod
	@rmdir ${OBJDIR}
	@rm -Rf ${BINDIR}/
	@rmdir ${BINDIR}

