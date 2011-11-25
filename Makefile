SUBDIRS = MEF90 m_VarFrac_Struct PrepVarFrac VarFracQS

all:
	for d in $(SUBDIRS); do $(MAKE) -C $$d DIM=${DIM}; done

