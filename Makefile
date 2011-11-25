all: PrepVarFrac VarFracQS
iallFilm: PrepFilm VarFilmQS 
allHeat: PrepHeat VarFracQS

include Makefile.mk

VarFilmQS::  MEF90 m_Film_Struct
	$(MAKE) -C VarFilm

PrepVarFilm:: MEF90 m_Film_Struct
	$(MAKE) -C PrepVarFilm

VarFracQS:: MEF90 m_VarFrac_Struct
	$(MAKE) -C VarFracQS DIM=2D
	$(MAKE) -C VarFracQS DIM=3D

PrepVarFrac:: MEF90 m_VarFrac_Struct
	$(MAKE) -C PrepVarFrac 

PrepVarFracNG:: MEF90 m_VarFrac_Struct
	$(MAKE) -C PrepVarFracNG

PrepHeat:: MEF90 m_VarFrac_Struct Heat
	$(MAKE) -C PrepVarFracNG DIM=2D WITH_HEAT=WITH_HEAT

Heat:: MEF90
	${Make} -C SimplePoissonnD DIM=2D

m_VarFrac_Struct:: 
	${MAKE} -C m_VarFrac_Struct

m_Film_Struct:: MEF90
	${MAKE} -C m_Film_Struct

MEF90::
	${MAKE} -C MEF90 


clean::
	@rm -Rf ${OBJDIR}/*
	@rmdir ${OBJDIR}
	@rm -Rf ${BINDIR}/*
	@rmdir ${BINDIR}

