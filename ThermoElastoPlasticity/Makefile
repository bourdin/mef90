all: ThermoElastoPlasticity
runall: runTestMaterials runTestMEF90Ctx 

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
include ${MEF90_DIR}/m_HeatXfer/Makefile.include
include ${MEF90_DIR}/m_DefMech/Makefile.include


ThermoElastoPlasticity.o: ${MEF90_DIR}/ThermoElastoPlasticity/ThermoElastoPlasticity.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $< 

ThermoElastoPlasticity: ThermoElastoPlasticity.o  ${MEF90_DEFMECH_OBJS}  ${MEF90_HEATXFER_OBJS}
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@${RM} $<


MEF90:
	make -C ${MEF90_DIR} MEF90
m_HeatXfer:
	make -C ${MEF90_DIR} m_HeatXfer
m_DEFMECH:
	make -C ${MEF90_DIR} m_DEFMECH

clean::
	-@rm -f TestMaterials TestMEF90Ctx testHeatXfer testLeak
allclean:
	make -C ${MEF90_DIR} clean
	
debug::
	@echo FFLAGS ${FFLAGS} 
	@echo FC_FLAGS ${FC_FLAGS}
	@echo FCFLAGS ${FCFLAGS}
	@echo FCPPFLAGS ${FCPPFLAGS}
	@echo FCPP_FLAGS ${FCPP_FLAGS}
