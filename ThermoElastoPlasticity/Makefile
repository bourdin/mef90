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

ThermoElastoPlasticity: ThermoElastoPlasticity.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@${RM} $<


runTestMaterials: TestMaterials
	@echo "         Running " $@
	-@cp ${MEF90_DIR}/TestMeshes/SquareNG-tri3.gen .
	-@ARGS='';\
	${MPIEXEC} -n 2 ./$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runTestMEF90Ctx: TestMEF90Ctx
	@echo "         Running " $@
	-@cp ${MEF90_DIR}/TestMeshes/SquareNG-tri3.gen .
	-@ARGS='-verbose 1 -prefix ${MEF90_DIR}/TestMeshes/SquareNG-tri3';\
	${MPIEXEC} -n 2 ./$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out SquareNG-tri3*
	
MEF90:
	make -C ${MEF90_DIR} MEF90
m_HeatXfer:
	make -C ${MEF90_DIR} m_HeatXfer
	
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
