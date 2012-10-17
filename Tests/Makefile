all: TestMaterials TestMEF90Ctx
runall: runTestMaterials runTestMEF90Ctx

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include

TestMaterials: TestMaterials.o
	-@echo "      $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS}  ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestMaterials.o: ${MEF90_DIR}/Tests/TestMaterials.F90 MEF90
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestMEF90Ctx: TestMEF90Ctx.o 
	-@echo "      $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS}  ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestMEF90Ctx.o: ${MEF90_DIR}/Tests/TestMEF90Ctx.F90 MEF90
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<


runTestMaterials: TestMaterials
	@echo "Running " $@
	-@ARGS='';\
	${MPIEXEC} -n 2 $< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out
	-@#echo ${MPIEXEC} -n 2 $< $${ARGS};\

runTestMEF90Ctx: TestMEF90Ctx
	@echo "Running " $@
	-@ARGS='-verbose 1';\
	${MPIEXEC} -n 2 $< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out
	-@#echo ${MPIEXEC} -n 2 $< $${ARGS};\

MEF90:
	make -C ${MEF90_DIR} MEF90

clean::
	-@rm TestMaterials TestMEF90Ctx