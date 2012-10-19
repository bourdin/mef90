all: TestMaterials TestMEF90Ctx testHeatXfer
runall: runTestMaterials runTestMEF90Ctx 

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
include ${MEF90_DIR}/m_HeatXfer/Makefile.include

TestHeatXfer: TestHeatXfer.o
	-@echo "      $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestHeatXfer.o: TestHeatXfer.F90 ${MEF90_HEATXFER_OBJS} MEF90 m_HeatXfer
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestMaterials: TestMaterials.o
	-@echo "      $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS}  ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestMaterials.o: ${MEF90_DIR}/Tests/TestMaterials.F90 MEF90
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestMEF90Ctx: TestMEF90Ctx.o 
	-@echo "      $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS}  ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestMEF90Ctx.o: ${MEF90_DIR}/Tests/TestMEF90Ctx.F90 MEF90
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<


runTestMaterials: TestMaterials
	@echo "         Running " $@
	-@ARGS='';\
	${MPIEXEC} -n 2 $< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runTestMEF90Ctx: TestMEF90Ctx
	@echo "         Running " $@
	-@ARGS='-verbose 1 -prefix ${MEF90_DIR}/TestMeshes/SquareNG-tri3';\
	${MPIEXEC} -n 2 ./$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

MEF90:
	make -C ${MEF90_DIR} MEF90
m_HeatXfer:
	make -C ${MEF90_DIR} m_HeatXfer
	
clean::
	-@rm TestMaterials TestMEF90Ctx
allclean:
	make -C ${MEF90_DIR} clean
	

