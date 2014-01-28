all: TestMaterials TestMEF90Ctx TestHeatXfer
runall: runTestMaterials runTestMEF90Ctx runTestHeatXfer

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
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestMaterials: TestMaterials.o
	-@echo "      $@"
	${FLINKER} -o $@ $< ${MEF90_OBJS}  ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestMaterials.o: ${MEF90_DIR}/Tests/TestMaterials.F90 MEF90
	${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

TestMEF90Ctx: TestMEF90Ctx.o 
	-@echo "      $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS}  ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestMEF90Ctx.o: ${MEF90_DIR}/Tests/TestMEF90Ctx.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestUpdateSets.o: ${MEF90_DIR}/Tests/TestUpdateSets.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestUpdateSets: TestUpdateSets.o
	-@echo "      $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS}  ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestNSP.o: ${MEF90_DIR}/Tests/TestNSP.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestNSP: TestNSP.o
	-@echo "      $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS}  ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestLeak.o: ${MEF90_DIR}/Tests/TestLeak.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestLeak: TestLeak.o
	-@echo "      $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS}  ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

HookeLaws.o: ${MEF90_DIR}/Tests/HookeLaws.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

HookeLaws: HookeLaws.o
	-@echo "      $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS}  ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
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
	
runTestHeatXfer: TestHeatXfer
	@echo "         Running " $@
	-@cp ${MEF90_DIR}/TestMeshes/SquareNG-tri3.gen .
	-@ARGS='-verbose 1 -prefix SquareNG-tri3 -time_numstep 51';\
	${MPIEXEC} -n 2 ./$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out SquareNG-tri3*
	
runSimplePoissonNG2D-1: 
	make -C ${MEF90_DIR} SimplePoissonNG2D
	@echo "         Running " $@
	-@cp ${MEF90_DIR}/TestMeshes/SquareNG-tri3.gen .
	-@ARGS='-verbose 1 -prefix SquareNG-tri3 -options_file ../SimplePoissonNG/opt1.txt';\
	${MPIEXEC} -n 2 ../bin/${PETSC_ARCH}/SimplePoissonNG2D $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out SquareNG-tri3*

runSimplePoissonNG2D-1TS: 
	make -C ${MEF90_DIR} SimplePoissonNG2D
	@echo "         Running " $@
	-@cp ${MEF90_DIR}/TestMeshes/SquareNG-tri3.gen .
	-@ARGS='-verbose 1 -prefix SquareNG-tri3 -options_file ../SimplePoissonNG/opt1TS.txt';\
	${MPIEXEC} -n 2 ../bin/${PETSC_ARCH}/SimplePoissonNG2D $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out SquareNG-tri3*

runSimplePoissonNG2D-2: 
	make -C ${MEF90_DIR} SimplePoissonNG2D
	@echo "         Running " $@
	-@cp ${MEF90_DIR}/TestMeshes/SquareNG-tri3.gen .
	-@ARGS='-verbose 1 -prefix SquareNG-tri3 -options_file ../SimplePoissonNG/opt2.txt';\
	${MPIEXEC} -n 2 ../bin/${PETSC_ARCH}/SimplePoissonNG2D $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out SquareNG-tri3*

runSimplePoissonNG2D-3: 
	make -C ${MEF90_DIR} SimplePoissonNG2D
	@echo "         Running " $@
	-@cp ${MEF90_DIR}/TestMeshes/SquareNG-tri3.gen .
	-@ARGS='-verbose 1 -prefix SquareNG-tri3 -options_file ../SimplePoissonNG/opt3.txt';\
	${MPIEXEC} -n 2 ../bin/${PETSC_ARCH}/SimplePoissonNG2D $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out SquareNG-tri3*


MEF90:
	make -C ${MEF90_DIR} MEF90
m_HeatXfer:
	make -C ${MEF90_DIR} m_HeatXfer
	
clean::
	-@rm TestMaterials TestMEF90Ctx testHeatXfer
allclean:
	make -C ${MEF90_DIR} clean
	
debug::
	@echo FFLAGS ${FFLAGS} 
	@echo FC_FLAGS ${FC_FLAGS}
	@echo FCFLAGS ${FCFLAGS}
	@echo FCPPFLAGS ${FCPPFLAGS}
	@echo FCPP_FLAGS ${FCPP_FLAGS}
