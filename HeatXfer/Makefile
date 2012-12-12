all: HeatXfer
runall: runTestMaterials runTestMEF90Ctx runTestHeatXfer

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
include ${MEF90_DIR}/m_HeatXfer/Makefile.include

HeatXfer:${MEF90_DIR}/bin/${PETSC_ARCH}/HeatXfer
	-@echo "      $@"

${MEF90_DIR}/bin/${PETSC_ARCH}/HeatXfer: HeatXfer.o
	-@echo "      $@"
	-@${FLINKER} -o $@ HeatXfer.o ${MEF90_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

HeatXfer.o: HeatXfer.F90 ${MEF90_HEATXFER_OBJS} MEF90 m_HeatXfer
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

runTestHeatXfer2D1: TestHeatXfer
	@echo "         Running " $@
	-@ARGS='-prefix data/SquareNG_X-tri3 -options_file data/test1.opts';\
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
