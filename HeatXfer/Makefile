all:
	-@make -C ${MEF90_DIR} HeatXfer
test: runTestHeatXfer2D1 runTestHeatXfer2D2 runTestHeatXfer2D3 runTestHeatXfer2D4 \
         runTestHeatXfer3D1 runTestHeatXfer3D2 runTestHeatXfer3D3 runTestHeatXfer3D4

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
include ${MEF90_DIR}/m_HeatXfer/Makefile.include
BASE_DIR=${MEF90_DIR}/HeatXfer

HeatXfer:${MEF90_DIR}/bin/${PETSC_ARCH}/HeatXfer
	-@echo "      $@"

${MEF90_DIR}/bin/${PETSC_ARCH}/HeatXfer: HeatXfer.o
	-@echo "      $@"
	-@${FLINKER} -o $@ HeatXfer.o ${MEF90_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

HeatXfer.o: ${BASE_DIR}/HeatXfer.F90 ${MEF90_HEATXFER_OBJS}
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

runTestHeatXfer2D1: HeatXfer
	@echo "         Running " $@
	-@ARGS='-prefix data/SquareNG_X-tri3 -options_file data/test1.opts';\
	${MPIEXEC} -n 2 ${MEF90_DIR}/bin/${PETSC_ARCH}/$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/HeatXfer/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} $@.out

runTestHeatXfer2D2: HeatXfer
	@echo "         Running " $@
	-@ARGS='-prefix data/SquareNG_Y-tri3 -options_file data/test1.opts';\
	${MPIEXEC} -n 2 ${MEF90_DIR}/bin/${PETSC_ARCH}/$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/HeatXfer/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} $@.out
	
runTestHeatXfer2D3: HeatXfer
	@echo "         Running " $@
	-@ARGS='-prefix data/SquareNG_X-tri6 -options_file data/test1.opts';\
	${MPIEXEC} -n 2 ${MEF90_DIR}/bin/${PETSC_ARCH}/$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/HeatXfer/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} $@.out
	
runTestHeatXfer2D4: HeatXfer
	@echo "         Running " $@
	-@ARGS='-prefix data/SquareNG_Y-tri6 -options_file data/test1.opts';\
	${MPIEXEC} -n 2 ${MEF90_DIR}/bin/${PETSC_ARCH}/$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/HeatXfer/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} $@.out
	
runTestHeatXfer3D1: HeatXfer
	@echo "         Running " $@
	-@ARGS='-prefix data/CubeNG_X-tetra4 -options_file data/test1.opts';\
	${MPIEXEC} -n 2 ${MEF90_DIR}/bin/${PETSC_ARCH}/$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/HeatXfer/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} $@.out

runTestHeatXfer3D2: HeatXfer
	@echo "         Running " $@
	-@ARGS='-prefix data/CubeNG_X-tetra4 -options_file data/test1.opts';\
	${MPIEXEC} -n 2 ${MEF90_DIR}/bin/${PETSC_ARCH}/$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/HeatXfer/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} $@.out
	
runTestHeatXfer3D3: HeatXfer
	@echo "         Running " $@
	-@ARGS='-prefix data/CubeNG_X-tetra10 -options_file data/test1.opts';\
	${MPIEXEC} -n 2 ${MEF90_DIR}/bin/${PETSC_ARCH}/$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/HeatXfer/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} $@.out
	
runTestHeatXfer3D4: HeatXfer
	@echo "         Running " $@
	-@ARGS='-prefix data/CubeNG_Y-tetra10 -options_file data/test1.opts';\
	${MPIEXEC} -n 2 ${MEF90_DIR}/bin/${PETSC_ARCH}/$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/HeatXfer/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} $@.out
	
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
