NP=2
all:
	-@make -C ${MEF90_DIR} HeatXfer
test: test2D test3D testTransient2D testTransient3D
test2D: ${MEF90_DIR}/bin/${PETSC_ARCH}/HeatXfer
	-@for OPTIONSFILE in test1.opts test2.opts test3.opts; do \
	make testHeatXfer PREFIX=SquareNG_X-tri3 OPTIONSFILE=$${OPTIONSFILE};\
	make testHeatXfer PREFIX=SquareNG_X-tri6 OPTIONSFILE=$${OPTIONSFILE};\
	make testHeatXfer PREFIX=SquareNG_Y-tri3 OPTIONSFILE=$${OPTIONSFILE};\
	make testHeatXfer PREFIX=SquareNG_Y-tri6 OPTIONSFILE=$${OPTIONSFILE};\
	done

test3D: ${MEF90_DIR}/bin/${PETSC_ARCH}/HeatXfer
	-@for OPTIONSFILE in test1.opts test2.opts test3.opts; do \
	make testHeatXfer PREFIX=CubeNG_X-tetra4  OPTIONSFILE=$${OPTIONSFILE};\
	make testHeatXfer PREFIX=CubeNG_X-tetra10 OPTIONSFILE=$${OPTIONSFILE};\
	make testHeatXfer PREFIX=CubeNG_Y-tetra4  OPTIONSFILE=$${OPTIONSFILE};\
	make testHeatXfer PREFIX=CubeNG_Y-tetra10 OPTIONSFILE=$${OPTIONSFILE};\
	make testHeatXfer PREFIX=CubeNG_Z-tetra4  OPTIONSFILE=$${OPTIONSFILE};\
	make testHeatXfer PREFIX=CubeNG_Z-tetra10 OPTIONSFILE=$${OPTIONSFILE};\
	done

testTransient2D: ${MEF90_DIR}/bin/${PETSC_ARCH}/HeatXfer
	-@for OPTIONSFILE in test5.opts test6.opts; do \
	make testHeatXfer PREFIX=Beam1x10-tri3 OPTIONSFILE=$${OPTIONSFILE};\
	make testHeatXfer PREFIX=Beam1x10-tri6 OPTIONSFILE=$${OPTIONSFILE};\
	done

testTransient3D: ${MEF90_DIR}/bin/${PETSC_ARCH}/HeatXfer
	-@for OPTIONSFILE in test5.opts test6.opts; do \
	make testHeatXfer PREFIX=Beam1x10-tetra4  OPTIONSFILE=$${OPTIONSFILE};\
	make testHeatXfer PREFIX=Beam1x10-tetra10 OPTIONSFILE=$${OPTIONSFILE};\
	done

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

testHeatXfer:
	@echo "Testing HeatXfer with $${PREFIX} $${OPTIONSFILE}"
	@ARGS="-prefix data/$${PREFIX} -options_file data/$${OPTIONSFILE}";\
	${MPIEXEC} -n ${NP} ${MEF90_DIR}/bin/${PETSC_ARCH}/HeatXfer $${ARGS}| grep -v -e ' \#' > $${PREFIX}-$${OPTIONSFILE}.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/HeatXfer/results/${PREFIX}-$${OPTIONSFILE}.out ${PREFIX}-$${OPTIONSFILE}.out) then echo "         ... passed";\
	else echo "Possible problem. Diffs above \n========================================="; fi
	-@${RM} $${PREFIX}-$${OPTIONSFILE}.out 
	
MEF90:
	make -C ${MEF90_DIR} MEF90
m_HeatXfer:
	make -C ${MEF90_DIR} m_HeatXfer
	
clean::
	-@rm HeatXfer
allclean:
	make -C ${MEF90_DIR} clean
	
debug::
	@echo FFLAGS ${FFLAGS} 
	@echo FC_FLAGS ${FC_FLAGS}
	@echo FCFLAGS ${FCFLAGS}
	@echo FCPPFLAGS ${FCPPFLAGS}
	@echo FCPP_FLAGS ${FCPP_FLAGS}
