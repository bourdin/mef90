NP=2
all:
	-@make -C ${MEF90_DIR} HeatXfer
test: test2D test3D testTransient2D testTransient3D
test2D:
	-@for OPTIONSFILE in test1.yaml test2.yaml test3.yaml; do \
	make -s testHeatXfer PREFIX=SquareNG_X-tri3 OPTIONSFILE=$${OPTIONSFILE};\
	make -s testHeatXfer PREFIX=SquareNG_X-tri6 OPTIONSFILE=$${OPTIONSFILE};\
	make -s testHeatXfer PREFIX=SquareNG_Y-tri3 OPTIONSFILE=$${OPTIONSFILE};\
	make -s testHeatXfer PREFIX=SquareNG_Y-tri6 OPTIONSFILE=$${OPTIONSFILE};\
	done

test3D:
	-@for OPTIONSFILE in test1.yaml test2.yaml test3.yaml; do \
	make -s testHeatXfer PREFIX=CubeNG_X-tetra4  OPTIONSFILE=$${OPTIONSFILE};\
	make -s testHeatXfer PREFIX=CubeNG_X-tetra10 OPTIONSFILE=$${OPTIONSFILE};\
	make -s testHeatXfer PREFIX=CubeNG_Y-tetra4  OPTIONSFILE=$${OPTIONSFILE};\
	make -s testHeatXfer PREFIX=CubeNG_Y-tetra10 OPTIONSFILE=$${OPTIONSFILE};\
	make -s testHeatXfer PREFIX=CubeNG_Z-tetra4  OPTIONSFILE=$${OPTIONSFILE};\
	make -s testHeatXfer PREFIX=CubeNG_Z-tetra10 OPTIONSFILE=$${OPTIONSFILE};\
	done

testTransient2D:
	-@for OPTIONSFILE in test5.yaml test6.yaml; do \
	make -s testHeatXfer PREFIX=Beam1x10-tri3 OPTIONSFILE=$${OPTIONSFILE};\
	make -s testHeatXfer PREFIX=Beam1x10-tri6 OPTIONSFILE=$${OPTIONSFILE};\
	done

testTransient3D:
	-@for OPTIONSFILE in test5.yaml test6.yaml; do \
	make -s testHeatXfer PREFIX=Beam1x10-tetra4  OPTIONSFILE=$${OPTIONSFILE};\
	make -s testHeatXfer PREFIX=Beam1x10-tetra10 OPTIONSFILE=$${OPTIONSFILE};\
	done

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
include ${MEF90_DIR}/m_HeatXfer/Makefile.include
BASE_DIR=${MEF90_DIR}/HeatXfer

HeatXfer:${MEF90_DIR}/${PETSC_ARCH}/bin/HeatXfer
	-@echo "      $@"

${MEF90_DIR}/${PETSC_ARCH}/bin/HeatXfer: HeatXfer.o
	-@echo "      $@"
	-@${FLINKER} -o $@ HeatXfer.o ${MEF90_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

HeatXfer.o: ${BASE_DIR}/HeatXfer.F90 ${MEF90_HEATXFER_OBJS}
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

testHeatXfer:
	@echo ;echo "Testing HeatXfer with $${PREFIX} $${OPTIONSFILE}"
	@ARGS="-prefix data/$${PREFIX} -options_file_yaml data/$${OPTIONSFILE}";\
	${MPIEXEC} -n ${NP} ${MEF90_DIR}/${PETSC_ARCH}/bin/HeatXfer $${ARGS}| grep -v -e ' \#' > $${PREFIX}-$${OPTIONSFILE}.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/HeatXfer/results/${PREFIX}-$${OPTIONSFILE}.out ${PREFIX}-$${OPTIONSFILE}.out) then echo "Passed";\
	else echo "Possible problem. Diffs above"; fi
	-@${RM} $${PREFIX}-$${OPTIONSFILE}.out 
	
MEF90:
	make -C ${MEF90_DIR} MEF90
m_HeatXfer:
	make -C ${MEF90_DIR} m_HeatXfer
	
testclean::
	-@rm -f data/*out.gen data/*log

allclean: testclean
	make -C ${MEF90_DIR} clean
	
debug::
	@echo FFLAGS ${FFLAGS} 
	@echo FC_FLAGS ${FC_FLAGS}
	@echo FCFLAGS ${FCFLAGS}
	@echo FCPPFLAGS ${FCPPFLAGS}
	@echo FCPP_FLAGS ${FCPP_FLAGS}
