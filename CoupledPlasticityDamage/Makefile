NP=2 
all:
	-@make -C ${MEF90_DIR} CoupledPlasticityDamage
test: test2D test3D
test2D: ${MEF90_DIR}/bin/${PETSC_ARCH}/CoupledPlasticityDamage
	-@for OPTIONSFILE in test7.yaml test8.yaml; do \
	make -s testCoupledPlasticityDamage PREFIX=Beam1x10-tri3 OPTIONSFILE=$${OPTIONSFILE};\
	make -s testCoupledPlasticityDamage PREFIX=Beam1x10-tri6 OPTIONSFILE=$${OPTIONSFILE};\
	done

test3D: ${MEF90_DIR}/bin/${PETSC_ARCH}/CoupledPlasticityDamage
	-@for OPTIONSFILE in test7.yaml test8.yaml; do \
	make -s testCoupledPlasticityDamage PREFIX=Beam1x10-tetra4  OPTIONSFILE=$${OPTIONSFILE};\
	done
	#make -s testCoupledPlasticityDamage PREFIX=Beam1x10-tetra10 OPTIONSFILE=$${OPTIONSFILE};\

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables

include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
include ${MEF90_DIR}/m_HeatXfer/Makefile.include
include ${MEF90_DIR}/m_DefMech/Makefile.include
BASE_DIR=${MEF90_DIR}/CoupledPlasticityDamage

CoupledPlasticityDamage:${MEF90_DIR}/bin/${PETSC_ARCH}/CoupledPlasticityDamage
	-@echo "      $@"

${MEF90_DIR}/bin/${PETSC_ARCH}/CoupledPlasticityDamage: CoupledPlasticityDamage.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS}  ${MEF90_HEATXFER_OBJS}
	-@echo "      $@"
	-@${FLINKER} -o $@ CoupledPlasticityDamage.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}

CoupledPlasticityDamage.o: ${BASE_DIR}/CoupledPlasticityDamage.F90 ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

testCoupledPlasticityDamage:
	@echo ;echo "Testing CoupledPlasticityDamage with $${PREFIX} $${OPTIONSFILE}"
	@ARGS="-prefix data/$${PREFIX} -options_file_yaml data/$${OPTIONSFILE}";\
	${MPIEXEC} -n ${NP} ${MEF90_DIR}/bin/${PETSC_ARCH}/CoupledPlasticityDamage $${ARGS} | grep -v -e ' \#' > $${PREFIX}-$${OPTIONSFILE}.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/CoupledPlasticityDamage/results/${PREFIX}-$${OPTIONSFILE}.out ${PREFIX}-$${OPTIONSFILE}.out) then echo "Passed";\
	else echo "Possible problem. Diffs above"; fi
	-@${RM} $${PREFIX}-$${OPTIONSFILE}.out 
	
MEF90:
	make -C ${MEF90_DIR} MEF90
m_HeatXfer:
	make -C ${MEF90_DIR} m_HeatXfer
m_DEFMECH:
	make -C ${MEF90_DIR} m_DEFMECH
testclean::
	-@rm -f data/*out.gen data/*log data/*ener*

allclean: testclean
	make -C ${MEF90_DIR} clean
	
debug::
	@echo FFLAGS ${FFLAGS} 
	@echo FC_FLAGS ${FC_FLAGS}
	@echo FCFLAGS ${FCFLAGS}
	@echo FCPPFLAGS ${FCPPFLAGS}
	@echo FCPP_FLAGS ${FCPP_FLAGS}
