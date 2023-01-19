NP=4 
all:
	-@make -C ${MEF90_DIR} ThermoElasticity
test: test2D test3D
test2D: ${MEF90_DIR}/${PETSC_ARCH}/bin/ThermoElasticity
	-@for i in {1..4}; do \
	make -s testThermoElasticity PREFIX=Beam1x10-tri.msh ORDER=1 OPTIONSFILE=testElastDisp$${i}.yaml;\
	make -s testThermoElasticity PREFIX=Beam1x10-tri.msh ORDER=2 OPTIONSFILE=testElastDisp$${i}.yaml;\
	done

test3D: ${MEF90_DIR}/${PETSC_ARCH}/bin/ThermoElasticity
	-@for i in {1..4}; do \
	make -s testThermoElasticity PREFIX=Beam3D1x10-tet.msh ORDER=1 OPTIONSFILE=testElastDisp$${i}.yaml;\
	make -s testThermoElasticity PREFIX=Beam3D1x10-tet.msh ORDER=2 OPTIONSFILE=testElastDisp$${i}.yaml;\
	done

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
include ${MEF90_DIR}/m_HeatXfer/Makefile.include
include ${MEF90_DIR}/m_DefMech/Makefile.include
BASE_DIR=${MEF90_DIR}/ThermoElasticity

ThermoElasticity:${MEF90_DIR}/${PETSC_ARCH}/bin/ThermoElasticity ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $@"

${MEF90_DIR}/${PETSC_ARCH}/bin/ThermoElasticity: ThermoElasticity.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} 
	-@echo "      $@"
	-@${FLINKER} -o $@ ThermoElasticity.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}

ThermoElasticity.o: ${BASE_DIR}/ThermoElasticity.F90 ${MEF90_OBJS}  ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} m_vDefDefault.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_vDefDefault.o: ${BASE_DIR}/m_vDefDefault.F90 ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<



testThermoElasticity:
	@echo ;echo "Testing ThermoElasticity with $${PREFIX} $${OPTIONSFILE} order $${ORDER}"
	# @ARGS="-geometry ${MEF90_DIR}/TestMeshes/$${PREFIX} -result $${OPTIONSFILE}.exo -element_order $${ORDER} -options_file data/$${OPTIONSFILE}";\
	# echo ${MPIEXEC} -n ${NP} ${MEF90_DIR}/${PETSC_ARCH}/bin/ThermoElasticity $${ARGS}
	ARGS="-geometry ${MEF90_DIR}/TestMeshes/$${PREFIX} -result $${PREFIX}-$${OPTIONSFILE}-$${ORDER}.exo -element_order $${ORDER} -options_file data/$${OPTIONSFILE}";\
	${MPIEXEC} -n ${NP} ${MEF90_DIR}/${PETSC_ARCH}/bin/ThermoElasticity $${ARGS} | grep -v -e ' \#' > $${PREFIX}-$${OPTIONSFILE}-$${ORDER}.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/ThermoElasticity/results/$${PREFIX}-$${OPTIONSFILE}-$${ORDER}.out $${PREFIX}-$${OPTIONSFILE}-$${ORDER}.out) then echo "Passed";\
	else echo "Possible problem. Diffs above"; fi
	-@${RM} $${PREFIX}-$${OPTIONSFILE}-$${ORDER}* 
	
MEF90:
	make -C ${MEF90_DIR} MEF90
m_DefMech:
	make -C ${MEF90_DIR} m_DefMech
	
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
