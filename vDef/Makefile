NP=2 
all:
	-@make -C ${MEF90_DIR} vDef
test: test2D test3D
test2D: ${MEF90_DIR}/bin/${PETSC_ARCH}/vDef
	-@for OPTIONSFILE in test7.yaml test8.yaml; do \
	make -s testvDef PREFIX=Beam1x10-tri3 OPTIONSFILE=$${OPTIONSFILE};\
	make -s testvDef PREFIX=Beam1x10-tri6 OPTIONSFILE=$${OPTIONSFILE};\
	done

test3D: ${MEF90_DIR}/bin/${PETSC_ARCH}/vDef
	-@for OPTIONSFILE in test7.yaml test8.yaml; do \
	make -s testvDef PREFIX=Beam1x10-tetra4  OPTIONSFILE=$${OPTIONSFILE};\
	done
	#make -s testvDef PREFIX=Beam1x10-tetra10 OPTIONSFILE=$${OPTIONSFILE};\

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables

include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
include ${MEF90_DIR}/m_HeatXfer/Makefile.include
include ${MEF90_DIR}/m_DefMech/Makefile.include
BASE_DIR=${MEF90_DIR}/vDef

vDef:${MEF90_DIR}/bin/${PETSC_ARCH}/vDef
	-@echo "      $@"

m_vDefDefault.o: ${BASE_DIR}/m_vDefDefault.F90 ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

${MEF90_DIR}/bin/${PETSC_ARCH}/vDef: vDef.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS}  ${MEF90_HEATXFER_OBJS}
	-@echo "      $@"
	-@${FLINKER} -o $@ vDef.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} 

vDef.o: ${BASE_DIR}/vDef.F90 m_vDefDefault.o ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

vDefP:${MEF90_DIR}/bin/${PETSC_ARCH}/vDefP
	-@echo "      $@"

vDefUpa:${MEF90_DIR}/bin/${PETSC_ARCH}/vDefUpa
	-@echo "      $@"

${MEF90_DIR}/bin/${PETSC_ARCH}/vDefP: vDefP.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS}  ${MEF90_HEATXFER_OBJS}
	-@echo "      $@"
	-@${FLINKER} -o $@ vDefP.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} 

vDefP.o: ${BASE_DIR}/vDefP.F90 m_vDefDefault.o ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

${MEF90_DIR}/bin/${PETSC_ARCH}/vDefUpa: vDefUpa.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS}  ${MEF90_HEATXFER_OBJS}
	-@echo "      $@"
	-@${FLINKER} -o $@ vDefUpa.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} 

vDefUpa.o: ${BASE_DIR}/vDefUpa.F90 m_vDefDefault.o ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

testvDef:
	@echo ;echo "Testing vDef with $${PREFIX} $${OPTIONSFILE}"
	@ARGS="-prefix data/$${PREFIX} -options_file_yaml data/$${OPTIONSFILE}";\
	${MPIEXEC} -n ${NP} ${MEF90_DIR}/bin/${PETSC_ARCH}/vDef $${ARGS} | grep -v -e ' \#' > $${PREFIX}-$${OPTIONSFILE}.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/vDef/results/${PREFIX}-$${OPTIONSFILE}.out ${PREFIX}-$${OPTIONSFILE}.out) then echo "Passed";\
	else echo "Possible problem. Diffs above"; fi
	-@${RM} $${PREFIX}-$${OPTIONSFILE}.out 
	
MEF90:
	make -C ${MEF90_DIR} MEF90
m_HeatXfer:
	make -C ${MEF90_DIR} m_HeatXfer
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
	@echo FLINKER ${FLINKER}
	@echo PETSC_LIB ${PETSC_LIB}
	@echo PETSC_FORTRAN_LIB ${PETSC_FORTRAN_LIB}
