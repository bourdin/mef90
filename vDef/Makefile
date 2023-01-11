NP=2 
all:
	-@make -C ${MEF90_DIR} vDef
test: test2D test3D
test2D: ${MEF90_DIR}/${PETSC_ARCH}/bin/vDef
	-@for OPTIONSFILE in test7.yaml test8.yaml; do \
	make -s testvDef PREFIX=Beam1x10-tri OPTIONSFILE=$${OPTIONSFILE};\
	done

test3D: ${MEF90_DIR}/${PETSC_ARCH}/bin/vDef
	-@for OPTIONSFILE in test7.yaml test8.yaml; do \
	make -s testvDef PREFIX=Beam1x10-tetra4  OPTIONSFILE=$${OPTIONSFILE};\
	done
	#make -s testvDef PREFIX=Beam1x10-tetra10 OPTIONSFILE=$${OPTIONSFILE};\

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
include ${MEF90_DIR}/m_HeatXfer/Makefile.include
include ${MEF90_DIR}/m_DefMech/Makefile.include
BASE_DIR=${MEF90_DIR}/vDef

m_vDefDefault.o: ${BASE_DIR}/m_vDefDefault.F90 ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

vDef:${MEF90_DIR}/${PETSC_ARCH}/bin/vDef
	-@echo "      $@"

${MEF90_DIR}/${PETSC_ARCH}/bin/vDef: vDef.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $@"
	-@${FLINKER} -o $@ vDef.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} 

vDef.o: ${BASE_DIR}/vDef.F90 m_vDefDefault.o ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

vDefTAO:${MEF90_DIR}/${PETSC_ARCH}/bin/vDefTAO
	-@echo "      $@"

${MEF90_DIR}/${PETSC_ARCH}/bin/vDefTAO: vDefTAO.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $@"
	-@${FLINKER} -o $@ vDefTAO.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} 

vDefTAO.o: ${BASE_DIR}/vDefTAO.F90 m_vDefDefault.o ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

vDefP:${MEF90_DIR}/${PETSC_ARCH}/bin/vDefP
	-@echo "      $@"

${MEF90_DIR}/${PETSC_ARCH}/bin/vDefP: vDefP.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $@"
	-@${FLINKER} -o $@ vDefP.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} 

vDefP.o: ${BASE_DIR}/vDefP.F90 m_vDefDefault.o ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

vDefHF:${MEF90_DIR}/${PETSC_ARCH}/bin/vDefHF
	-@echo "      $@"

${MEF90_DIR}/${PETSC_ARCH}/bin/vDefHF: vDefHF.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $@"
	-@${FLINKER} -o $@ vDefHF.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} 

vDefHF.o: ${BASE_DIR}/vDefHF.F90 m_vDefDefault.o ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

vDefUpa:${MEF90_DIR}/${PETSC_ARCH}/bin/vDefUpa
	-@echo "      $@"

${MEF90_DIR}/${PETSC_ARCH}/bin/vDefUpa: vDefUpa.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $@"
	-@${FLINKER} -o $@ vDefUpa.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} 

vDefUpa.o: ${BASE_DIR}/vDefUpa.F90 m_vDefDefault.o ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

vDefBT:${MEF90_DIR}/${PETSC_ARCH}/bin/vDefBT
	-@echo "      $@"

${MEF90_DIR}/${PETSC_ARCH}/bin/vDefBT: vDefBT.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $@"
	-@${FLINKER} -o $@ vDefBT.o m_vDefDefault.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} 

vDefBT.o: ${BASE_DIR}/vDefBT.F90 m_vDefDefault.o ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

testvDef:
	@echo ;echo "Testing vDef with $${PREFIX} $${OPTIONSFILE}"
	@ARGS="-prefix data/$${PREFIX} -options_file_yaml data/$${OPTIONSFILE}";\
	${MPIEXEC} -n ${NP} ${MEF90_DIR}/${PETSC_ARCH}/bin/vDef $${ARGS} | grep -v -e ' \#' > $${PREFIX}-$${OPTIONSFILE}.out 2>&1;\
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
