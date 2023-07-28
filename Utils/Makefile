all: gmsh2exo viewSets viewDAG YAMLValidator HookeLaws 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
include ${MEF90_DIR}/m_HeatXfer/Makefile.include
include ${MEF90_DIR}/m_DefMech/Makefile.include

### Tools
${MEF90_DIR}/${PETSC_ARCH}/bin/viewDAG: viewDAG.o
	-@${CLINKER} -o viewDAG viewDAG.o ${PETSC_LIB}
	-@${RM} -f viewDAG.o

viewDAG:${MEF90_DIR}/${PETSC_ARCH}/bin/viewDAG
	-@echo "      $@"

${MEF90_DIR}/${PETSC_ARCH}/bin/viewSets: viewSets.o
	-@${CLINKER} -o viewSets viewSets.o ${PETSC_LIB}
	-@${RM} -f viewSets.o

viewSets:${MEF90_DIR}/${PETSC_ARCH}/bin/viewSets
	-@echo "      $@"

gmsh2exo:${MEF90_DIR}/${PETSC_ARCH}/bin/gmsh2exo
	-@echo "      $@"

gmsh2exo.o: ${MEF90_DIR}/Utils/gmsh2exo.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

${MEF90_DIR}/${PETSC_ARCH}/bin/gmsh2exo: gmsh2exo.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestOrientation:${MEF90_DIR}/${PETSC_ARCH}/bin/TestOrientation
	-@echo "      $@"

TestOrientation.o: ${MEF90_DIR}/Utils/TestOrientation.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

${MEF90_DIR}/${PETSC_ARCH}/bin/TestOrientation: TestOrientation.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

YAMLValidator.o: ${MEF90_DIR}/Utils/YAMLValidator.F90 
	${FC} ${FC_FLAGS} ${FCPPFLAGS} -c -o $@ $<

YAMLValidator:${MEF90_DIR}/${PETSC_ARCH}/bin/YAMLValidator
	-@echo "      $@"

${MEF90_DIR}/${PETSC_ARCH}/bin/YAMLValidator: YAMLValidator.o
	-@echo "      $@"
	-@${FLINKER} -o $@ YAMLValidator.o ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

HookeLaws:${MEF90_DIR}/${PETSC_ARCH}/bin/HookeLaws
	-@echo "      $@"

HookeLaws.o: ${MEF90_DIR}/Utils/HookeLaws.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

${MEF90_DIR}/${PETSC_ARCH}/bin/HookeLaws: HookeLaws.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

MEF90:
	make -C ${MEF90_DIR} MEF90

clean::
	-@rm -f gmsh2exo viewSets viewDAG YAMLValidator HookeLaws 
allclean:
	make -C ${MEF90_DIR} clean

debug::
	@echo FFLAGS              ${FFLAGS} 
	@echo FC_FLAGS            ${FC_FLAGS}
	@echo FCFLAGS             ${FCFLAGS}
	@echo FCPPFLAGS           ${FCPPFLAGS}
	@echo FCPP_FLAGS          ${FCPP_FLAGS}
	@echo LDFLAGS             ${LDFLAGS}
	@echo CLINKER             ${CLINKER}
	@echo FLINKER             ${FLINKER}
	@echo PETSC_LIB           ${PETSC_LIB}
	@echo MEF90_HEATXFER_OBJS ${MEF90_HEATXFER_OBJS}
	@echo MEF90_INCLUDE       ${MEF90_INCLUDE}
	@echo DIFF                ${DIFF}
