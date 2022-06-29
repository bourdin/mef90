all: TestMaterials TestMEF90Ctx 
runall: runTestMaterials runTestMEF90Ctx 

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
#include ${MEF90_DIR}/m_HeatXfer/Makefile.include
#include ${MEF90_DIR}/m_DefMech/Makefile.include

TestConstraintIO.o: ${MEF90_DIR}/Tests/TestConstraintIO.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestConstraintIO: TestConstraintIO.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<
	
viewDAG: viewDAG.o
	-@${CLINKER} -o viewDAG viewDAG.o ${PETSC_LIB}
	${RM} -f viewDAG.o

TestLabels.o: ${MEF90_DIR}/Tests/TestLabels.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestLabels: TestLabels.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestSection.o: ${MEF90_DIR}/Tests/TestSection.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestSection: TestSection.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestConstraints.o: ${MEF90_DIR}/Tests/TestConstraints.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestConstraints: TestConstraints.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestDMPlexComputeGeometryTri.o: ${MEF90_DIR}/Tests/TestDMPlexComputeGeometryTri.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestDMPlexComputeGeometryTri: TestDMPlexComputeGeometryTri.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestDMPlexComputeGeometryQuad.o: ${MEF90_DIR}/Tests/TestDMPlexComputeGeometryQuad.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestDMPlexComputeGeometryQuad: TestDMPlexComputeGeometryQuad.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestLinAlg.o: ${MEF90_DIR}/Tests/TestLinAlg.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestLinAlg: TestLinAlg.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

YAMLValidator.o: ${MEF90_DIR}/Tests/YAMLValidator.F90 
	${FC} ${FC_FLAGS} ${FCPPFLAGS} -c -o $@ $<

#YAMLValidator: YAMLValidator.o
#	-@echo " $@"
#	-@${FLINKER} -o $@ $< ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
#	-@${RM} $<
YAMLValidator:${MEF90_DIR}/bin/${PETSC_ARCH}/YAMLValidator
	echo "      $@"

${MEF90_DIR}/bin/${PETSC_ARCH}/YAMLValidator: YAMLValidator.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS}  ${MEF90_HEATXFER_OBJS}
	-@echo "      $@"
	${FLINKER} -o $@ YAMLValidator.o ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} 

TestMaterials: TestMaterials.o
	-@echo " $@"
	${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestMaterials.o: ${MEF90_DIR}/Tests/TestMaterials.F90 MEF90
	${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

TestMEF90Ctx: TestMEF90Ctx.o 
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestMEF90Ctx.o: ${MEF90_DIR}/Tests/TestMEF90Ctx.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestUpdateSets.o: ${MEF90_DIR}/Tests/TestUpdateSets.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestUpdateSets: TestUpdateSets.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestNSP.o: ${MEF90_DIR}/Tests/TestNSP.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestNSP: TestNSP.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestLeak.o: ${MEF90_DIR}/Tests/TestLeak.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestLeak: TestLeak.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

HookeLaws.o: ${MEF90_DIR}/Tests/HookeLaws.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

HookeLaws: HookeLaws.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestHookeLaws.o: ${MEF90_DIR}/Tests/TestHookeLaws.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestHookeLaws: TestHookeLaws.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestQuadrature.o: ${MEF90_DIR}/Tests/TestQuadrature.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestQuadrature: TestQuadrature.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestSpectralDecomposition.o: ${MEF90_DIR}/Tests/TestSpectralDecomposition.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestSpectralDecomposition: TestSpectralDecomposition.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestATClass.o: ${MEF90_DIR}/Tests/TestATClass.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestATClass: TestATClass.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestMasonry.o: ${MEF90_DIR}/Tests/TestMasonry.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

TestMasonry: TestMasonry.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestHD.o: ${MEF90_DIR}/Tests/TestHD.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

TestHD: TestHD.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestSNLP.o: ${MEF90_DIR}/Tests/TestSNLP.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestSNLP: TestSNLP.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestVonMises3D.o: ${MEF90_DIR}/Tests/TestVonMises3D.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestVonMises3D: TestVonMises3D.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestTresca3D.o: ${MEF90_DIR}/Tests/TestTresca3D.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestTresca3D: TestTresca3D.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestTresca2D.o: ${MEF90_DIR}/Tests/TestTresca2D.F90 MEF90
	${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestTresca2D: TestTresca2D.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestVonMises2D.o: ${MEF90_DIR}/Tests/TestVonMises2D.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestVonMises2D: TestVonMises2D.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestInelasticStrainCell.o: ${MEF90_DIR}/Tests/TestInelasticStrainCell.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $< 

TestInelasticStrainCell: TestInelasticStrainCell.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB} ${MEF90_DEFMECH_OBJS} ${MEF90_HEATXFER_OBJS}
	-@${RM} $<

ExodusBuild: ExodusBuild.o  chkopts
	${CLINKER} -o ExodusBuild ExodusBuild.o ${PETSC_LIB}
	${RM} ExodusBuild.o


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

MEF90:
	make -C ${MEF90_DIR} MEF90
m_HeatXfer:
	make -C ${MEF90_DIR} m_HeatXfer

clean::
	-@rm -f TestMaterials TestMEF90Ctx testHeatXfer testLeak
allclean:
	make -C ${MEF90_DIR} clean

debug::
	@echo FFLAGS     ${FFLAGS} 
	@echo FC_FLAGS   ${FC_FLAGS}
	@echo FCFLAGS    ${FCFLAGS}
	@echo FCPPFLAGS  ${FCPPFLAGS}
	@echo FCPP_FLAGS ${FCPP_FLAGS}
	@echo LDFLAGS    ${LDFLAGS}
	@echo CLINKER    ${CLINKER}
	@echo FLINKER    ${FLINKER}
	@echo PETSC_LIB  ${PETSC_LIB}
