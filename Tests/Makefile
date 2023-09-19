all: TestDotProduct TestMassMatrix TestHeatXferCtx TestMEF90Ctx TestSection \
     TestConstraints TestConstraintIO TestVec TestDMPlexVecGetClosure \
	 TestVecSetBCFromOptions TestLabels TestATClass TestSNLP\
	TestHookeLaws \
	 viewDAG YAMLValidator HookeLaws 
runall: runTestDotProduct runTestMassMatrix runTestHeatXferCtx runTestMEF90Ctx \
        runTestSection runTestConstraints runTestConstraintIO runTestVec runTestDMPlexVecGetClosure \
		runTestVecSetBCFromOptions runTestLabels runTestATClass runTestSNLP

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
include ${MEF90_DIR}/m_HeatXfer/Makefile.include
include ${MEF90_DIR}/m_DefMech/Makefile.include

TestOrientation.o: ${MEF90_DIR}/Tests/TestOrientation.F90 MEF90 m_HeatXfer
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestOrientation: TestOrientation.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestDotProduct.o: ${MEF90_DIR}/Tests/TestDotProduct.F90 MEF90 m_HeatXfer
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestDotProduct: TestDotProduct.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestMassMatrix.o: ${MEF90_DIR}/Tests/TestMassMatrix.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestMassMatrix: TestMassMatrix.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestHeatXferCtx.o: ${MEF90_DIR}/Tests/TestHeatXferCtx.F90 MEF90 m_HeatXfer
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestHeatXferCtx: TestHeatXferCtx.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestMEF90Ctx.o: ${MEF90_DIR}/Tests/TestMEF90Ctx.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestMEF90Ctx: TestMEF90Ctx.o 
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestSection.o: ${MEF90_DIR}/Tests/TestSection.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestSection: TestSection.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestConstraints.o: ${MEF90_DIR}/Tests/TestConstraints.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestConstraints: TestConstraints.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestConstraintIO.o: ${MEF90_DIR}/Tests/TestConstraintIO.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestConstraintIO: TestConstraintIO.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestConstraintIO2.o: ${MEF90_DIR}/Tests/TestConstraintIO2.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestConstraintIO2: TestConstraintIO2.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

# TestMaterials is a false positive	
# TestMaterials: TestMaterials.o
# 	-@echo " $@"
# 	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
# 	-@${RM} $<

# TestMaterials.o: ${MEF90_DIR}/Tests/TestMaterials.F90 MEF90
# 	${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

TestVec.o: ${MEF90_DIR}/Tests/TestVec.F90 MEF90
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestVec: TestVec.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestDMPlexVecGetClosure.o: ${MEF90_DIR}/Tests/TestDMPlexVecGetClosure.F90 MEF90
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestDMPlexVecGetClosure: TestDMPlexVecGetClosure.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestVecSetBCFromOptions.o: ${MEF90_DIR}/Tests/TestVecSetBCFromOptions.F90 MEF90
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestVecSetBCFromOptions: TestVecSetBCFromOptions.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestLabels.o: ${MEF90_DIR}/Tests/TestLabels.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestLabels: TestLabels.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestATClass.o: ${MEF90_DIR}/Tests/TestATClass.F90 MEF90 m_DefMech
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestATClass: TestATClass.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestSNLP.o: ${MEF90_DIR}/Tests/TestSNLP.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestSNLP: TestSNLP.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestVonMises2D.o: ${MEF90_DIR}/Tests/TestVonMises2D.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestVonMises2D: TestVonMises2D.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

TestVonMises3D.o: ${MEF90_DIR}/Tests/TestVonMises3D.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestVonMises3D: TestVonMises3D.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${SNLP_LIB} ${SNLP_FORTRAN_LIB}
	-@${RM} $<

### Old stuff
# TestDofOrdering was used to figure out petsc dof ordering. It is not needed on a regular basis
TestDofOrdering.o: ${MEF90_DIR}/Tests/TestDofOrdering.F90 MEF90 m_HeatXfer
	-@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestDofOrdering: TestDofOrdering.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_HEATXFER_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestDMPlexComputeCellGeometryAffineFEM: TestDMPlexComputeCellGeometryAffineFEM.o
	-@${CLINKER} -o TestDMPlexComputeCellGeometryAffineFEM TestDMPlexComputeCellGeometryAffineFEM.o ${PETSC_LIB}
	${RM} -f TestDMPlexComputeCellGeometryAffineFEM.o

TestDMPlexComputeCellGeometryAffineFEMF90: TestDMPlexComputeCellGeometryAffineFEMF90.o
	-@${FLINKER} -o TestDMPlexComputeCellGeometryAffineFEMF90 TestDMPlexComputeCellGeometryAffineFEMF90.o ${PETSC_LIB} ${PETSC_FORTRAN_LIB}
	${RM} -f TestDMPlexComputeCellGeometryAffineFEMF90.o

TestEXOReadPart.o: ${MEF90_DIR}/Tests/TestEXOReadPart.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestEXOReadPart: TestEXOReadPart.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB}
	-@${RM} $<

TestDMPlexComputeGeometryTri.o: ${MEF90_DIR}/Tests/TestDMPlexComputeGeometryTri.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestDMPlexComputeGeometryTri: TestDMPlexComputeGeometryTri.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB}
	-@${RM} $<

TestDMPlexComputeGeometryQuad.o: ${MEF90_DIR}/Tests/TestDMPlexComputeGeometryQuad.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestDMPlexComputeGeometryQuad: TestDMPlexComputeGeometryQuad.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB}
	-@${RM} $<

TestQuadrature.o: ${MEF90_DIR}/Tests/TestQuadrature.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestQuadrature: TestQuadrature.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestSpectralDecomposition.o: ${MEF90_DIR}/Tests/TestSpectralDecomposition.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestSpectralDecomposition: TestSpectralDecomposition.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestMasonry.o: ${MEF90_DIR}/Tests/TestMasonry.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

TestMasonry: TestMasonry.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestHD.o: ${MEF90_DIR}/Tests/TestHD.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

TestHD: TestHD.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestTresca3D.o: ${MEF90_DIR}/Tests/TestTresca3D.F90 MEF90
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestTresca3D: TestTresca3D.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

TestTresca2D.o: ${MEF90_DIR}/Tests/TestTresca2D.F90 MEF90
	${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

TestTresca2D: TestTresca2D.o
	-@echo " $@"
	-@${FLINKER} -o $@ $< ${MEF90_OBJS} ${MEF90_DEFMECH_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	-@${RM} $<

###
### TESTS
###

runTestDotProduct: TestDotProduct
	@echo "         Running " $@
	-@ARGS='-geometry ${MEF90_DIR}/TestMeshes/Beam1x10-tri.gen -result $@.exo -verbose 1';\
	${MPIEXEC} -n 2 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestMassMatrix: TestMassMatrix
	@echo "         Running " $@
	-@ARGS='-geometry ${MEF90_DIR}/TestMeshes/Beam1x10-tri.gen -result $@.exo -verbose 1';\
	${MPIEXEC} -n 2 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestHeatXferCtx: TestHeatXferCtx
	@echo "         Running " $@
	-@ARGS='-geometry ${MEF90_DIR}/TestMeshes/SquareFaceSet.gen -result $@.exo -verbose 1 -time_numstep 3';\
	${MPIEXEC} -n 2 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestMEF90Ctx: TestMEF90Ctx
	@echo "         Running " $@
	-@ARGS='-geometry ${MEF90_DIR}/TestMeshes/SquareFaceSet.gen -result $@.exo -verbose 1';\
	${MPIEXEC} -n 2 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestSection: TestSection
	@echo "         Running " $@
	-@ARGS='-geometry ${MEF90_DIR}/TestMeshes/SquareFaceSet.gen -result $@.exo -mef90dm_view -sectionU_view -sectionU0_view -sectionV_view -sectionV0_view';\
	${MPIEXEC} -n 2 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestConstraints: TestConstraints
	@echo "         Running " $@
	-@ARGS='-geometry ${MEF90_DIR}/TestMeshes/SquareFaceSet.gen -result $@.exo -mef90section_view';\
	${MPIEXEC} -n 2 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestConstraintIO: TestConstraintIO
	@echo "         Running " $@
	-@ARGS='-geometry ${MEF90_DIR}/TestMeshes/SquareFaceSet.gen -result $@.exo -sectionU_view -sectionU0_view -sectionSigma_view -sectionSigma0_view';\
	${MPIEXEC} -n 2 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

# TestMaterials is a false positive	
# runTestMaterials: TestMaterials
# 	@echo "         Running " $@
# 	-@ARGS='';\
# 	${MPIEXEC} -n 2 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
# 	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
# 	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
# 	${RM} -f $@.out

runTestVec: TestVec
	@echo "         Running " $@
	-@ARGS='-geometry ${MEF90_DIR}/TestMeshes/SquareFaceSet.gen -result $@.exo -mef90section_view -mef90vec_view';\
	${MPIEXEC} -n 2 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestDMPlexVecGetClosure: TestDMPlexVecGetClosure
	@echo "         Running " $@
	-@ARGS='-geometry ${MEF90_DIR}/TestMeshes/SquareFaceSet.gen -result $@.exo -mef90section_view -mef90vec_view';\
	${MPIEXEC} -n 2 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestVecSetBCFromOptions: TestVecSetBCFromOptions
	@echo "         Running " $@
	-@ARGS='-geometry ${MEF90_DIR}/TestMeshes/SquareFaceSet.gen -result $@.exo -temperature_section_view    \
     -temperature_vec_view -sdim 2 -fs0021_temperatureBC yes,no -fs0021_boundaryTemperature 10,20 \
     -vs0010_temperatureBC yes,yes -vs0010_boundaryTemperature 100,101 -cs0001_temperature 1000   \
     -vs0010_temperature -2 -fs0020_temperature 20 -fs0021_temperature 21 -vs0010_temperature 100 ';\
	${MPIEXEC} -n 2 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestLabels: TestLabels
	@echo "         Running " $@
	-@ARGS='-geometry ${MEF90_DIR}/TestMeshes/SquareFaceSet.gen -result $@.exo -mef90section_view -mef90vec_view';\
	${MPIEXEC} -n 2 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestATClass: TestATClass
	@echo "         Running " $@
	-@ARGS='';\
	${MPIEXEC} -n 1 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestSNLP: TestSNLP
	@echo "         Running " $@
	-@ARGS='';\
	${MPIEXEC} -n 1 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestVonMises2D: TestVonMises2D
	@echo "         Running " $@
	-@ARGS='';\
	${MPIEXEC} -n 1 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

runTestVonMises3D: TestVonMises3D
	@echo "         Running " $@
	-@ARGS='';\
	${MPIEXEC} -n 1 ./$< $${ARGS} | grep -v -e ' \#' > $@.out 2>&1;\
	if (${DIFF} -B ${MEF90_DIR}/Tests/results/$@.out $@.out) then echo "         ... passed"; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above "; fi;\
	${RM} -f $@.out $@.exo

MEF90:
	make -C ${MEF90_DIR} MEF90
m_HeatXfer:
	make -C ${MEF90_DIR} m_HeatXfer
m_DefMech:
	make -C ${MEF90_DIR} m_DefMech

clean::
	-@rm -f TestMaterials TestMEF90Ctx testHeatXfer testLeak
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
