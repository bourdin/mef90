MEF90_DIM=2
all:
	-@make -C ${MEF90_DIR} m_DefMech

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
BASE_DIR=${MEF90_DIR}/m_DefMech

m_DefMech:m_MEF90_DefMech.o m_MEF90_DefMechCtx.o

m_MEF90_DefMech.o: ${BASE_DIR}/m_MEF90_DefMech.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechAT.o m_MEF90_DefMechSplit2D.o  m_MEF90_DefMechSplit3D.o m_MEF90_DefMechPlasticity2D.o  m_MEF90_DefMechPlasticity3D.o #m_MEF90_DefMechAssembly2D.o  m_MEF90_DefMechAssembly3D.o 
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_DefMechATClass.o: ${BASE_DIR}/m_MEF90_DefMechATClass.F90
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_DefMechAT1.o: ${BASE_DIR}/m_MEF90_DefMechAT1.F90 m_MEF90_DefMechATClass.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_DefMechAT1Exp.o: ${BASE_DIR}/m_MEF90_DefMechAT1Exp.F90 m_MEF90_DefMechATClass.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_DefMechAT2.o: ${BASE_DIR}/m_MEF90_DefMechAT2.F90 m_MEF90_DefMechATClass.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_DefMechATKKL.o: ${BASE_DIR}/m_MEF90_DefMechATKKL.F90 m_MEF90_DefMechATClass.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_DefMechATLinSoft.o: ${BASE_DIR}/m_MEF90_DefMechATLinSoft.F90 m_MEF90_DefMechATClass.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_DefMechAT.o: ${BASE_DIR}/m_MEF90_DefMechAT.F90 m_MEF90_DefMechAT1.o m_MEF90_DefMechAT1Exp.o m_MEF90_DefMechATKKL.o m_MEF90_DefMechAT2.o m_MEF90_DefMechATKKL.o m_MEF90_DefMechATLinSoft.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_DefMechSplitClass2D.o: ${BASE_DIR}/m_MEF90_DefMechSplitClass.F90
	@echo "      $< MEF90_DIM=2"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechSplitClass3D.o: ${BASE_DIR}/m_MEF90_DefMechSplitClass.F90
	@echo "      $< MEF90_DIM=3"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechSplit2D.o: ${BASE_DIR}/m_MEF90_DefMechSplit.F90 m_MEF90_DefMechSplitHD2D.o m_MEF90_DefMechSplitNone2D.o m_MEF90_DefMechSplitDeviatoric2D.o m_MEF90_DefMechSplitHydrostatic2D.o #m_MEF90_DefMechSplitMasonry2D.o 
	@echo "      $< MEF90_DIM=2"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechSplit3D.o: ${BASE_DIR}/m_MEF90_DefMechSplit.F90 m_MEF90_DefMechSplitHD3D.o m_MEF90_DefMechSplitNone3D.o m_MEF90_DefMechSplitDeviatoric3D.o m_MEF90_DefMechSplitHydrostatic3D.o #m_MEF90_DefMechSplitMasonry3D.o 
	@echo "      $< MEF90_DIM=3"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechSplitNone2D.o: ${BASE_DIR}/m_MEF90_DefMechSplitNone.F90 m_MEF90_DefMechSplitClass2D.o
	@echo "      $< MEF90_DIM=2"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechSplitNone3D.o: ${BASE_DIR}/m_MEF90_DefMechSplitNone.F90 m_MEF90_DefMechSplitClass3D.o
	@echo "      $< MEF90_DIM=3"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechSplitMasonry2D.o: ${BASE_DIR}/m_MEF90_DefMechSplitMasonry.F90 m_MEF90_DefMechSplitClass2D.o
	@echo "      $< MEF90_DIM=2"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechSplitMasonry3D.o: ${BASE_DIR}/m_MEF90_DefMechSplitMasonry.F90 m_MEF90_DefMechSplitClass3D.o
	@echo "      $< MEF90_DIM=3"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechSplitHD2D.o: ${BASE_DIR}/m_MEF90_DefMechSplitHD.F90 m_MEF90_DefMechSplitClass2D.o
	@echo "      $< MEF90_DIM=2"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechSplitHD3D.o: ${BASE_DIR}/m_MEF90_DefMechSplitHD.F90 m_MEF90_DefMechSplitClass3D.o
	@echo "      $< MEF90_DIM=3"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechSplitDeviatoric2D.o: ${BASE_DIR}/m_MEF90_DefMechSplitDeviatoric.F90 m_MEF90_DefMechSplitClass2D.o
	@echo "      $< MEF90_DIM=2"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechSplitDeviatoric3D.o: ${BASE_DIR}/m_MEF90_DefMechSplitDeviatoric.F90 m_MEF90_DefMechSplitClass3D.o
	@echo "      $< MEF90_DIM=3"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechSplitHydrostatic2D.o: ${BASE_DIR}/m_MEF90_DefMechSplitHydrostatic.F90 m_MEF90_DefMechSplitClass2D.o
	@echo "      $< MEF90_DIM=2"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechSplitHydrostatic3D.o: ${BASE_DIR}/m_MEF90_DefMechSplitHydrostatic.F90 m_MEF90_DefMechSplitClass3D.o
	@echo "      $< MEF90_DIM=3"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechAssembly2D.o: ${BASE_DIR}/m_MEF90_DefMechAssembly.F90 m_MEF90_DefMechSplit2D.o m_MEF90_DefMechPlasticity2D.o m_MEF90_DefMechAT.o m_MEF90_DefMechCtx.o
	@echo "      $< MEF90_DIM=2"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechAssembly3D.o: ${BASE_DIR}/m_MEF90_DefMechAssembly.F90 m_MEF90_DefMechSplit3D.o m_MEF90_DefMechPlasticity3D.o m_MEF90_DefMechAT.o m_MEF90_DefMechCtx.o
	@echo "      $< MEF90_DIM=3"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechCtx.o: ${BASE_DIR}/m_MEF90_DefMechCtx.F90
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<
	
m_MEF90_DefMechPlasticity2D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticity.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityNone2D.o \
                   m_MEF90_DefMechPlasticityCap2D.o m_MEF90_DefMechPlasticityCrystalBCC2D.o m_MEF90_DefMechPlasticityCrystalSingleSlip2D.o m_MEF90_DefMechPlasticityDruckerPragerCap2D.o m_MEF90_DefMechPlasticityGreen2D.o \
				   m_MEF90_DefMechPlasticityGurson2D.o m_MEF90_DefMechPlasticityHill2D.o \
				   m_MEF90_DefMechPlasticityTresca2D.o \
				   m_MEF90_DefMechPlasticityVonMises2D.o
	@echo "      $< MEF90_DIM=2"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechPlasticity3D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticity.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityNone3D.o \
                   m_MEF90_DefMechPlasticityCap3D.o m_MEF90_DefMechPlasticityCrystalBCC3D.o m_MEF90_DefMechPlasticityCrystalSingleSlip3D.o m_MEF90_DefMechPlasticityDruckerPragerCap3D.o m_MEF90_DefMechPlasticityGreen3D.o \
				   m_MEF90_DefMechPlasticityGurson3D.o m_MEF90_DefMechPlasticityHill3D.o \
				   m_MEF90_DefMechPlasticityTresca3D.o \
				   m_MEF90_DefMechPlasticityVonMises3D.o
	@echo "      $< MEF90_DIM=3"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechPlasticityCtx2D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityCtx.F90 m_MEF90_DefMechCtx.o
	@echo "      $< MEF90_DIM=2"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechPlasticityCtx3D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityCtx.F90 m_MEF90_DefMechCtx.o
	@echo "      $< MEF90_DIM=3"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechPlasticityNone2D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityNone.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx2D.o
	@echo "      $< MEF90_DIM=2"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechPlasticityNone3D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityNone.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx3D.o
	@echo "      $< MEF90_DIM=3"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechPlasticityCap2D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityCap.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx2D.o
	@echo "      $< MEF90_DIM=2"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechPlasticityCap3D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityCap.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx3D.o
	@echo "      $< MEF90_DIM=3"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechPlasticityCrystalBCC2D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityCrystalBCC.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx2D.o
	@echo "      $< MEF90_DIM=2"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechPlasticityCrystalBCC3D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityCrystalBCC.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx3D.o
	@echo "      $< MEF90_DIM=3"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechPlasticityCrystalSingleSlip2D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityCrystalSingleSlip.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx2D.o
	@echo "      $< MEF90_DIM=2"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechPlasticityCrystalSingleSlip3D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityCrystalSingleSlip.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx3D.o
	@echo "      $< MEF90_DIM=3"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechPlasticityDruckerPragerCap2D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityDruckerPragerCap.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx2D.o
	@echo "      $< MEF90_DIM=2"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechPlasticityDruckerPragerCap3D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityDruckerPragerCap.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx3D.o
	@echo "      $< MEF90_DIM=3"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechPlasticityHill2D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityHill.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx2D.o
	@echo "      $< MEF90_DIM=2"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechPlasticityHill3D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityHill.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx3D.o
	@echo "      $< MEF90_DIM=3"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechPlasticityGreen2D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityGreen.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx2D.o
	@echo "      $< MEF90_DIM=2"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechPlasticityGreen3D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityGreen.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx3D.o
	@echo "      $< MEF90_DIM=3"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechPlasticityGurson2D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityGurson.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx2D.o
	@echo "      $< MEF90_DIM=2"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechPlasticityGurson3D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityGurson.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx3D.o
	@echo "      $< MEF90_DIM=3"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechPlasticityTresca2D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityTresca.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx2D.o
	@echo "      $< MEF90_DIM=2"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechPlasticityTresca3D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityTresca.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx3D.o
	@echo "      $< MEF90_DIM=3"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DefMechPlasticityVonMises2D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityVonMises.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx2D.o
	@echo "      $< MEF90_DIM=2"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DefMechPlasticityVonMises3D.o: ${BASE_DIR}/m_MEF90_DefMechPlasticityVonMises.F90 m_MEF90_DefMechCtx.o m_MEF90_DefMechPlasticityCtx3D.o
	@echo "      $< MEF90_DIM=3"
	@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<


MEF90:
	make -C ${MEF90_DIR} MEF90
allclean:
	make -C ${MEF90_DIR} clean