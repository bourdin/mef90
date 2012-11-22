MEF90_DIM=2
all:
	-@make -C ${MEF90_DIR} m_Elasticity

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
BASE_DIR=${MEF90_DIR}/m_Elasticity

m_Elasticity:m_MEF90_Elasticity.o

m_MEF90_Elasticity.o: ${BASE_DIR}/m_MEF90_Elasticity.F90 m_MEF90_ElasticityAssembly2D.o  m_MEF90_ElasticityAssembly3D.o m_MEF90_ElasticityCtx.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_ElasticityAssembly2D.o: ${BASE_DIR}/m_MEF90_ElasticityAssembly.F90 m_MEF90_ElasticityCtx.o
	@echo "      $< MEF90_DIM=2"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_ElasticityAssembly3D.o: ${BASE_DIR}/m_MEF90_ElasticityAssembly.F90 m_MEF90_ElasticityCtx.o
	@echo "      $< MEF90_DIM=3"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_ElasticityCtx.o: ${BASE_DIR}/m_MEF90_ElasticityCtx.F90
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

MEF90:
	make -C ${MEF90_DIR} MEF90
allclean:
	make -C ${MEF90_DIR} clean