MEF90_DIM=2
all:
	-@make -C ${MEF90_DIR} m_HeatXfer

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
BASE_DIR=${MEF90_DIR}/m_HeatXfer

m_HeatXfer:m_MEF90_HeatXfer.o

m_MEF90_HeatXfer.o: ${BASE_DIR}/m_MEF90_HeatXfer.F90 m_MEF90_HeatXferCtx.o m_MEF90_HeatXferAssembly2D.o  m_MEF90_HeatXferAssembly3D.o 
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_HeatXferAssembly2D.o: ${BASE_DIR}/m_MEF90_HeatXferAssembly.F90 m_MEF90_HeatXferCtx.o
	@echo "      $< MEF90_DIM=2"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_HeatXferAssembly3D.o: ${BASE_DIR}/m_MEF90_HeatXferAssembly.F90 m_MEF90_HeatXferCtx.o
	@echo "      $< MEF90_DIM=3"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_HeatXferCtx.o: ${BASE_DIR}/m_MEF90_HeatXferCtx.F90
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

MEF90:
	make -C ${MEF90_DIR} MEF90
allclean:
	make -C ${MEF90_DIR} clean