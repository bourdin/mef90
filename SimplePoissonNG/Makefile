MEF90_DIM=2
#all: MEF90 ${MEF90_DIR}/bin/${PETSC_ARCH}/SimplePoissonNG${MEF90_DIM}D# SimplePoissonNG3D
all:
	-@make -C ${MEF90_DIR} SimplePoissonNG

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables
include ${MEF90_DIR}/Makefile.include
include ${MEF90_DIR}/MEF90/Makefile.include
BASE_DIR=${MEF90_DIR}/SimplePoissonNG

SimplePoissonNG:${MEF90_DIR}/bin/${PETSC_ARCH}/SimplePoissonNG${MEF90_DIM}D
	-@echo "      $@"

${MEF90_DIR}/bin/${PETSC_ARCH}/SimplePoissonNG${MEF90_DIM}D: SimplePoissonNG${MEF90_DIM}D.o
	-@echo "      $@"
	-@${FLINKER} -o $@ SimplePoissonNG${MEF90_DIM}D.o m_Poisson${MEF90_DIM}D.o m_PoissonGlobalCtx${MEF90_DIM}D.o m_PoissonCellSetCtx${MEF90_DIM}D.o m_PoissonVertexSetCtx${MEF90_DIM}D.o m_PoissonAssembly${MEF90_DIM}D.o ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

SimplePoissonNG${MEF90_DIM}D.o: ${BASE_DIR}/SimplePoissonNG.F90 m_Poisson${MEF90_DIM}D.o m_PoissonGlobalCtx${MEF90_DIM}D.o m_PoissonCellSetCtx${MEF90_DIM}D.o m_PoissonVertexSetCtx${MEF90_DIM}D.o m_PoissonAssembly${MEF90_DIM}D.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=${MEF90_DIM} -c -o $@ $<

m_PoissonGlobalCtx${MEF90_DIM}D.o: ${BASE_DIR}/m_PoissonGlobalCtx.F90
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=${MEF90_DIM} -c -o $@ $<

m_PoissonCellSetCtx${MEF90_DIM}D.o: ${BASE_DIR}/m_PoissonCellSetCtx.F90
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=${MEF90_DIM} -c -o $@ $<

m_PoissonVertexSetCtx${MEF90_DIM}D.o: ${BASE_DIR}/m_PoissonVertexSetCtx.F90
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=${MEF90_DIM} -c -o $@ $<

m_PoissonAssembly${MEF90_DIM}D.o: ${BASE_DIR}/m_PoissonAssembly.F90 m_PoissonGlobalCtx${MEF90_DIM}D.o m_PoissonCellSetCtx${MEF90_DIM}D.o m_PoissonVertexSetCtx${MEF90_DIM}D.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=${MEF90_DIM} -c -o $@ $<

m_Poisson${MEF90_DIM}D.o: ${BASE_DIR}/m_Poisson.F90 m_PoissonGlobalCtx${MEF90_DIM}D.o m_PoissonCellSetCtx${MEF90_DIM}D.o m_PoissonVertexSetCtx${MEF90_DIM}D.o m_PoissonAssembly${MEF90_DIM}D.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=${MEF90_DIM} -c -o $@ $<

MEF90:
	make -C ${MEF90_DIR} MEF90
allclean:
	make -C ${MEF90_DIR} clean
