all:
	-@make -C ${MEF90_DIR} MEF90

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${MEF90_DIR}/Makefile.include

MEF90: m_MEF90.o

m_MEF90.o: ${MEF90_DIR}/MEF90/m_MEF90.F90 m_MEF_Ctx.o m_MEF_Diffusion.o m_MEF_Elements.o  \
           m_MEF_EXO.o m_MEF_LinAlg.o m_MEF_MassMatrix.o m_MEF_Materials.o m_MEF_MPI.o m_MEF_Norm.o \
           m_MEF_Parameters.o m_MEF_Sieve.o m_MEF_Utils.o 
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Ctx.o: ${MEF90_DIR}/MEF90/m_MEF_Ctx.F90 m_MEF_Parameters.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Diffusion.o: ${MEF90_DIR}/MEF90/m_MEF_Diffusion.F90 m_MEF_Elements.o m_MEF_Parameters.o m_MEF_Sieve.o m_MEF_Utils.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Elements.o: ${MEF90_DIR}/MEF90/m_MEF_Elements.F90 m_MEF_LinAlg.o m_MEF_Parameters.o m_MEF_Utils.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_EXO.o: ${MEF90_DIR}/MEF90/m_MEF_EXO.F90 m_MEF_Ctx.o m_MEF_Elements.o m_MEF_Parameters.o m_MEF_Sieve.o m_MEF_Utils.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_LinAlg.o: ${MEF90_DIR}/MEF90/m_MEF_LinAlg.F90 m_MEF_Parameters.o m_MEF_Utils.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_MassMatrix.o: ${MEF90_DIR}/MEF90/m_MEF_MassMatrix.F90 m_MEF_Elements.o m_MEF_Parameters.o m_MEF_Utils.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Materials.o: ${MEF90_DIR}/MEF90/m_MEF_Materials.F90 m_MEF_LinAlg.o m_MEF_Parameters.o m_MEF_Utils.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_MPI.o: ${MEF90_DIR}/MEF90/m_MEF_MPI.F90 m_MEF_Parameters.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Norm.o: ${MEF90_DIR}/MEF90/m_MEF_Norm.F90 m_MEF_Elements.o m_MEF_Parameters.o m_MEF_Sieve.o m_MEF_Utils.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Parameters.o: ${MEF90_DIR}/MEF90/m_MEF_Parameters.F90
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Sieve.o: ${MEF90_DIR}/MEF90/m_MEF_Sieve.F90 m_MEF_LinAlg.o m_MEF_Parameters.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Utils.o: ${MEF90_DIR}/MEF90/m_MEF_Utils.F90 m_MEF_MPI.o m_MEF_Parameters.o
	@echo "      $<"
	-@${PETSC_MAKE_STOP_ON_ERROR}${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

