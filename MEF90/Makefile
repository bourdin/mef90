all:
	-@make -C ${MEF90_DIR} MEF90

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ${MEF90_DIR}/Makefile.include

MEF90: m_MEF90.o

m_MEF90.o: ${MEF90_DIR}/MEF90/m_MEF90.F90 m_MEF_Ctx.o \
           m_MEF_DiffusionImplementation_2D.o \
           m_MEF_DiffusionImplementation_3D.o \
           m_MEF_DiffusionInterface.o \
           m_MEF_Elements.o \
           m_MEF_EXO.o \
           m_MEF_LinAlg.o \
           m_MEF_MassMatrixImplementation_MEF90Element2D_Scal.o \
           m_MEF_MassMatrixImplementation_MEF90Element2D_Vect.o \
           m_MEF_MassMatrixImplementation_MEF90Element2D_Elast.o \
           m_MEF_MassMatrixImplementation_MEF90Element3D_Scal.o \
           m_MEF_MassMatrixImplementation_MEF90Element3D_Vect.o \
           m_MEF_MassMatrixImplementation_MEF90Element3D_Elast.o \
           m_MEF_MassMatrixInterface.o \
           m_MEF_Materials.o \
           m_MEF_MPI.o \
           m_MEF_Norm.o \
           m_MEF_Parameters.o \
           m_MEF_Sieve.o \
           m_MEF_Utils.o 
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Ctx.o: ${MEF90_DIR}/MEF90/m_MEF_Ctx.F90 \
        m_MEF_Parameters.o
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_DiffusionImplementation_2D.o: ${MEF90_DIR}/MEF90/m_MEF_DiffusionImplementation.F90 \
        m_MEF_Elements.o \
        m_MEF_Parameters.o \
        m_MEF_Sieve.o \
        m_MEF_Utils.o
	@echo "      $< 2D / MEF90Element2D_Scal"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2D_Scal -DMEF90_DIM=2 -c -o $@ $<

m_MEF_DiffusionImplementation_3D.o: ${MEF90_DIR}/MEF90/m_MEF_DiffusionImplementation.F90 \
        m_MEF_Elements.o \
        m_MEF_Parameters.o \
        m_MEF_Sieve.o \
        m_MEF_Utils.o
	@echo "      $< 3D / MEF90Element3D_Scal"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3D_Scal -DMEF90_DIM=3 -c -o $@ $<

m_MEF_DiffusionInterface.o: ${MEF90_DIR}/MEF90/m_MEF_DiffusionInterface.F90 \
        m_MEF_DiffusionImplementation_2D.o \
        m_MEF_DiffusionImplementation_3D.o \
        m_MEF_Elements.o \
        m_MEF_Parameters.o \
        m_MEF_Sieve.o \
        m_MEF_Utils.o
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Elements.o: ${MEF90_DIR}/MEF90/m_MEF_Elements.F90 \
        m_MEF_LinAlg.o \
        m_MEF_Parameters.o \
        m_MEF_Utils.o
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_EXO.o: ${MEF90_DIR}/MEF90/m_MEF_EXO.F90 m_MEF_Ctx.o \
        m_MEF_Elements.o \
        m_MEF_Parameters.o \
        m_MEF_Sieve.o \
        m_MEF_Utils.o
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_LinAlg.o: ${MEF90_DIR}/MEF90/m_MEF_LinAlg.F90 \
        m_MEF_Parameters.o \
        m_MEF_Utils.o
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

###
### MassMatrix for each element type
###
m_MEF_MassMatrixImplementation_MEF90Element2D_Scal.o: ${MEF90_DIR}/MEF90/m_MEF_MassMatrixImplementation.F90 \
        m_MEF_Elements.o \
        m_MEF_Parameters.o \
        m_MEF_Utils.o
	@echo "      $< MEF90Element2D_Scal"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2D_Scal -c -o $@ $<
m_MEF_MassMatrixImplementation_MEF90Element2D_Vect.o: ${MEF90_DIR}/MEF90/m_MEF_MassMatrixImplementation.F90 \
        m_MEF_Elements.o \
        m_MEF_Parameters.o \
        m_MEF_Utils.o
	@echo "      $< MEF90Element2D_Vect"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2D_Vect -c -o $@ $<
m_MEF_MassMatrixImplementation_MEF90Element2D_Elast.o: ${MEF90_DIR}/MEF90/m_MEF_MassMatrixImplementation.F90 \
        m_MEF_Elements.o \
        m_MEF_Parameters.o \
        m_MEF_Utils.o
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2D_Elast -c -o $@ $<
m_MEF_MassMatrixImplementation_MEF90Element3D_Scal.o: ${MEF90_DIR}/MEF90/m_MEF_MassMatrixImplementation.F90 \
        m_MEF_Elements.o \
        m_MEF_Parameters.o \
        m_MEF_Utils.o
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3D_Scal -c -o $@ $<
m_MEF_MassMatrixImplementation_MEF90Element3D_Vect.o: ${MEF90_DIR}/MEF90/m_MEF_MassMatrixImplementation.F90 \
        m_MEF_Elements.o \
        m_MEF_Parameters.o \
        m_MEF_Utils.o
	@echo "      $< MEF90Element3D_Vect"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3D_Vect -c -o $@ $<
m_MEF_MassMatrixImplementation_MEF90Element3D_Elast.o: ${MEF90_DIR}/MEF90/m_MEF_MassMatrixImplementation.F90 \
        m_MEF_Elements.o \
        m_MEF_Parameters.o \
        m_MEF_Utils.o
	@echo "      $< MEF90Element3D_Elast"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3D_Elast -c -o $@ $<

m_MEF_MassMatrixInterface.o: ${MEF90_DIR}/MEF90/m_MEF_MassMatrixInterface.F90 m_MEF_LinAlg.o \
        m_MEF_MassMatrixImplementation_MEF90Element2D_Scal.o \
        m_MEF_MassMatrixImplementation_MEF90Element2D_Vect.o \
        m_MEF_MassMatrixImplementation_MEF90Element2D_Elast.o \
        m_MEF_MassMatrixImplementation_MEF90Element3D_Scal.o \
        m_MEF_MassMatrixImplementation_MEF90Element3D_Vect.o \
        m_MEF_MassMatrixImplementation_MEF90Element3D_Elast.o \
        m_MEF_Parameters.o m_MEF_Utils.o
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Materials.o: ${MEF90_DIR}/MEF90/m_MEF_Materials.F90 \
        m_MEF_LinAlg.o \
        m_MEF_Parameters.o \
        m_MEF_Utils.o
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_MPI.o: ${MEF90_DIR}/MEF90/m_MEF_MPI.F90 \
        m_MEF_Parameters.o
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Norm.o: ${MEF90_DIR}/MEF90/m_MEF_Norm.F90 \
        m_MEF_Elements.o \
        m_MEF_Parameters.o \
        m_MEF_Sieve.o \
        m_MEF_Utils.o
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Parameters.o: ${MEF90_DIR}/MEF90/m_MEF_Parameters.F90
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Sieve.o: ${MEF90_DIR}/MEF90/m_MEF_Sieve.F90 \
        m_MEF_Elements.o \
        m_MEF_LinAlg.o \
        m_MEF_Parameters.o
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF_Utils.o: ${MEF90_DIR}/MEF90/m_MEF_Utils.F90 \
        m_MEF_MPI.o \
        m_MEF_Parameters.o
	@echo "      $<"
	-@${FC} ${FCFLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

