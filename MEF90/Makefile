all:
	-@make -C ${MEF90_DIR} MEF90

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${MEF90_DIR}/Makefile.include

m_MEF90_Parameters.F90: ${MEF90_DIR}/mef90version.h

MEF90: m_MEF90.o

m_MEF90.o: ${MEF90_DIR}/MEF90/m_MEF90.F90 \
	       m_MEF90_Ctx.o \
           m_MEF90_LinAlg.o \
           m_MEF90_MPI.o \
           m_MEF90_Parameters.o \
           m_MEF90_Materials.o \
           m_MEF90_Utils.o 
#           m_MEF90_DiffusionImplementation_2D.o \
#           m_MEF90_DiffusionImplementation_3D.o \
#           m_MEF90_DiffusionInterface.o \
#           m_MEF90_ElasticityImplementation_2D.o \
#           m_MEF90_ElasticityImplementation_3D.o \
#           m_MEF90_ElasticityInterface.o \
#           m_MEF90_GradDamageImplementation_2D.o \
#           m_MEF90_GradDamageImplementation_3D.o \
#           m_MEF90_GradDamageInterface.o \
#           m_MEF90_Elements.o \
#           m_MEF90_EXO.o \
#           m_MEF90_MassMatrixImplementation_MEF90Element2D_Scal.o \
#           m_MEF90_MassMatrixImplementation_MEF90Element2D_Vect.o \
#           m_MEF90_MassMatrixImplementation_MEF90Element2D_Elast.o \
#           m_MEF90_MassMatrixImplementation_MEF90Element3D_Scal.o \
#           m_MEF90_MassMatrixImplementation_MEF90Element3D_Vect.o \
#           m_MEF90_MassMatrixImplementation_MEF90Element3D_Elast.o \
#           m_MEF90_MassMatrixInterface.o \
#           m_MEF90_Norm.o \
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_Ctx.o: ${MEF90_DIR}/MEF90/m_MEF90_Ctx.F90 \
        m_MEF90_Parameters.o m_MEF90_Utils.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_DiffusionImplementation_2D.o: ${MEF90_DIR}/MEF90/m_MEF90_DiffusionImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $< 2D / MEF90Element2D_Scal"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2D_Scal -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_DiffusionImplementation_3D.o: ${MEF90_DIR}/MEF90/m_MEF90_DiffusionImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $< 3D / MEF90Element3D_Scal"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3D_Scal -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_DiffusionInterface.o: ${MEF90_DIR}/MEF90/m_MEF90_DiffusionInterface.F90 \
        m_MEF90_DiffusionImplementation_2D.o \
        m_MEF90_DiffusionImplementation_3D.o \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_ElasticityImplementation_2D.o: ${MEF90_DIR}/MEF90/m_MEF90_ElasticityImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Materials.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $< 2D"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2D_Elast -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_ElasticityImplementation_3D.o: ${MEF90_DIR}/MEF90/m_MEF90_ElasticityImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Materials.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $< 3D"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3D_Elast -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_ElasticityInterface.o: ${MEF90_DIR}/MEF90/m_MEF90_ElasticityInterface.F90 \
        m_MEF90_ElasticityImplementation_2D.o \
        m_MEF90_ElasticityImplementation_3D.o \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_GradDamageImplementation_2D.o: ${MEF90_DIR}/MEF90/m_MEF90_GradDamageImplementation.F90\
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $< MEF90_DIM=2"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_GradDamageImplementation_3D.o: ${MEF90_DIR}/MEF90/m_MEF90_GradDamageImplementation.F90\
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $< MEF90_DIM=3"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_GradDamageInterface.o: ${MEF90_DIR}/MEF90/m_MEF90_GradDamageInterface.F90 \
                                   m_MEF90_GradDamageImplementation_2D.o \
                                   m_MEF90_GradDamageImplementation_3D.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_HookesLaw.o: ${MEF90_DIR}/MEF90/m_MEF90_HookesLaw.F90 m_MEF90_LinAlg.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_Elements.o: ${MEF90_DIR}/MEF90/m_MEF90_Elements.F90 \
        m_MEF90_LinAlg.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_EXO.o: ${MEF90_DIR}/MEF90/m_MEF90_EXO.F90 m_MEF90_Ctx.o \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_LinAlg.o: ${MEF90_DIR}/MEF90/m_MEF90_LinAlg.F90 \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

###
### MassMatrix for each element type
###
m_MEF90_MassMatrixImplementation_MEF90Element2D_Scal.o: ${MEF90_DIR}/MEF90/m_MEF90_MassMatrixImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $< MEF90Element2D_Scal"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2D_Scal -c -o $@ $<
m_MEF90_MassMatrixImplementation_MEF90Element2D_Vect.o: ${MEF90_DIR}/MEF90/m_MEF90_MassMatrixImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $< MEF90Element2D_Vect"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2D_Vect -c -o $@ $<
m_MEF90_MassMatrixImplementation_MEF90Element2D_Elast.o: ${MEF90_DIR}/MEF90/m_MEF90_MassMatrixImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $< MEF90Element2D_Elast"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2D_Elast -c -o $@ $<
m_MEF90_MassMatrixImplementation_MEF90Element3D_Scal.o: ${MEF90_DIR}/MEF90/m_MEF90_MassMatrixImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $< MEF90Element3D_Scal"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3D_Scal -c -o $@ $<
m_MEF90_MassMatrixImplementation_MEF90Element3D_Vect.o: ${MEF90_DIR}/MEF90/m_MEF90_MassMatrixImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $< MEF90Element3D_Vect"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3D_Vect -c -o $@ $<
m_MEF90_MassMatrixImplementation_MEF90Element3D_Elast.o: ${MEF90_DIR}/MEF90/m_MEF90_MassMatrixImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $< MEF90Element3D_Elast"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3D_Elast -c -o $@ $<

m_MEF90_MassMatrixInterface.o: ${MEF90_DIR}/MEF90/m_MEF90_MassMatrixInterface.F90 m_MEF90_LinAlg.o \
        m_MEF90_MassMatrixImplementation_MEF90Element2D_Scal.o \
        m_MEF90_MassMatrixImplementation_MEF90Element2D_Vect.o \
        m_MEF90_MassMatrixImplementation_MEF90Element2D_Elast.o \
        m_MEF90_MassMatrixImplementation_MEF90Element3D_Scal.o \
        m_MEF90_MassMatrixImplementation_MEF90Element3D_Vect.o \
        m_MEF90_MassMatrixImplementation_MEF90Element3D_Elast.o \
        m_MEF90_Parameters.o m_MEF90_Utils.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_Materials.o: ${MEF90_DIR}/MEF90/m_MEF90_Materials.F90 \
        m_MEF90_Ctx.o \
        m_MEF90_LinAlg.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_MPI.o: ${MEF90_DIR}/MEF90/m_MEF90_MPI.F90 \
        m_MEF90_Parameters.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_Norm.o: ${MEF90_DIR}/MEF90/m_MEF90_Norm.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_Parameters.o: ${MEF90_DIR}/MEF90/m_MEF90_Parameters.F90
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_Utils.o: ${MEF90_DIR}/MEF90/m_MEF90_Utils.F90 \
        m_MEF90_MPI.o \
        m_MEF90_Parameters.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

