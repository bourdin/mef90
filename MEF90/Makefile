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
        m_MEF90_Elements.o \
        m_MEF90_EXO.o \
        m_MEF90_DMPlex.o \
        m_MEF90_Utils.o \
        m_MEF90_MassMatrixImplementation_MEF90Element2D_Scal.o \
        m_MEF90_MassMatrixImplementation_MEF90Element2D_Vect.o \
        m_MEF90_MassMatrixImplementation_MEF90Element3D_Scal.o \
        m_MEF90_MassMatrixImplementation_MEF90Element3D_Vect.o \
        m_MEF90_MassMatrixInterface.o \
        m_MEF90_NormsImplementation_MEF90Element2D_Scal.o \
        m_MEF90_NormsImplementation_MEF90Element2D_Vect.o \
        m_MEF90_NormsImplementation_MEF90Element3D_Scal.o \
        m_MEF90_NormsImplementation_MEF90Element3D_Vect.o \
        m_MEF90_NormsInterface.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_Ctx.o: ${MEF90_DIR}/MEF90/m_MEF90_Ctx.F90 \
        m_MEF90_Parameters.o m_MEF90_Utils.o m_MEF90_Elements.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_DMPlex.o: ${MEF90_DIR}/MEF90/m_MEF90_DMPlex.F90 \
        m_MEF90_Elements.o m_MEF90_Ctx.o
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
        m_MEF90_Utils.o \
        m_MEF90_DMPlex.o
	@echo "      $< MEF90Element2DScal"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2DScal -c -o $@ $<
m_MEF90_MassMatrixImplementation_MEF90Element2D_Vect.o: ${MEF90_DIR}/MEF90/m_MEF90_MassMatrixImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o \
        m_MEF90_DMPlex.o
	@echo "      $< MEF90Element2DVect"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2DVect -c -o $@ $<
m_MEF90_MassMatrixImplementation_MEF90Element3D_Scal.o: ${MEF90_DIR}/MEF90/m_MEF90_MassMatrixImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o \
        m_MEF90_DMPlex.o
	@echo "      $< MEF90Element3DScal"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3DScal -c -o $@ $<
m_MEF90_MassMatrixImplementation_MEF90Element3D_Vect.o: ${MEF90_DIR}/MEF90/m_MEF90_MassMatrixImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o \
        m_MEF90_DMPlex.o
	@echo "      $< MEF90Element3DVect"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3DVect -c -o $@ $<

m_MEF90_MassMatrixInterface.o: ${MEF90_DIR}/MEF90/m_MEF90_MassMatrixInterface.F90 m_MEF90_LinAlg.o \
        m_MEF90_MassMatrixImplementation_MEF90Element2D_Scal.o \
        m_MEF90_MassMatrixImplementation_MEF90Element2D_Vect.o \
        m_MEF90_MassMatrixImplementation_MEF90Element3D_Scal.o \
        m_MEF90_MassMatrixImplementation_MEF90Element3D_Vect.o \
        m_MEF90_Parameters.o m_MEF90_Utils.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

###
### Norms for each element type
###
m_MEF90_NormsImplementation_MEF90Element2D_Scal.o: ${MEF90_DIR}/MEF90/m_MEF90_NormsImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o \
        m_MEF90_DMPlex.o
	@echo "      $< MEF90Element2DScal"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2DScal -DMEF90_ELEMENTTYPE_SCALAR -DMEF90_DIM=2 -c -o $@ $<

m_MEF90_NormsImplementation_MEF90Element2D_Vect.o: ${MEF90_DIR}/MEF90/m_MEF90_NormsImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o \
        m_MEF90_DMPlex.o
	@echo "      $< MEF90Element2DVect"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element2DVect -DMEF90_ELEMENTTYPE_VECT  -DMEF90_DIM=2 -c -o $@ $<
m_MEF90_NormsImplementation_MEF90Element3D_Scal.o: ${MEF90_DIR}/MEF90/m_MEF90_NormsImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o \
        m_MEF90_DMPlex.o
	@echo "      $< MEF90Element3DScal"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3DScal -DMEF90_ELEMENTTYPE_SCALAR -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_NormsImplementation_MEF90Element3D_Vect.o: ${MEF90_DIR}/MEF90/m_MEF90_NormsImplementation.F90 \
        m_MEF90_Elements.o \
        m_MEF90_Parameters.o \
        m_MEF90_Utils.o \
        m_MEF90_DMPlex.o
	@echo "      $< MEF90Element3DVect"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -DMEF90_ELEMENTTYPE=MEF90Element3DVect -DMEF90_ELEMENTTYPE_VECT  -DMEF90_DIM=3 -c -o $@ $<

m_MEF90_NormsInterface.o: ${MEF90_DIR}/MEF90/m_MEF90_NormsInterface.F90 m_MEF90_LinAlg.o \
        m_MEF90_NormsImplementation_MEF90Element2D_Scal.o \
        m_MEF90_NormsImplementation_MEF90Element2D_Vect.o \
        m_MEF90_NormsImplementation_MEF90Element3D_Scal.o \
        m_MEF90_NormsImplementation_MEF90Element3D_Vect.o \
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

m_MEF90_Parameters.o: ${MEF90_DIR}/MEF90/m_MEF90_Parameters.F90
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

m_MEF90_Utils.o: ${MEF90_DIR}/MEF90/m_MEF90_Utils.F90 \
        m_MEF90_MPI.o \
        m_MEF90_Parameters.o
	@echo "      $<"
	-@${FC} ${FC_FLAGS} ${FCPPFLAGS} -I${MEF90_INCLUDE} -c -o $@ $<

