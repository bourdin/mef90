#include "../MEF90/mef90.inc"
module m_MEF90_DefMech_class
#include "petsc/finclude/petsc.h"
   use petscsys
   use m_MEF90_BaseClass
   use m_MEF90_Ctx
   use m_MEF90_DMPlex
   use, intrinsic :: iso_c_binding
   implicit none(type)

!!!
!!!
!!!  MEF90DefMech_Type: The class holding the state and the options of a defect mechanics problem
!!!
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu
!!!      2026    Blaise Bourdin bourdin@mcmaster.ca
!!!

   enum, bind(c)
      enumerator :: MEF90DefMech_SolverTypeAltMin = 0, &
         MEF90DefMech_SolverTypeQuasiNewton1, &
         MEF90DefMech_SolverTypeQuasiNewton2
   end enum
   character(len=MEF90MXSTRLEN), dimension(6), protected   :: MEF90DefMech_SolverTypeList = [character(len=MEF90MXSTRLEN) :: &
      'AltMin', &
      'QuasiNewton1', &
      'QuasiNewton2', &
      'MEF90DefMech_SolverType', &
      '_MEF90DefMech_SolverType', &
      '']

   enum, bind(c)
      enumerator :: MEF90DefMech_TimeSteppingTypeNULL = 0, &
         MEF90DefMech_TimeSteppingTypeQuasiStatic
   end enum
   character(len=MEF90MXSTRLEN), dimension(5), protected   :: MEF90DefMech_TimeSteppingTypeList = [character(len=MEF90MXSTRLEN) :: &
      'Null', &
      'QuasiStatic', &
      'MEF90DefMech_TimeSteppingType', &
      '_MEF90DefMech_TimeSteppingType', &
      '']

   enum, bind(c)
      enumerator :: MEF90DefMech_DamageSolverTypeSNES = 0, &
         MEF90DefMech_DamageSolverTypeTAO
   end enum
   character(len=MEF90MXSTRLEN), dimension(5), protected   :: MEF90DefMech_DamageSolverTypeList = [character(len=MEF90MXSTRLEN) :: &
      'SNES', &
      'Tao', &
      'MEF90DefMech_DamageSolverType', &
      '_MEF90DefMech_DamageSolverType', &
      '']

   enum, bind(c)
      enumerator :: MEF90DefMech_BTTypeNULL = 0, &
         MEF90DefMech_BTTypeBackward, &
         MEF90DefMech_BTTypeForward
   end enum
   character(len=MEF90MXSTRLEN), dimension(6), protected   :: MEF90DefMech_BTTypeList = [character(len=MEF90MXSTRLEN) :: &
      'Null', &
      'Backward', &
      'Forward', &
      'MEF90DefMech_BTType', &
      '_MEF90DefMech_BTType', &
      '']

   enum, bind(c)
      enumerator :: MEF90DefMech_plasticityTypeNone = 0, &
         MEF90DefMech_plasticityTypeTresca, &
         MEF90DefMech_plasticityTypeVonMises, &
         MEF90DefMech_plasticityTypeVonMisesPlaneTheory, &
         MEF90DefMech_plasticityTypeCapModel, &
         MEF90DefMech_plasticityTypeDruckerPragerCapModel, &
         MEF90DefMech_plasticityTypeVonMises1D, &
         MEF90DefMech_plasticityTypeHillPlaneTheory, &
         MEF90DefMech_PlasticityTypeGreen, &
         MEF90DefMech_PlasticityTypeGurson
   end enum
   character(len=MEF90MXSTRLEN), dimension(13), protected   :: MEF90DefMech_plasticityTypeList = [character(len=MEF90MXSTRLEN) :: &
      'None', &
      'Tresca', &
      'VonMises', &
      'VonMisesPlaneTheory', &
      'CapModel', &
      'DruckerPragerCapModel', &
      'VonMises1D', &
      'HillPlaneTheory', &
      'Green', &
      'Gurson', &
      'MEF90DefMech_plasticityType', &
      '_MEF90DefMech_plasticityType', &
      '']

   enum, bind(c)
      enumerator :: MEF90DefMech_unilateralContactTypeNone = 0, &
         MEF90DefMech_unilateralContactTypeHydrostaticDeviatoric, &
         MEF90DefMech_unilateralContactTypeHydrostatic, &
         MEF90DefMech_unilateralContactTypeDeviatoric, &
         MEF90DefMech_unilateralContactTypePrincipalStrains, &
         MEF90DefMech_unilateralContactTypeMasonry
   end enum
   character(len=MEF90MXSTRLEN), dimension(9), protected   :: MEF90DefMech_unilateralContactTypeList = [character(len=MEF90MXSTRLEN) :: &
      'None', &
      'HydrostaticDeviatoric', &
      'Hydrostatic', &
      'Deviatoric', &
      'PrincipalStrains', &
      'Masonry', &
      'MEF90DefMech_unilateralContactTypeList', &
      '_MEF90DefMech_unilateralContactTypeList', &
      '']

!!!
!!!  MEF90DefMechGlobalOptions_Type: the problem-wide options of a MEF90DefMech_Type.
!!!                                  The values given here are the defaults used by setFromOptions
!!!
   type MEF90DefMechGlobalOptions_Type
      PetscEnum                              :: timeSteppingType = MEF90DefMech_TimeSteppingTypeQuasiStatic
      PetscEnum                              :: solverType = MEF90DefMech_SolverTypeAltMin
      PetscEnum                              :: damageSolverType = MEF90DefMech_DamageSolverTypeSNES

      !!! scaling = time (step) scaling law currently CST, Linear, or File
      PetscEnum                              :: boundaryDisplacementScaling = MEF90Scaling_Linear
      PetscEnum                              :: displacementLowerBoundScaling = MEF90Scaling_CST
      PetscEnum                              :: displacementUpperBoundScaling = MEF90Scaling_CST
      PetscEnum                              :: cohesiveDisplacementScaling = MEF90Scaling_CST
      PetscEnum                              :: boundaryDamageScaling = MEF90Scaling_CST
      PetscEnum                              :: bodyForceScaling = MEF90Scaling_Linear
      PetscEnum                              :: boundaryForceScaling = MEF90Scaling_Linear
      PetscEnum                              :: pressureForceScaling = MEF90Scaling_Linear
      PetscEnum                              :: CrackPressureScaling = MEF90Scaling_Linear

      PetscReal                              :: damageATol = 1.0e-3_kr
      PetscInt                               :: damageMaxIt = 1000_ki
      PetscInt                               :: PCLag = 10_ki
      PetscReal                              :: SOROmega = 1.0_kr
      PetscReal                              :: irrevthres = 0.0_kr
      PetscEnum                              :: BTType = MEF90DefMech_BTTypeNULL
      PetscInt                               :: BTInterval = -1_ki
      PetscInt                               :: BTScope = -1_ki
      PetscReal                              :: BTTol = 1.0e-2_kr
      PetscReal                              :: plasticStrainATol = 1.0e-4_kr
      PetscReal                              :: InjectedVolumeATol = 1.0e-3_kr
      PetscReal                              :: dampingCoefficientDisplacement = 0.0_kr
      PetscReal                              :: dampingCoefficientDamage = 0.0_kr
      PetscBool                              :: temperatureExport = PETSC_FALSE
      PetscBool                              :: displacementExport = PETSC_TRUE
      PetscBool                              :: damageExport = PETSC_TRUE
      PetscBool                              :: stressExport = PETSC_TRUE
      PetscBool                              :: plasticStrainExport = PETSC_FALSE
      PetscBool                              :: cumulatedPlasticDissipationExport = PETSC_FALSE
   end type MEF90DefMechGlobalOptions_Type

   type MEF90DefMechCellSetOptions_Type
      PetscReal, dimension(3)                 :: bodyforce = [0.0_kr, 0.0_kr, 0.0_kr]
      PetscReal                               :: crackPressure = 0.0_kr
      PetscEnum                               :: plasticityType = MEF90DefMech_plasticityTypeNone
      PetscReal, dimension(3)                 :: cohesiveDisplacement = [0.0_kr, 0.0_kr, 0.0_kr]
      PetscBool, dimension(3)                 :: Has_displacementBC = [PETSC_FALSE, PETSC_FALSE, PETSC_FALSE]
      PetscReal, dimension(3)                 :: boundaryDisplacement = [0.0_kr, 0.0_kr, 0.0_kr]
      PetscReal, dimension(3)                 :: displacementLowerBound = [MEF90NINFINITY, MEF90NINFINITY, MEF90NINFINITY]
      PetscReal, dimension(3)                 :: displacementUpperBound = [MEF90INFINITY, MEF90INFINITY, MEF90INFINITY]
      PetscBool                               :: Has_damageBC = PETSC_FALSE
      PetscReal                               :: boundaryDamage = 0.0_kr
      PetscBool                               :: CrackVolumeControlled = PETSC_FALSE
      PetscBool                               :: WorkControlled = PETSC_FALSE
   end type MEF90DefMechCellSetOptions_Type

   type MEF90DefMechFaceSetOptions_Type
      PetscReal, dimension(3)                 :: boundaryforce = [0.0_kr, 0.0_kr, 0.0_kr]
      PetscReal                               :: pressureForce = 0.0_kr
      PetscBool, dimension(3)                 :: Has_displacementBC = [PETSC_FALSE, PETSC_FALSE, PETSC_FALSE]
      PetscReal, dimension(3)                 :: boundaryDisplacement = [0.0_kr, 0.0_kr, 0.0_kr]
      PetscReal, dimension(3)                 :: displacementLowerBound = [MEF90NINFINITY, MEF90NINFINITY, MEF90NINFINITY]
      PetscReal, dimension(3)                 :: displacementUpperBound = [MEF90INFINITY, MEF90INFINITY, MEF90INFINITY]
      PetscBool                               :: Has_damageBC = PETSC_FALSE
      PetscReal                               :: boundaryDamage = 0.0_kr
   end type MEF90DefMechFaceSetOptions_Type

   type MEF90DefMechVertexSetOptions_Type
      PetscBool, dimension(3)                 :: Has_displacementBC = [PETSC_FALSE, PETSC_FALSE, PETSC_FALSE]
      PetscReal, dimension(3)                 :: boundaryDisplacement = [0.0_kr, 0.0_kr, 0.0_kr]
      PetscReal, dimension(3)                 :: displacementLowerBound = [MEF90NINFINITY, MEF90NINFINITY, MEF90NINFINITY]
      PetscReal, dimension(3)                 :: displacementUpperBound = [MEF90INFINITY, MEF90INFINITY, MEF90INFINITY]
      PetscBool                               :: Has_damageBC = PETSC_FALSE
      PetscReal                               :: boundaryDamage = 0.0_kr
   end type MEF90DefMechVertexSetOptions_Type

   type, extends(MEF90Object) :: MEF90DefMech_Type
      type(MEF90Ctx_Type), pointer             :: MEF90Ctx => null()
      type(tDM)                                :: megaDM
      PetscInt                                 :: dim
      PetscReal                                :: analysisTime, timeStep

      !!!  vertex based vec
      type(tVec), pointer                       :: displacementLocal => null()
      type(tVec), pointer                       :: displacementPreviousStepLocal => null()
      type(tVec), pointer                       :: damageLocal => null()
      type(tVec), pointer                       :: damagePreviousStepLocal => null()
      type(tVec), pointer                       :: displacementLowerBoundLocal => null()
      type(tVec), pointer                       :: displacementUpperBoundLocal => null()
      type(tVec), pointer                       :: temperatureLocal => null()

      !!! cell based vec
      type(tVec), pointer                       :: bodyForce => null()
      type(tVec), pointer                       :: boundaryForce => null()
      type(tVec), pointer                       :: pressureForce => null()
      !Type(tVec),pointer                       :: crackPressure => null()
      type(tVec), pointer                       :: cohesiveDisplacement => null()
      type(tVec), pointer                       :: plasticStrain => null()
      type(tVec), pointer                       :: cumulatedPlasticDissipation => null()
      type(tVec), pointer                       :: stress => null()

      type(tPetscSF)                            :: displacementToIOSF, IOToDisplacementSF
      type(tPetscSF)                            :: cohesiveDisplacementToIOSF, IOToCohesiveDisplacementSF
      type(tPetscSF)                            :: displacementConstraintsSF
      type(tPetscSF)                            :: damageToIOSF, IOToDamageSF
      type(tPetscSF)                            :: damageConstraintsSF
      type(tPetscSF)                            :: temperatureToIOSF, IOToTemperatureSF
      type(tPetscSF)                            :: bodyForceToIOSF, IOToBodyForceSF
      type(tPetscSF)                            :: boundaryForceToIOSF, IOToBoundaryForceSF
      type(tPetscSF)                            :: pressureForceToIOSF, IOToPressureForceSF
      type(tPetscSF)                            :: stressToIOSF, IOToStressSF
      type(tPetscSF)                            :: plasticStrainToIOSF, IOToPlasticStrainSF
      type(tPetscSF)                            :: cumulatedPlasticDissToIOSF, IOToCumulatedPlasticDissSF

      type(MEF90DefMechGlobalOptions_Type)      :: globalOptions
      !!! Per-set options are not stored: they are read from the options database where they are needed,
      !!! with MEF90DefMech[Cell,Face,Vertex]SetOptionsSetFromOptions, the same way the AT model, the
      !!! energy split, and the Hooke's law are obtained. Only the number of sets is kept here.
      PetscInt                                  :: numCellSet = 0
      PetscInt                                  :: numFaceSet = 0
      PetscInt                                  :: numVertexSet = 0
      type(tPetscBag), dimension(:), pointer    :: MaterialPropertiesBag => null()

      type(tPetscViewer)                        :: globalEnergyViewer
      type(tPetscViewer), dimension(:), pointer :: setEnergyViewer => null()

      PetscBool                                 :: hasDisplacementBounds = PETSC_FALSE
      PetscBool                                 :: hasUnilateralContact = PETSC_FALSE

      !!! Handle on self, to be used as the application context of SNES, TAO, etc. Since MEF90DefMech_Type
      !!! extends MEF90Object, it has type-bound procedures and cannot be passed to the assumed-type (type(*))
      !!! context arguments of the PETSc Fortran API. PETSc callbacks recover the context with
      !!! call c_f_pointer(PETScCtx, MEF90DefMechCtx)
      type(c_ptr)                               :: PETScCtx = C_NULL_PTR
   contains
      procedure, pass(self) :: setFromOptions => MEF90DefMechSetFromOptions
      procedure, pass(self) :: view => MEF90DefMechView
   end type MEF90DefMech_Type

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCreate"
!!!
!!!
!!!  MEF90DefMechCreate: allocates all the fields of a MEF90DefMech_Type
!!!
!!!  (c) 2012-20 Blaise Bourdin bourdin@lsu.edu
!!!      2026    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechCreate(DefMech, dm, MEF90Ctx, prefix, ierr)
      !!! DefMech has the target attribute so that DefMech%PETScCtx can be handed over to PETSc as an
      !!! application context. The actual argument MUST therefore also have the target attribute.
      type(MEF90DefMech_Type), target, intent(OUT)              :: DefMech
      type(tDM), target, intent(IN)                             :: dm
      type(MEF90Ctx_Type), target, intent(IN)                   :: MEF90Ctx
      character(len=*), intent(IN)                              :: prefix
      PetscErrorCode, intent(INOUT)                             :: ierr

      type(MEF90CtxGlobalOptions_Type), pointer                 :: MEF90CtxGlobalOptions
      type(tIS)                                                 :: setIS
      PetscInt                                                  :: set
      PetscInt, dimension(:), pointer                           :: setID
      character(len=MEF90MXSTRLEN)                              :: filename, IOBuffer
      character(len=MEF90MXSTRLEN)                              :: vecName
      type(tDM), dimension(:), pointer                          :: dmList
      type(tPetscSF)                                            :: dummySF

      DefMech%MEF90Ctx => MEF90Ctx
      DefMech%comm = MEF90Ctx%comm
      DefMech%prefix = prefix
      DefMech%PETScCtx = c_loc(DefMech)

      !!!
      !!! Count the sets. The per-set options themselves are read on demand, see MEF90DefMech_Type
      !!!
      PetscCall(DMGetLabelIdIS(dm, MEF90CellSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetLocalSize(setIS, DefMech%numCellSet, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(dm, MEF90FaceSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetLocalSize(setIS, DefMech%numFaceSet, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(dm, MEF90VertexSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetLocalSize(setIS, DefMech%numVertexSet, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      !!!
      !!! Create energy viewers
      !!!
      filename = trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.ener'
      PetscCall(PetscViewerASCIIOpen(MEF90Ctx%comm, filename, DefMech%globalEnergyViewer, ierr))
      PetscCall(PetscViewerASCIIPrintf(DefMech%globalEnergyViewer, "# step     load            elastic energy  work            cohesive energy surface energy  total energy   plastic dissipation \n", ierr))
      PetscCall(PetscViewerFlush(DefMech%globalEnergyViewer, ierr))

      PetscCall(DMGetLabelIdIS(dm, MEF90CellSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))

      allocate (DefMech%setEnergyViewer(size(setID)), stat=ierr)
      do set = 1, size(setID)
         write (filename, 101) trim(MEF90FilePrefix(MEF90Ctx%resultFile)), setID(set)
         PetscCall(PetscViewerASCIIOpen(MEF90Ctx%comm, filename, DefMech%setEnergyViewer(set), ierr))
         write (IOBuffer, 102) setID(set)
         PetscCall(PetscViewerASCIIPrintf(DefMech%setEnergyViewer(set), IOBuffer, ierr))
         PetscCall(PetscViewerASCIIPrintf(DefMech%setEnergyViewer(set), "# step     load            elastic energy  work            cohesive energy surface energy  total energy   plastic dissipation\n", ierr))
         PetscCall(PetscViewerFlush(DefMech%setEnergyViewer(set), ierr))
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
101   format(A, '-', I4.4, '.enerblk')
102   format("# cell set ", I4, "\n")
      DefMech%analysisTime = 0.0_kr
      DefMech%timeStep = 0.0_kr

      !!! Create Vecs and SF
      PetscCall(DMGetDimension(dm, DefMech%dim, ierr))

      PetscCall(PetscBagGetDataMEF90CtxGlobalOptions(DefMech%MEF90Ctx%GlobalOptionsBag, MEF90CtxGlobalOptions, ierr))

      vecName = "Displacement"
      allocate (DefMech%displacementLocal, stat=ierr)
      PetscCall(MEF90CreateLocalVector(dm, MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, DefMech%dim, vecName, DefMech%displacementLocal, ierr))
      allocate (DefMech%displacementPreviousStepLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMech%displacementLocal, DefMech%displacementPreviousStepLocal, ierr))
      PetscCall(PetscObjectSetName(DefMech%displacementPreviousStepLocal, "DisplacementPreviousStep", ierr))
      allocate (DefMech%displacementLowerBoundLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMech%displacementLocal, DefMech%displacementLowerBoundLocal, ierr))
      PetscCall(PetscObjectSetName(DefMech%displacementLowerBoundLocal, "DisplacementLowerBound", ierr))
      allocate (DefMech%displacementUpperBoundLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMech%displacementLocal, DefMech%displacementUpperBoundLocal, ierr))
      PetscCall(PetscObjectSetName(DefMech%displacementUpperBoundLocal, "DisplacementUpperBound", ierr))

      vecName = "Damage"
      allocate (DefMech%damageLocal, stat=ierr)
      PetscCall(MEF90CreateLocalVector(dm, MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, 1_ki, vecName, DefMech%damageLocal, ierr))
      allocate (DefMech%damagePreviousStepLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMech%damageLocal, DefMech%damagePreviousStepLocal, ierr))
      PetscCall(PetscObjectSetName(DefMech%damagePreviousStepLocal, "damagePreviousStep", ierr))

      allocate (DefMech%TemperatureLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMech%damageLocal, DefMech%TemperatureLocal, ierr))
      PetscCall(PetscObjectSetName(DefMech%TemperatureLocal, "Temperature", ierr))

      vecName = "cohesiveDisplacement"
      allocate (DefMech%cohesiveDisplacement, stat=ierr)
      PetscCall(MEF90CreateLocalVector(dm, MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, DefMech%dim, vecName, DefMech%cohesiveDisplacement, ierr))
      vecName = "bodyForce"
      allocate (DefMech%bodyForce, stat=ierr)
      PetscCall(MEF90CreateCellVector(dm, DefMech%dim, vecName, DefMech%bodyForce, ierr))
      vecName = "boundaryForce"
      allocate (DefMech%boundaryForce, stat=ierr)
      PetscCall(MEF90CreateBoundaryCellVector(dm, DefMech%dim, vecName, DefMech%boundaryForce, ierr))
      vecName = "pressureForce"
      allocate (DefMech%pressureForce, stat=ierr)
      PetscCall(MEF90CreateBoundaryCellVector(dm, 1_ki, vecName, DefMech%pressureForce, ierr))

      vecName = "plasticStrain"
      allocate (DefMech%plasticStrain, stat=ierr)
      PetscCall(MEF90CreateCellVector(dm, (DefMech%dim * (DefMech%dim + 1_ki)) / 2_ki, vecName, DefMech%plasticStrain, ierr))
      allocate (DefMech%cumulatedPlasticDissipation, stat=ierr)
      PetscCall(VecDuplicate(DefMech%plasticStrain, DefMech%cumulatedPlasticDissipation, ierr))
      PetscCall(PetscObjectSetName(DefMech%cumulatedPlasticDissipation, "cumulatedPlasticDissipation", ierr))
      allocate (DefMech%stress, stat=ierr)
      PetscCall(VecDuplicate(DefMech%plasticStrain, DefMech%stress, ierr))
      PetscCall(PetscObjectSetName(DefMech%stress, "Stress", ierr))

      !!! Create megaDM
      allocate (dmList(7))
      PetscCall(VecGetDM(DefMech%displacementLocal, dmList(1), ierr))
      PetscCall(VecGetDM(DefMech%damageLocal, dmList(2), ierr))
      PetscCall(VecGetDM(DefMech%temperatureLocal, dmList(3), ierr))
      PetscCall(VecGetDM(DefMech%bodyForce, dmList(4), ierr))
      PetscCall(VecGetDM(DefMech%boundaryForce, dmList(5), ierr))
      PetscCall(VecGetDM(DefMech%pressureForce, dmList(6), ierr))
      PetscCall(VecGetDM(DefMech%plasticStrain, dmList(7), ierr))
      PetscCall(DMCreateSuperDM(dmList, 7_ki, PETSC_NULL_IS_POINTER, DefMech%megaDM, ierr))
      deallocate (dmList)

      !!! Create the IO SF for all fields
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMech%displacementLocal, DefMech%displacementToIOSF, DefMech%IOTodisplacementSF, ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMech%damageLocal, DefMech%damageToIOSF, DefMech%IOTodamageSF, ierr))

      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMech%cohesiveDisplacement, DefMech%cohesiveDisplacementToIOSF, DefMech%IOToCohesiveDisplacementSF, ierr))

      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMech%temperatureLocal, DefMech%temperatureToIOSF, DefMech%IOTotemperatureSF, ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMech%bodyForce, DefMech%bodyForceToIOSF, DefMech%IOTobodyForceSF, ierr))
      ! PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMech%boundaryForce,DefMech%boundaryForceToIOSF,DefMech%IOToboundaryForceSF,ierr))
      ! PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMech%pressureForce,DefMech%pressureForceToIOSF,DefMech%IOTopressureForceSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMech%stress, DefMech%stressToIOSF, DefMech%IOToStressSF, ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMech%plasticStrain, DefMech%plasticStrainToIOSF, DefMech%IOToplasticStrainSF, ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMech%cumulatedPlasticDissipation, DefMech%cumulatedPlasticDissToIOSF, DefMech%IOToCumulatedPlasticDissSF, ierr))

      !!! Create the SF to exchange boundary values of the displacement and damage.
      PetscCall(MEF90ConstraintSFCreate(DefMech%MEF90Ctx, DefMech%displacementLocal, DefMech%displacementLocal, DefMech%displacementConstraintsSF, dummySF, ierr))
      PetscCall(PetscSFDestroy(dummySF, ierr))
      PetscCall(MEF90ConstraintSFCreate(DefMech%MEF90Ctx, DefMech%damageLocal, DefMech%damageLocal, DefMech%damageConstraintsSF, dummySF, ierr))
      PetscCall(PetscSFDestroy(dummySF, ierr))
   end subroutine MEF90DefMechCreate

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechDestroy"
!!!
!!!
!!!  MEF90DefMechDestroy: destroys a MEF90DefMech_Type
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2026    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechDestroy(DefMech, ierr)
      type(MEF90DefMech_Type), intent(INOUT)           :: DefMech
      PetscErrorCode, intent(INOUT)                    :: ierr

      PetscInt                                        :: set

      DefMech%PETScCtx = C_NULL_PTR

      !!!
      !!! Close energy viewers
      !!!
      PetscCall(PetscViewerDestroy(DefMech%globalEnergyViewer, ierr))
      do set = 1, size(DefMech%setEnergyViewer)
         PetscCall(PetscViewerDestroy(DefMech%setEnergyViewer(set), ierr))
      end do
      deallocate (DefMech%setEnergyViewer)

      !!! Destroy Vecs and SF and deAllocate them
      if (associated(DefMech%displacementLocal)) then
         PetscCall(VecDestroy(DefMech%displacementLocal, ierr))
         deallocate (DefMech%displacementLocal)
         nullify (DefMech%displacementLocal)
      end if
      if (associated(DefMech%displacementPreviousStepLocal)) then
         PetscCall(VecDestroy(DefMech%displacementPreviousStepLocal, ierr))
         deallocate (DefMech%displacementPreviousStepLocal)
         nullify (DefMech%displacementPreviousStepLocal)
      end if
      if (associated(DefMech%displacementLowerBoundLocal)) then
         PetscCall(VecDestroy(DefMech%displacementLowerBoundLocal, ierr))
         deallocate (DefMech%displacementLowerBoundLocal)
         nullify (DefMech%displacementLowerBoundLocal)
      end if
      if (associated(DefMech%displacementUpperBoundLocal)) then
         PetscCall(VecDestroy(DefMech%displacementUpperBoundLocal, ierr))
         deallocate (DefMech%displacementUpperBoundLocal)
         nullify (DefMech%displacementUpperBoundLocal)
      end if

      if (associated(DefMech%damageLocal)) then
         PetscCall(VecDestroy(DefMech%damageLocal, ierr))
         deallocate (DefMech%damageLocal)
         nullify (DefMech%damageLocal)
      end if
      if (associated(DefMech%damagePreviousStepLocal)) then
         PetscCall(VecDestroy(DefMech%damagePreviousStepLocal, ierr))
         deallocate (DefMech%damagePreviousStepLocal)
         nullify (DefMech%damagePreviousStepLocal)
      end if

      if (associated(DefMech%temperatureLocal)) then
         PetscCall(VecDestroy(DefMech%temperatureLocal, ierr))
         deallocate (DefMech%temperatureLocal)
         nullify (DefMech%temperatureLocal)
      end if

      if (associated(DefMech%cohesiveDisplacement)) then
         PetscCall(VecDestroy(DefMech%cohesiveDisplacement, ierr))
         deallocate (DefMech%cohesiveDisplacement)
         nullify (DefMech%cohesiveDisplacement)
      end if

      if (associated(DefMech%bodyForce)) then
         PetscCall(VecDestroy(DefMech%bodyForce, ierr))
         deallocate (DefMech%bodyForce)
         nullify (DefMech%bodyForce)
      end if
      if (associated(DefMech%boundaryForce)) then
         PetscCall(VecDestroy(DefMech%boundaryForce, ierr))
         deallocate (DefMech%boundaryForce)
         nullify (DefMech%boundaryForce)
      end if
      if (associated(DefMech%pressureForce)) then
         PetscCall(VecDestroy(DefMech%pressureForce, ierr))
         deallocate (DefMech%pressureForce)
         nullify (DefMech%pressureForce)
      end if

      if (associated(DefMech%plasticStrain)) then
         PetscCall(VecDestroy(DefMech%plasticStrain, ierr))
         deallocate (DefMech%plasticStrain)
         nullify (DefMech%plasticStrain)
      end if
      if (associated(DefMech%cumulatedPlasticDissipation)) then
         PetscCall(VecDestroy(DefMech%cumulatedPlasticDissipation, ierr))
         deallocate (DefMech%cumulatedPlasticDissipation)
         nullify (DefMech%cumulatedPlasticDissipation)
      end if
      if (associated(DefMech%Stress)) then
         PetscCall(VecDestroy(DefMech%Stress, ierr))
         deallocate (DefMech%Stress)
         nullify (DefMech%Stress)
      end if

      !!! Destroy all PetscSF
      PetscCall(PetscSFDestroy(DefMech%displacementToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%IOTodisplacementSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%cohesiveDisplacementToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%IOToCohesiveDisplacementSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%damageToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%IOTodamageSF, ierr))

      PetscCall(PetscSFDestroy(DefMech%temperatureToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%IOTotemperatureSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%bodyForceToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%IOTobodyForceSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%boundaryForceToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%IOToboundaryForceSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%pressureForceToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%IOTopressureForceSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%stressToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%IOToStressSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%plasticStrainToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%IOToplasticStrainSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%cumulatedPlasticDissToIOSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%IOToCumulatedPlasticDissSF, ierr))

      !!! Destroy the SF to exchange boundary values of the displacement and damage.
      PetscCall(PetscSFDestroy(DefMech%displacementConstraintsSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%damageConstraintsSF, ierr))

      !!! Destroy the megaDM
      PetscCall(DMDestroy(DefMech%megaDM, ierr))
   end subroutine MEF90DefMechDestroy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechGlobalOptionsSetFromOptions_Private"
!!!
!!!
!!!  MEF90DefMechGlobalOptionsSetFromOptions_Private: reads the problem-wide options of a MEF90DefMech_Type
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechGlobalOptionsSetFromOptions_Private(comm, prefix, options, ierr)
      MPIU_Comm, intent(IN)                                  :: comm
      character(len=*), intent(IN)                           :: prefix
      type(MEF90DefMechGlobalOptions_Type), intent(INOUT)    :: options
      PetscErrorCode, intent(INOUT)                          :: ierr

      PetscCall(PetscOptionsBegin(comm, prefix, "Options for MEF90DefMech_Type", "mef90DefMech", ierr))
         PetscCall(PetscOptionsEnum('-DefMech_TimeStepping_Type', 'Type of defect mechanics Time steping', 'mef90DefMech', MEF90DefMech_TimeSteppingTypeList, options%timeSteppingType, options%timeSteppingType, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-DefMech_solver_Type', 'Type of defect mechanics solver', 'mef90DefMech', MEF90DefMech_SolverTypeList, options%solverType, options%solverType, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-DefMech_damageSolver_Type', 'Type of defect mechanics damage solver', 'mef90DefMech', MEF90DefMech_DamageSolverTypeList, options%damageSolverType, options%damageSolverType, PETSC_NULL_BOOL, ierr))

         PetscCall(PetscOptionsBool('-temperature_export', 'Export temperature in result file', 'mef90DefMech', options%temperatureExport, options%temperatureExport, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-displacement_export', 'Export displacement in result file', 'mef90DefMech', options%displacementExport, options%displacementExport, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-damage_export', 'Export damage in result file', 'mef90DefMech', options%damageExport, options%damageExport, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-stress_export', 'Export stress in result file', 'mef90DefMech', options%stressExport, options%stressExport, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-plasticstrain_export', 'Export plastic strain in result file', 'mef90DefMech', options%plasticStrainExport, options%plasticStrainExport, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-cumulatedplasticdissipation_export', 'Export cumulated plastic dissipation in result file', 'mef90DefMech', options%cumulatedPlasticDissipationExport, options%cumulatedPlasticDissipationExport, PETSC_NULL_BOOL, ierr))

         PetscCall(PetscOptionsEnum('-boundaryDisplacement_scaling', 'Boundary displacement scaling', 'mef90DefMech', MEF90ScalingList, options%boundaryDisplacementScaling, options%boundaryDisplacementScaling, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-displacementlowerbound_scaling', 'Displacement lower bound scaling', 'mef90DefMech', MEF90ScalingList, options%displacementLowerBoundScaling, options%displacementLowerBoundScaling, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-displacementupperbound_scaling', 'Displacement upper bound scaling', 'mef90DefMech', MEF90ScalingList, options%displacementUpperBoundScaling, options%displacementUpperBoundScaling, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-cohesiveDisplacement_scaling', 'Cohesive displacement scaling', 'mef90DefMech', MEF90ScalingList, options%cohesiveDisplacementScaling, options%cohesiveDisplacementScaling, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-boundaryDamage_scaling', 'Boundary damage scaling', 'mef90DefMech', MEF90ScalingList, options%boundaryDamageScaling, options%boundaryDamageScaling, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-bodyforce_scaling', 'Body force scaling', 'mef90DefMech', MEF90ScalingList, options%bodyForceScaling, options%bodyForceScaling, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-boundaryforce_scaling', 'Boundary force scaling', 'mef90DefMech', MEF90ScalingList, options%boundaryForceScaling, options%boundaryForceScaling, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-pressureForce_scaling', 'Pressure force scaling', 'mef90DefMech', MEF90ScalingList, options%pressureForceScaling, options%pressureForceScaling, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-crackPressure_scaling', 'Crack Pressure scaling', 'mef90DefMech', MEF90ScalingList, options%CrackPressureScaling, options%CrackPressureScaling, PETSC_NULL_BOOL, ierr))

         PetscCall(PetscOptionsReal('-defmech_damage_atol', 'Absolute tolerance on damage error', 'mef90DefMech', options%damageATol, options%damageATol, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-defmech_damage_maxit', 'Maximum number of alternate minimizations for damage', 'mef90DefMech', options%damageMaxIt, options%damageMaxIt, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-defmech_pclag', 'Interval at which the PC is recomputed during alternate minimization', 'mef90DefMech', options%PCLag, options%PCLag, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-defmech_SOR_Omega', 'Alterate Minimization over relaxation factor (>0 for limited, <0 for projected)', 'mef90DefMech', options%SOROmega, options%SOROmega, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-defmech_irrevThres', 'Threshold above which irreversibility is enforced (0 for monotonicity, .99 for equality)', 'mef90DefMech', options%irrevthres, options%irrevthres, PETSC_NULL_BOOL, ierr))

         PetscCall(PetscOptionsEnum('-BT_Type', 'Backtracking type', 'mef90DefMech', MEF90DefMech_BTTypeList, options%BTType, options%BTType, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-BT_Interval', 'Interval at which Backtracking is run in inner loop (0 for outer loop)', 'mef90DefMech', options%BTInterval, options%BTInterval, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-BT_Scope', 'Backtracking scope (0 for unlimited)', 'mef90DefMech', options%BTScope, options%BTScope, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-BT_Tol', 'Backtracking relative tolerance', 'mef90DefMech', options%BTTol, options%BTTol, PETSC_NULL_BOOL, ierr))

         PetscCall(PetscOptionsReal('-defmech_plasticstrain_atol', 'Absolute tolerance on plastic error', 'mef90DefMech', options%plasticStrainATol, options%plasticStrainATol, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-defmech_InjectedVolume_atol', 'Absolute tolerance on injected volume error', 'mef90DefMech', options%InjectedVolumeATol, options%InjectedVolumeATol, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-defmech_dampingCoefficient_displacement', 'Damping coefficient on displacement field (0 for minimization, 1 for semi-implicit gradient flow)', 'mef90DefMech', options%dampingCoefficientDisplacement, options%dampingCoefficientDisplacement, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-defmech_dampingCoefficient_damage', 'Damping coefficient on damage field (0 for minimization, 1 for semi-implicit gradient flow)', 'mef90DefMech', options%dampingCoefficientDamage, options%dampingCoefficientDamage, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))
   end subroutine MEF90DefMechGlobalOptionsSetFromOptions_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCellSetOptionsSetFromOptions"
!!!
!!!
!!!  MEF90DefMechCellSetOptionsSetFromOptions: reads the options of a single cell set
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechCellSetOptionsSetFromOptions(comm, prefix, options, ierr)
      MPIU_Comm, intent(IN)                                  :: comm
      character(len=*), intent(IN)                           :: prefix
      type(MEF90DefMechCellSetOptions_Type), intent(OUT)     :: options
      PetscErrorCode, intent(INOUT)                          :: ierr

      PetscInt                                               :: nOpt

      PetscCall(PetscOptionsBegin(comm, prefix, "Options for a MEF90DefMech_Type cell set", "mef90DefMech", ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-bodyForce', '[N.m^(-3) / N.m^(-2)] (f): body force', 'mef90DefMech', options%bodyForce, nOpt, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-CrackPressure', 'without unit: internal crack pressure', 'mef90DefMech', options%CrackPressure, options%CrackPressure, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-plasticity_type', 'Type of plasticity law', 'mef90DefMech', MEF90DefMech_plasticityTypeList, options%plasticityType, options%plasticityType, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-cohesiveDisplacement', '[m] (U): Cohesive displacement value', 'mef90DefMech', options%cohesiveDisplacement, nOpt, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsBoolArray('-DisplacementBC', 'Displacement has Dirichlet boundary Condition (Y/N)', 'mef90DefMech', options%Has_displacementBC, nOpt, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-boundaryDisplacement', '[m] (U): Displacement boundary value', 'mef90DefMech', options%boundaryDisplacement, nOpt, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-displacementLowerBound', '[m] (U): Displacement lower bound', 'mef90DefMech', options%displacementLowerBound, nOpt, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-displacementUpperBound', '[m] (U): Displacement upper bound', 'mef90DefMech', options%displacementUpperBound, nOpt, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-DamageBC', 'Damage has Dirichlet boundary Condition (Y/N)', 'mef90DefMech', options%Has_DamageBC, options%Has_DamageBC, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-CrackVolumeControlled', 'Crack Pressure controlled by the crack volume in this block (Y/N)', 'mef90DefMech', options%CrackVolumeControlled, options%CrackVolumeControlled, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-WorkControlled', 'Force magnitude controlled by its work in this block (Y/N)', 'mef90DefMech', options%WorkControlled, options%WorkControlled, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-boundaryDamage', '[unit-less] (alpha): Damage boundary value', 'mef90DefMech', options%boundaryDamage, options%boundaryDamage, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))
   end subroutine MEF90DefMechCellSetOptionsSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechFaceSetOptionsSetFromOptions"
!!!
!!!
!!!  MEF90DefMechFaceSetOptionsSetFromOptions: reads the options of a single face set
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechFaceSetOptionsSetFromOptions(comm, prefix, options, ierr)
      MPIU_Comm, intent(IN)                                  :: comm
      character(len=*), intent(IN)                           :: prefix
      type(MEF90DefMechFaceSetOptions_Type), intent(OUT)     :: options
      PetscErrorCode, intent(INOUT)                          :: ierr

      PetscInt                                               :: nOpt

      PetscCall(PetscOptionsBegin(comm, prefix, "Options for a MEF90DefMech_Type face set", "mef90DefMech", ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-boundaryForce', '[N.m^(-2) / N.m^(-1)] (f): boundary force', 'mef90DefMech', options%boundaryForce, nOpt, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-pressureForce', 'without unit: internal crack pressure', 'mef90DefMech', options%pressureForce, options%pressureForce, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsBoolArray('-DisplacementBC', 'Displacement has Dirichlet boundary Condition (Y/N)', 'mef90DefMech', options%Has_displacementBC, nOpt, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-boundaryDisplacement', '[m] (U): Displacement boundary value', 'mef90DefMech', options%boundaryDisplacement, nOpt, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-displacementLowerBound', '[m] (U): Displacement lower bound', 'mef90DefMech', options%displacementLowerBound, nOpt, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-displacementUpperBound', '[m] (U): Displacement upper bound', 'mef90DefMech', options%displacementUpperBound, nOpt, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-DamageBC', 'Damage has Dirichlet boundary Condition (Y/N)', 'mef90DefMech', options%Has_DamageBC, options%Has_DamageBC, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-boundaryDamage', '[unit-less] (alpha): Damage boundary value', 'mef90DefMech', options%boundaryDamage, options%boundaryDamage, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))
   end subroutine MEF90DefMechFaceSetOptionsSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechVertexSetOptionsSetFromOptions"
!!!
!!!
!!!  MEF90DefMechVertexSetOptionsSetFromOptions: reads the options of a single vertex set
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechVertexSetOptionsSetFromOptions(comm, prefix, options, ierr)
      MPIU_Comm, intent(IN)                                  :: comm
      character(len=*), intent(IN)                           :: prefix
      type(MEF90DefMechVertexSetOptions_Type), intent(OUT)   :: options
      PetscErrorCode, intent(INOUT)                          :: ierr

      PetscInt                                               :: nOpt

      PetscCall(PetscOptionsBegin(comm, prefix, "Options for a MEF90DefMech_Type vertex set", "mef90DefMech", ierr))
         nOpt = 3
         PetscCall(PetscOptionsBoolArray('-DisplacementBC', 'Displacement has Dirichlet boundary Condition (Y/N)', 'mef90DefMech', options%Has_displacementBC, nOpt, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-boundaryDisplacement', '[m] (U): Displacement boundary value', 'mef90DefMech', options%boundaryDisplacement, nOpt, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-displacementLowerBound', '[m] (U): Displacement lower bound', 'mef90DefMech', options%displacementLowerBound, nOpt, PETSC_NULL_BOOL, ierr))
         nOpt = 3
         PetscCall(PetscOptionsRealArray('-displacementUpperBound', '[m] (U): Displacement upper bound', 'mef90DefMech', options%displacementUpperBound, nOpt, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-DamageBC', 'Damage has Dirichlet boundary Condition (Y/N)', 'mef90DefMech', options%Has_DamageBC, options%Has_DamageBC, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-boundaryDamage', '[unit-less] (alpha): boundaryDamage', 'mef90DefMech', options%boundaryDamage, options%boundaryDamage, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))
   end subroutine MEF90DefMechVertexSetOptionsSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechSetFromOptions"
!!!
!!!
!!!  MEF90DefMechSetFromOptions: initializes a MEF90DefMech_Type from options
!!!
!!!  (c) 2012-14 Blaise Bourdin bourdin@lsu.edu
!!!      2026    Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechSetFromOptions(self, ierr)
      class(MEF90DefMech_Type), intent(INOUT)                :: self
      PetscErrorCode, intent(INOUT)                          :: ierr

      type(tIS)                                              :: setIS
      PetscInt, dimension(:), pointer                        :: setID
      PetscInt                                               :: set
      character(len=MEF90MXSTRLEN)                           :: setPrefix
      type(MEF90DefMechCellSetOptions_Type)                  :: cellSetOptions
      type(MEF90DefMechFaceSetOptions_Type)                  :: faceSetOptions
      PetscBool                                              :: printHelp

      !!!
      !!! Problem-wide options
      !!!
      PetscCall(MEF90DefMechGlobalOptionsSetFromOptions_Private(self%comm, trim(self%prefix), self%globalOptions, ierr))

      !!!
      !!! The per-set options are read where they are needed, but they have to be visited once here so
      !!! that they are registered with the options database, and so that the displacement bounds of all
      !!! sets are known before the solvers are set up
      !!!
      self%hasDisplacementBounds = PETSC_FALSE

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90CellSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'cs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90DefMechCellSetOptionsSetFromOptions(self%comm, trim(setPrefix), cellSetOptions, ierr))
         self%hasDisplacementBounds = any(cellSetOptions%displacementLowerBound /= MEF90NINFINITY) .or. self%hasDisplacementBounds
         self%hasDisplacementBounds = any(cellSetOptions%displacementUpperBound /= MEF90INFINITY) .or. self%hasDisplacementBounds
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90FaceSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'fs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90DefMechFaceSetOptionsSetFromOptions(self%comm, trim(setPrefix), faceSetOptions, ierr))
         self%hasDisplacementBounds = any(faceSetOptions%displacementLowerBound /= MEF90NINFINITY) .or. self%hasDisplacementBounds
         self%hasDisplacementBounds = any(faceSetOptions%displacementUpperBound /= MEF90INFINITY) .or. self%hasDisplacementBounds
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", printHelp, PETSC_NULL_BOOL, ierr))
      if (printHelp) then
         call self%view(PETSC_VIEWER_STDOUT_WORLD, ierr)
      end if
   end subroutine MEF90DefMechSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechView"
!!!
!!!
!!!  MEF90DefMechView: the default viewer for a MEF90DefMech_Type
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechView(self, viewer, ierr)
      class(MEF90DefMech_Type), intent(IN)                   :: self
      type(tPetscViewer), intent(IN)                         :: viewer
      PetscErrorCode, intent(INOUT)                          :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)              :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)              :: viewerType
      type(tIS)                                              :: setIS
      PetscInt, dimension(:), pointer                        :: setID
      PetscInt                                               :: set
      character(len=MEF90MXSTRLEN)                           :: setPrefix
      type(MEF90DefMechCellSetOptions_Type)                  :: cellSetOptions
      type(MEF90DefMechFaceSetOptions_Type)                  :: faceSetOptions
      type(MEF90DefMechVertexSetOptions_Type)                :: vertexSetOptions

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType /= 'ascii') return

      write (IOBuffer, "(A,': Options for MEF90DefMech_Type\n')") trim(self%prefix)//"defmech"
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         time stepping type: ',A,'\n')") trim(MEF90DefMech_TimeSteppingTypeList(self%globalOptions%timeSteppingType + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         solver type: ',A,'\n')") trim(MEF90DefMech_SolverTypeList(self%globalOptions%solverType + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         damage solver type: ',A,'\n')") trim(MEF90DefMech_DamageSolverTypeList(self%globalOptions%damageSolverType + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         damage atol / maxit: ',ES12.5,' / ',I6,'\n')") self%globalOptions%damageATol, self%globalOptions%damageMaxIt
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         SOR omega: ',ES12.5,' irreversibility threshold: ',ES12.5,'\n')") self%globalOptions%SOROmega, self%globalOptions%irrevthres
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         backtracking type: ',A,' interval: ',I6,' scope: ',I6,' tol: ',ES12.5,'\n')") &
         trim(MEF90DefMech_BTTypeList(self%globalOptions%BTType + 1)), self%globalOptions%BTInterval, self%globalOptions%BTScope, self%globalOptions%BTTol
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         export displacement / damage / stress: ',3(L1,' '),'\n')") &
         self%globalOptions%displacementExport, self%globalOptions%damageExport, self%globalOptions%stressExport
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         export temperature / plastic strain / cumulated plastic dissipation: ',3(L1,' '),'\n')") &
         self%globalOptions%temperatureExport, self%globalOptions%plasticStrainExport, self%globalOptions%cumulatedPlasticDissipationExport
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90CellSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'cs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90DefMechCellSetOptionsSetFromOptions(self%comm, trim(setPrefix), cellSetOptions, ierr))
         write (IOBuffer, "(A,'cs',I4.4,': body force: ',3(ES12.5,' '),'plasticity: ',A,'\n')") &
            trim(self%prefix), setID(set), cellSetOptions%bodyForce, trim(MEF90DefMech_plasticityTypeList(cellSetOptions%plasticityType + 1))
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
         write (IOBuffer, "('         displacement BC: ',3(L1,' '),'value: ',3(ES12.5,' '),'\n')") &
            cellSetOptions%Has_displacementBC, cellSetOptions%boundaryDisplacement
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
         write (IOBuffer, "('         damage BC: ',L1,' value: ',ES12.5,'\n')") &
            cellSetOptions%Has_damageBC, cellSetOptions%boundaryDamage
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90FaceSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'fs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90DefMechFaceSetOptionsSetFromOptions(self%comm, trim(setPrefix), faceSetOptions, ierr))
         write (IOBuffer, "(A,'fs',I4.4,': boundary force: ',3(ES12.5,' '),'pressure force: ',ES12.5,'\n')") &
            trim(self%prefix), setID(set), faceSetOptions%boundaryForce, faceSetOptions%pressureForce
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
         write (IOBuffer, "('         displacement BC: ',3(L1,' '),'value: ',3(ES12.5,' '),'\n')") &
            faceSetOptions%Has_displacementBC, faceSetOptions%boundaryDisplacement
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
         write (IOBuffer, "('         damage BC: ',L1,' value: ',ES12.5,'\n')") &
            faceSetOptions%Has_damageBC, faceSetOptions%boundaryDamage
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90VertexSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'vs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90DefMechVertexSetOptionsSetFromOptions(self%comm, trim(setPrefix), vertexSetOptions, ierr))
         write (IOBuffer, "(A,'vs',I4.4,': displacement BC: ',3(L1,' '),'value: ',3(ES12.5,' '),'\n')") &
            trim(self%prefix), setID(set), vertexSetOptions%Has_displacementBC, vertexSetOptions%boundaryDisplacement
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
         write (IOBuffer, "('         damage BC: ',L1,' value: ',ES12.5,'\n')") &
            vertexSetOptions%Has_damageBC, vertexSetOptions%boundaryDamage
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end subroutine MEF90DefMechView
end module m_MEF90_DefMech_class
