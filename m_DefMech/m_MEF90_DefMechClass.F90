#include "../MEF90/mef90.inc"
module m_MEF90_DefMech_class
#include "petsc/finclude/petsc.h"
   use petscsys
   use m_MEF90_BaseClass
   use m_MEF90_Ctx
   use m_MEF90_DMPlex
   use m_MEF90_LinAlg
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
         MEF90DefMech_unilateralContactTypeDeviatoric
   end enum
   character(len=MEF90MXSTRLEN), dimension(7), protected   :: MEF90DefMech_unilateralContactTypeList = [character(len=MEF90MXSTRLEN) :: &
      'None', &
      'HydrostaticDeviatoric', &
      'Hydrostatic', &
      'Deviatoric', &
      'MEF90DefMech_unilateralContactTypeList', &
      '_MEF90DefMech_unilateralContactTypeList', &
      '']

!!!
!!!  MEF90DefMechGlobalOptions_Type: the problem-wide options of a MEF90DefMech_Type.
!!!                                  The values given here are the defaults used by setFromOptions
!!!
   type, extends(MEF90Object) :: MEF90DefMechGlobalOptions_Type
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
      PetscBool                              :: multiPhaseField = PETSC_FALSE
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
   contains
      procedure, pass(self) :: view_internal => MEF90DefMechGlobalOptionsView
   end type MEF90DefMechGlobalOptions_Type

   type, extends(MEF90Object) :: MEF90DefMechCellSetOptions_Type
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
      !!! Material properties. linearThermalExpansion is a symmetric matrix, so it is dimension
      !!! dependent and allocated by MEF90DefMechCellSetOptionsSetFromOptions
      class(mef90Mat), allocatable            :: linearThermalExpansion
      PetscReal                               :: cohesiveStiffness = 0.0_kr
   contains
      procedure, pass(self) :: view_internal => MEF90DefMechCellSetOptionsView
   end type MEF90DefMechCellSetOptions_Type

   type, extends(MEF90Object) :: MEF90DefMechFaceSetOptions_Type
      PetscReal, dimension(3)                 :: boundaryforce = [0.0_kr, 0.0_kr, 0.0_kr]
      PetscReal                               :: pressureForce = 0.0_kr
      PetscBool, dimension(3)                 :: Has_displacementBC = [PETSC_FALSE, PETSC_FALSE, PETSC_FALSE]
      PetscReal, dimension(3)                 :: boundaryDisplacement = [0.0_kr, 0.0_kr, 0.0_kr]
      PetscReal, dimension(3)                 :: displacementLowerBound = [MEF90NINFINITY, MEF90NINFINITY, MEF90NINFINITY]
      PetscReal, dimension(3)                 :: displacementUpperBound = [MEF90INFINITY, MEF90INFINITY, MEF90INFINITY]
      PetscBool                               :: Has_damageBC = PETSC_FALSE
      PetscReal                               :: boundaryDamage = 0.0_kr
   contains
      procedure, pass(self) :: view_internal => MEF90DefMechFaceSetOptionsView
   end type MEF90DefMechFaceSetOptions_Type

   type, extends(MEF90Object) :: MEF90DefMechVertexSetOptions_Type
      PetscBool, dimension(3)                 :: Has_displacementBC = [PETSC_FALSE, PETSC_FALSE, PETSC_FALSE]
      PetscReal, dimension(3)                 :: boundaryDisplacement = [0.0_kr, 0.0_kr, 0.0_kr]
      PetscReal, dimension(3)                 :: displacementLowerBound = [MEF90NINFINITY, MEF90NINFINITY, MEF90NINFINITY]
      PetscReal, dimension(3)                 :: displacementUpperBound = [MEF90INFINITY, MEF90INFINITY, MEF90INFINITY]
      PetscBool                               :: Has_damageBC = PETSC_FALSE
      PetscReal                               :: boundaryDamage = 0.0_kr
   contains
      procedure, pass(self) :: view_internal => MEF90DefMechVertexSetOptionsView
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
      type(tVec), dimension(:), pointer         :: partialDamageLocal => null()
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

      !!! Neither the problem-wide nor the per-set options are stored: they are read from the options
      !!! database where they are needed, with MEF90DefMechGlobalOptionsSetFromOptions and
      !!! MEF90DefMech[Cell,Face,Vertex]SetOptionsSetFromOptions, the same way the AT model, the
      !!! energy split, and the Hooke's law are obtained. The number of sets is obtained from the
      !!! megaDM with MEF90DMGetNumSets.

      type(tPetscViewer)                        :: globalEnergyViewer
      type(tPetscViewer), dimension(:), pointer :: setEnergyViewer => null()

      PetscBool                                 :: hasDisplacementBounds = PETSC_FALSE
      PetscBool                                 :: hasUnilateralContact = PETSC_FALSE

      PetscInt                               :: currentSet = 0_Ki ! used to pass the PF number in multiPhaseField
      !!! Handle on self, set once in MEF90DefMechCreate with PETScCtx = c_loc(DefMech), and handed to
      !!! PETSc wherever an application context is expected: SNESSetFunction, SNESSetJacobian,
      !!! TAOSetObjective, ... The callbacks in m_MEF90_DefMech recover the context with
      !!!    call c_f_pointer(PETScCtx, MEF90DefMechCtx)
      !!!
      !!! MEF90DefMech_Type cannot be passed to PETSc directly, the way MEF90DefMechCtx_Type was up to
      !!! mef90 0.5.2, because extending MEF90Object gives it type-bound procedures and PETSc declares
      !!! its context arguments as assumed-type. This has to be a component rather than a local or an
      !!! inline c_loc(), for lifetime reasons that are easy to get wrong: see the note at the top of
      !!! MEF90/m_MEF90_BaseClass.F90 before changing any of this.
      type(c_ptr)                               :: PETScCtx = C_NULL_PTR
   contains
      procedure, pass(self) :: setFromOptions => MEF90DefMechSetFromOptions
      procedure, pass(self) :: view_internal => MEF90DefMechView
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

   subroutine MEF90DefMechCreate(DefMechCtx, dm, MEF90Ctx, prefix, ierr)
      !!! DefMech has the target attribute so that DefMech%PETScCtx can be handed over to PETSc as an
      !!! application context. The actual argument MUST therefore also have the target attribute.
      type(MEF90DefMech_Type), target, intent(OUT)              :: DefMechCtx
      type(tDM), target, intent(IN)                             :: dm
      type(MEF90Ctx_Type), target, intent(IN)                   :: MEF90Ctx
      character(len=*), intent(IN)                              :: prefix
      PetscErrorCode, intent(INOUT)                             :: ierr

      type(MEF90CtxGlobalOptions_Type)                          :: MEF90CtxGlobalOptions
      type(MEF90DefMechGlobalOptions_Type)                      :: DefMechGlobalOptions
      type(tIS)                                                 :: setIS
      PetscInt                                                  :: set
      PetscInt, dimension(:), pointer                           :: setID
      character(len=MEF90MXSTRLEN)                              :: filename, IOBuffer
      character(len=MEF90MXSTRLEN)                              :: vecName
      type(tDM), dimension(:), pointer                          :: dmList
      type(tPetscSF)                                            :: dummySF

      DefMechCtx%MEF90Ctx => MEF90Ctx
      DefMechCtx%comm = MEF90Ctx%comm
      DefMechCtx%prefix = prefix
      DefMechCtx%name = trim(prefix)//"DefMech"
      DefMechCtx%PETScCtx = c_loc(DefMechCtx)

      !!
      !! Create energy viewers
      !!
      filename = trim(MEF90FilePrefix(MEF90Ctx%resultFile))//'.ener'
      PetscCall(PetscViewerASCIIOpen(MEF90Ctx%comm, filename, DefMechCtx%globalEnergyViewer, ierr))
      PetscCall(PetscViewerASCIIPrintf(DefMechCtx%globalEnergyViewer, "# step     load            elastic energy  work            cohesive energy surface energy  total energy   plastic dissipation \n", ierr))
      PetscCall(PetscViewerFlush(DefMechCtx%globalEnergyViewer, ierr))

      PetscCall(DMGetLabelIdIS(dm, MEF90CellSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(MEF90Ctx%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))

      allocate (DefMechCtx%setEnergyViewer(size(setID)), stat=ierr)
      do set = 1, size(setID)
         write (filename, 101) trim(MEF90FilePrefix(MEF90Ctx%resultFile)), setID(set)
         PetscCall(PetscViewerASCIIOpen(MEF90Ctx%comm, filename, DefMechCtx%setEnergyViewer(set), ierr))
         write (IOBuffer, 102) setID(set)
         PetscCall(PetscViewerASCIIPrintf(DefMechCtx%setEnergyViewer(set), IOBuffer, ierr))
         PetscCall(PetscViewerASCIIPrintf(DefMechCtx%setEnergyViewer(set), "# step     load            elastic energy  work            cohesive energy surface energy  total energy   plastic dissipation\n", ierr))
         PetscCall(PetscViewerFlush(DefMechCtx%setEnergyViewer(set), ierr))
      end do
101   format(A, '-', I4.4, '.enerblk')
102   format("# cell set ", I4, "\n")
      DefMechCtx%analysisTime = 0.0_kr
      DefMechCtx%timeStep = 0.0_kr

      !! Create Vecs and SF
      PetscCall(DMGetDimension(dm, DefMechCtx%dim, ierr))

      PetscCall(MEF90CtxGlobalOptionsSetFromOptions(DefMechCtx%MEF90Ctx%comm, trim(DefMechCtx%MEF90Ctx%prefix), MEF90CtxGlobalOptions, ierr))
      PetscCall(MEF90DefMechGlobalOptionsSetFromOptions(DefMechCtx%comm, trim(DefMechCtx%prefix), DefMechGlobalOptions, ierr))
      vecName = "Displacement"
      allocate (DefMechCtx%displacementLocal, stat=ierr)
      PetscCall(MEF90CreateLocalVector(dm, MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, DefMechCtx%dim, vecName, DefMechCtx%displacementLocal, ierr))
      allocate (DefMechCtx%displacementPreviousStepLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%displacementLocal, DefMechCtx%displacementPreviousStepLocal, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%displacementPreviousStepLocal, "DisplacementPreviousStep", ierr))
      allocate (DefMechCtx%displacementLowerBoundLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%displacementLocal, DefMechCtx%displacementLowerBoundLocal, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%displacementLowerBoundLocal, "DisplacementLowerBound", ierr))
      allocate (DefMechCtx%displacementUpperBoundLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%displacementLocal, DefMechCtx%displacementUpperBoundLocal, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%displacementUpperBoundLocal, "DisplacementUpperBound", ierr))

      allocate (DefMechCtx%damageLocal, stat=ierr)
      vecName = "Damage"
      PetscCall(MEF90CreateLocalVector(dm, MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, 1_ki, vecName, DefMechCtx%damageLocal, ierr))

      if (DefMechGlobalOptions%multiPhaseField) then
         allocate (DefMechCtx%partialDamageLocal(size(setID)), stat=ierr)
         do set = 1, size(DefMechCtx%partialDamageLocal)
            !! I need to cheat here:
            !! I create the Vec with name "Damage" so that it inherits the proper damage BC, then change its name to the proper value
            vecName = "Damage"
            PetscCall(MEF90CreateLocalVector(dm, MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, 1_ki, vecName, DefMechCtx%partialDamageLocal(set), ierr))
            write(Vecname,'("partialDamage-", I4.4)') set
            PetscCall(PetscObjectSetName(DefMechCtx%partialDamageLocal(set), vecName, ierr))
         end do
      end if
 
      allocate (DefMechCtx%TemperatureLocal, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%damageLocal, DefMechCtx%TemperatureLocal, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%TemperatureLocal, "Temperature", ierr))

      vecName = "cohesiveDisplacement"
      allocate (DefMechCtx%cohesiveDisplacement, stat=ierr)
      PetscCall(MEF90CreateLocalVector(dm, MEF90CtxGlobalOptions%elementFamily, MEF90CtxGlobalOptions%elementOrder, DefMechCtx%dim, vecName, DefMechCtx%cohesiveDisplacement, ierr))
      vecName = "bodyForce"
      allocate (DefMechCtx%bodyForce, stat=ierr)
      PetscCall(MEF90CreateCellVector(dm, DefMechCtx%dim, vecName, DefMechCtx%bodyForce, ierr))
      vecName = "boundaryForce"
      allocate (DefMechCtx%boundaryForce, stat=ierr)
      PetscCall(MEF90CreateBoundaryCellVector(dm, DefMechCtx%dim, vecName, DefMechCtx%boundaryForce, ierr))
      vecName = "pressureForce"
      allocate (DefMechCtx%pressureForce, stat=ierr)
      PetscCall(MEF90CreateBoundaryCellVector(dm, 1_ki, vecName, DefMechCtx%pressureForce, ierr))

      vecName = "plasticStrain"
      allocate (DefMechCtx%plasticStrain, stat=ierr)
      PetscCall(MEF90CreateCellVector(dm, (DefMechCtx%dim * (DefMechCtx%dim + 1_ki)) / 2_ki, vecName, DefMechCtx%plasticStrain, ierr))
      allocate (DefMechCtx%cumulatedPlasticDissipation, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%plasticStrain, DefMechCtx%cumulatedPlasticDissipation, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%cumulatedPlasticDissipation, "cumulatedPlasticDissipation", ierr))
      allocate (DefMechCtx%stress, stat=ierr)
      PetscCall(VecDuplicate(DefMechCtx%plasticStrain, DefMechCtx%stress, ierr))
      PetscCall(PetscObjectSetName(DefMechCtx%stress, "Stress", ierr))

      !! Create megaDM
      !! This needs to be modified to add the individual damage fields if needed
      allocate (dmList(7))
      PetscCall(VecGetDM(DefMechCtx%displacementLocal, dmList(1), ierr))
      PetscCall(VecGetDM(DefMechCtx%damageLocal, dmList(2), ierr))
      PetscCall(VecGetDM(DefMechCtx%temperatureLocal, dmList(3), ierr))
      PetscCall(VecGetDM(DefMechCtx%bodyForce, dmList(4), ierr))
      PetscCall(VecGetDM(DefMechCtx%boundaryForce, dmList(5), ierr))
      PetscCall(VecGetDM(DefMechCtx%pressureForce, dmList(6), ierr))
      PetscCall(VecGetDM(DefMechCtx%plasticStrain, dmList(7), ierr))
      PetscCall(DMCreateSuperDM(dmList, 7_ki, PETSC_NULL_IS_POINTER, DefMechCtx%megaDM, ierr))
      deallocate (dmList)

      !! Create the IO SF for all fields
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%displacementLocal, DefMechCtx%displacementToIOSF, DefMechCtx%IOTodisplacementSF, ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%damageLocal, DefMechCtx%damageToIOSF, DefMechCtx%IOTodamageSF, ierr))

      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%cohesiveDisplacement, DefMechCtx%cohesiveDisplacementToIOSF, DefMechCtx%IOToCohesiveDisplacementSF, ierr))

      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%temperatureLocal, DefMechCtx%temperatureToIOSF, DefMechCtx%IOTotemperatureSF, ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%bodyForce, DefMechCtx%bodyForceToIOSF, DefMechCtx%IOTobodyForceSF, ierr))
      ! PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMech%boundaryForce,DefMech%boundaryForceToIOSF,DefMech%IOToboundaryForceSF,ierr))
      ! PetscCall(MEF90IOSFCreate(MEF90Ctx,DefMech%pressureForce,DefMech%pressureForceToIOSF,DefMech%IOTopressureForceSF,ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%stress, DefMechCtx%stressToIOSF, DefMechCtx%IOToStressSF, ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%plasticStrain, DefMechCtx%plasticStrainToIOSF, DefMechCtx%IOToplasticStrainSF, ierr))
      PetscCall(MEF90IOSFCreate(MEF90Ctx, DefMechCtx%cumulatedPlasticDissipation, DefMechCtx%cumulatedPlasticDissToIOSF, DefMechCtx%IOToCumulatedPlasticDissSF, ierr))

      !! Create the SF to exchange boundary values of the displacement and damage.
      PetscCall(MEF90ConstraintSFCreate(DefMechCtx%MEF90Ctx, DefMechCtx%displacementLocal, DefMechCtx%displacementLocal, DefMechCtx%displacementConstraintsSF, dummySF, ierr))
      PetscCall(PetscSFDestroy(dummySF, ierr))
      PetscCall(MEF90ConstraintSFCreate(DefMechCtx%MEF90Ctx, DefMechCtx%damageLocal, DefMechCtx%damageLocal, DefMechCtx%damageConstraintsSF, dummySF, ierr))
      PetscCall(PetscSFDestroy(dummySF, ierr))

      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

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

      PetscInt                                         :: set

      DefMech%PETScCtx = C_NULL_PTR

      !!
      !! Close energy viewers
      !!
      PetscCall(PetscViewerDestroy(DefMech%globalEnergyViewer, ierr))
      do set = 1, size(DefMech%setEnergyViewer)
         PetscCall(PetscViewerDestroy(DefMech%setEnergyViewer(set), ierr))
      end do
      deallocate (DefMech%setEnergyViewer)

      !! Destroy Vecs and SF and deAllocate them
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

      if (associated(DefMech%partialDamageLocal)) then
         do set = 1, size(DefMech%partialDamageLocal)
            PetscCall(VecDestroy(DefMech%partialDamageLocal(set), ierr))
         end do
         deallocate (DefMech%partialDamageLocal)
         nullify (DefMech%partialDamageLocal)
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

      !! Destroy all PetscSF
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

      !! Destroy the SF to exchange boundary values of the displacement and damage.
      PetscCall(PetscSFDestroy(DefMech%displacementConstraintsSF, ierr))
      PetscCall(PetscSFDestroy(DefMech%damageConstraintsSF, ierr))

      !! Destroy the megaDM
      PetscCall(DMDestroy(DefMech%megaDM, ierr))
   end subroutine MEF90DefMechDestroy

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechGlobalOptionsSetFromOptions"
!!!
!!!
!!!  MEF90DefMechGlobalOptionsSetFromOptions: reads the problem-wide options of a MEF90DefMech_Type
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechGlobalOptionsSetFromOptions(comm, prefix, options, ierr)
      MPIU_Comm, intent(IN)                                  :: comm
      character(len=*), intent(IN)                           :: prefix
      type(MEF90DefMechGlobalOptions_Type), intent(OUT)      :: options
      PetscErrorCode, intent(INOUT)                          :: ierr

      PetscInt                                               :: verbose = 0

      options%comm = comm
      options%prefix = prefix
      options%name = trim(prefix)//"DefMechGlobalOptions"

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
         PetscCall(PetscOptionsBool('-multiPhaseField', 'Use one damage variable per cell set', 'mef90DefMech', options%multiPhaseField, options%multiPhaseField, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsEnum('-BT_Type', 'Backtracking type', 'mef90DefMech', MEF90DefMech_BTTypeList, options%BTType, options%BTType, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-BT_Interval', 'Interval at which Backtracking is run in inner loop (0 for outer loop)', 'mef90DefMech', options%BTInterval, options%BTInterval, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsInt('-BT_Scope', 'Backtracking scope (0 for unlimited)', 'mef90DefMech', options%BTScope, options%BTScope, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-BT_Tol', 'Backtracking relative tolerance', 'mef90DefMech', options%BTTol, options%BTTol, PETSC_NULL_BOOL, ierr))

         PetscCall(PetscOptionsReal('-defmech_plasticstrain_atol', 'Absolute tolerance on plastic error', 'mef90DefMech', options%plasticStrainATol, options%plasticStrainATol, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-defmech_InjectedVolume_atol', 'Absolute tolerance on injected volume error', 'mef90DefMech', options%InjectedVolumeATol, options%InjectedVolumeATol, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-defmech_dampingCoefficient_displacement', 'Damping coefficient on displacement field (0 for minimization, 1 for semi-implicit gradient flow)', 'mef90DefMech', options%dampingCoefficientDisplacement, options%dampingCoefficientDisplacement, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-defmech_dampingCoefficient_damage', 'Damping coefficient on damage field (0 for minimization, 1 for semi-implicit gradient flow)', 'mef90DefMech', options%dampingCoefficientDamage, options%dampingCoefficientDamage, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))

      PetscCall(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", verbose, PETSC_NULL_BOOL, ierr))
      if (verbose > 0) then
         !! PETSC_VIEWER_STDOUT_WORLD and not PetscViewerASCIIGetStdout: the latter is collective, and the
         !! options are read where they are needed, which is not always on all ranks. view itself only
         !! prints the first time a given set is read, so it is collective only then.
         call options%view(PETSC_VIEWER_STDOUT_WORLD, ierr)
      end if
   end subroutine MEF90DefMechGlobalOptionsSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechGlobalOptionsView"
!!!
!!!
!!!  MEF90DefMechGlobalOptionsView: the default viewer for a MEF90DefMechGlobalOptions_Type
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechGlobalOptionsView(self, viewer, ierr)
      class(MEF90DefMechGlobalOptions_Type), intent(IN)      :: self
      type(tPetscViewer), intent(IN)                         :: viewer
      PetscErrorCode, intent(INOUT)                          :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)              :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)              :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType /= 'ascii') return

      write (IOBuffer, "(A,': Options for MEF90DefMechGlobalOptions_Type\n')") trim(self%prefix)//"defmech"
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         time stepping type: ',A,'\n')") trim(MEF90DefMech_TimeSteppingTypeList(self%timeSteppingType + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         solver type: ',A,'\n')") trim(MEF90DefMech_SolverTypeList(self%solverType + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         damage solver type: ',A,'\n')") trim(MEF90DefMech_DamageSolverTypeList(self%damageSolverType + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         multi phase field: ',L1,'\n')") self%multiPhaseField
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         damage atol / maxit: ',ES12.5,' / ',I6,'\n')") self%damageATol, self%damageMaxIt
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         PC lag: ',I6,' SOR omega: ',ES12.5,' irreversibility threshold: ',ES12.5,'\n')") &
         self%PCLag, self%SOROmega, self%irrevthres
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         backtracking type: ',A,' interval: ',I6,' scope: ',I6,' tol: ',ES12.5,'\n')") &
         trim(MEF90DefMech_BTTypeList(self%BTType + 1)), self%BTInterval, self%BTScope, self%BTTol
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         plastic strain atol / injected volume atol: ',ES12.5,' / ',ES12.5,'\n')") &
         self%plasticStrainATol, self%InjectedVolumeATol
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         damping coefficient displacement / damage: ',ES12.5,' / ',ES12.5,'\n')") &
         self%dampingCoefficientDisplacement, self%dampingCoefficientDamage
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         scaling boundary displacement / lower bound / upper bound: ',A,' / ',A,' / ',A,'\n')") &
         trim(MEF90ScalingList(self%boundaryDisplacementScaling + 1)), &
         trim(MEF90ScalingList(self%displacementLowerBoundScaling + 1)), &
         trim(MEF90ScalingList(self%displacementUpperBoundScaling + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         scaling cohesive displacement / boundary damage: ',A,' / ',A,'\n')") &
         trim(MEF90ScalingList(self%cohesiveDisplacementScaling + 1)), trim(MEF90ScalingList(self%boundaryDamageScaling + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         scaling body force / boundary force / pressure force / crack pressure: ',A,' / ',A,' / ',A,' / ',A,'\n')") &
         trim(MEF90ScalingList(self%bodyForceScaling + 1)), trim(MEF90ScalingList(self%boundaryForceScaling + 1)), &
         trim(MEF90ScalingList(self%pressureForceScaling + 1)), trim(MEF90ScalingList(self%CrackPressureScaling + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         export displacement / damage / stress: ',3(L1,' '),'\n')") &
         self%displacementExport, self%damageExport, self%stressExport
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         export temperature / plastic strain / cumulated plastic dissipation: ',3(L1,' '),'\n')") &
         self%temperatureExport, self%plasticStrainExport, self%cumulatedPlasticDissipationExport
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
   end subroutine MEF90DefMechGlobalOptionsView

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCellSetOptionsSetFromOptions"
!!!
!!!
!!!  MEF90DefMechCellSetOptionsSetFromOptions: reads the options of a single cell set
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechCellSetOptionsSetFromOptions(comm, prefix, dim, options, ierr)
      MPIU_Comm, intent(IN)                                  :: comm
      character(len=*), intent(IN)                           :: prefix
      PetscInt, intent(IN)                                   :: dim
      type(MEF90DefMechCellSetOptions_Type), intent(OUT)     :: options
      PetscErrorCode, intent(INOUT)                          :: ierr

      PetscInt                                               :: nOpt
      PetscReal, dimension(:), allocatable                   :: tmpArray
      PetscInt                                               :: verbose = 0

      options%comm = comm
      options%prefix = prefix
      options%name = trim(prefix)//"DefMechCellSetOptions"

      select case (dim)
      case (2)
         options%linearThermalExpansion = MEF90MatS2DIdentity
      case (3)
         options%linearThermalExpansion = MEF90MatS3DIdentity
      case default
         SETERRQ(comm, PETSC_ERR_ARG_OUTOFRANGE, "dim must be 2 or 3 in " // __FUNCT__)
      end select

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
         PetscCall(PetscOptionsReal('-cohesiveStiffness', '[N.m^(-4)] (k) cohesive stiffness in Winkler-type models', 'mef90DefMech', options%cohesiveStiffness, options%cohesiveStiffness, PETSC_NULL_BOOL, ierr))
         select type (alpha => options%linearThermalExpansion)
         type is (MatS2D)
            nOpt = 3
            allocate (tmpArray(nOpt))
            tmpArray = alpha
            PetscCall(PetscOptionsRealArray('-LinearThermalExpansion', '[K^(-1)] (alpha) Linear thermal expansion matrix', 'mef90DefMech', tmpArray, nOpt, PETSC_NULL_BOOL, ierr))
            alpha = tmpArray
            deallocate (tmpArray)
         type is (MatS3D)
            nOpt = 6
            allocate (tmpArray(nOpt))
            tmpArray = alpha
            PetscCall(PetscOptionsRealArray('-LinearThermalExpansion', '[K^(-1)] (alpha) Linear thermal expansion matrix', 'mef90DefMech', tmpArray, nOpt, PETSC_NULL_BOOL, ierr))
            alpha = tmpArray
            deallocate (tmpArray)
         end select
      PetscCall(PetscOptionsEnd(ierr))

      PetscCall(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", verbose, PETSC_NULL_BOOL, ierr))
      if (verbose > 0) then
         !! PETSC_VIEWER_STDOUT_WORLD and not PetscViewerASCIIGetStdout: the latter is collective, and the
         !! options are read where they are needed, which is not always on all ranks. view itself only
         !! prints the first time a given set is read, so it is collective only then.
         call options%view(PETSC_VIEWER_STDOUT_WORLD, ierr)
      end if
   end subroutine MEF90DefMechCellSetOptionsSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechCellSetOptionsView"
!!!
!!!
!!!  MEF90DefMechCellSetOptionsView: the default viewer for a MEF90DefMechCellSetOptions_Type
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechCellSetOptionsView(self, viewer, ierr)
      class(MEF90DefMechCellSetOptions_Type), intent(IN)     :: self
      type(tPetscViewer), intent(IN)                         :: viewer
      PetscErrorCode, intent(INOUT)                          :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)              :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)              :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType /= 'ascii') return

      write (IOBuffer, "(A,': Options for MEF90DefMechCellSetOptions_Type\n')") trim(self%prefix)
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         body force: ',3(ES12.5,' '),'crack pressure: ',ES12.5,'\n')") &
         self%bodyForce, self%crackPressure
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         plasticity type: ',A,'\n')") trim(MEF90DefMech_plasticityTypeList(self%plasticityType + 1))
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         cohesive displacement: ',3(ES12.5,' '),'stiffness: ',ES12.5,'\n')") &
         self%cohesiveDisplacement, self%cohesiveStiffness
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         displacement BC: ',3(L1,' '),'value: ',3(ES12.5,' '),'\n')") &
         self%Has_displacementBC, self%boundaryDisplacement
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         displacement lower / upper bound: ',3(ES12.5,' '),' / ',3(ES12.5,' '),'\n')") &
         self%displacementLowerBound, self%displacementUpperBound
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         damage BC: ',L1,' value: ',ES12.5,'\n')") &
         self%Has_damageBC, self%boundaryDamage
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         crack volume controlled / work controlled: ',L1,' / ',L1,'\n')") &
         self%CrackVolumeControlled, self%WorkControlled
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      if (allocated(self%linearThermalExpansion)) then
         select type (alpha => self%linearThermalExpansion)
         type is (MatS2D)
            write (IOBuffer, "('         linear thermal expansion (alpha): ',2(ES12.5,', '),ES12.5,' [K^(-1)]\n')") alpha
         type is (MatS3D)
            write (IOBuffer, "('         linear thermal expansion (alpha): ',5(ES12.5,', '),ES12.5,' [K^(-1)]\n')") alpha
         class default
            write (IOBuffer, *) 'somehow linear thermal expansion is neither MatS2D nor MatS3D. This is wrong\n'
         end select
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      end if
   end subroutine MEF90DefMechCellSetOptionsView

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
      PetscInt                                               :: verbose = 0

      options%comm = comm
      options%prefix = prefix
      options%name = trim(prefix)//"DefMechFaceSetOptions"

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

      PetscCall(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", verbose, PETSC_NULL_BOOL, ierr))
      if (verbose > 0) then
         !! PETSC_VIEWER_STDOUT_WORLD and not PetscViewerASCIIGetStdout: the latter is collective, and the
         !! options are read where they are needed, which is not always on all ranks. view itself only
         !! prints the first time a given set is read, so it is collective only then.
         call options%view(PETSC_VIEWER_STDOUT_WORLD, ierr)
      end if
   end subroutine MEF90DefMechFaceSetOptionsSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechFaceSetOptionsView"
!!!
!!!
!!!  MEF90DefMechFaceSetOptionsView: the default viewer for a MEF90DefMechFaceSetOptions_Type
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechFaceSetOptionsView(self, viewer, ierr)
      class(MEF90DefMechFaceSetOptions_Type), intent(IN)     :: self
      type(tPetscViewer), intent(IN)                         :: viewer
      PetscErrorCode, intent(INOUT)                          :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)              :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)              :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType /= 'ascii') return

      write (IOBuffer, "(A,': Options for MEF90DefMechFaceSetOptions_Type\n')") trim(self%prefix)
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         boundary force: ',3(ES12.5,' '),'pressure force: ',ES12.5,'\n')") &
         self%boundaryForce, self%pressureForce
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         displacement BC: ',3(L1,' '),'value: ',3(ES12.5,' '),'\n')") &
         self%Has_displacementBC, self%boundaryDisplacement
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         displacement lower / upper bound: ',3(ES12.5,' '),' / ',3(ES12.5,' '),'\n')") &
         self%displacementLowerBound, self%displacementUpperBound
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         damage BC: ',L1,' value: ',ES12.5,'\n')") &
         self%Has_damageBC, self%boundaryDamage
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
   end subroutine MEF90DefMechFaceSetOptionsView

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
      PetscInt                                               :: verbose = 0

      options%comm = comm
      options%prefix = prefix
      options%name = trim(prefix)//"DefMechVertexSetOptions"

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

      PetscCall(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", verbose, PETSC_NULL_BOOL, ierr))
      if (verbose > 0) then
         !! PETSC_VIEWER_STDOUT_WORLD and not PetscViewerASCIIGetStdout: the latter is collective, and the
         !! options are read where they are needed, which is not always on all ranks. view itself only
         !! prints the first time a given set is read, so it is collective only then.
         call options%view(PETSC_VIEWER_STDOUT_WORLD, ierr)
      end if
   end subroutine MEF90DefMechVertexSetOptionsSetFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechVertexSetOptionsView"
!!!
!!!
!!!  MEF90DefMechVertexSetOptionsView: the default viewer for a MEF90DefMechVertexSetOptions_Type
!!!
!!!  (c) 2026 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90DefMechVertexSetOptionsView(self, viewer, ierr)
      class(MEF90DefMechVertexSetOptions_Type), intent(IN)   :: self
      type(tPetscViewer), intent(IN)                         :: viewer
      PetscErrorCode, intent(INOUT)                          :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)              :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)              :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType /= 'ascii') return

      write (IOBuffer, "(A,': Options for MEF90DefMechVertexSetOptions_Type\n')") trim(self%prefix)
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         displacement BC: ',3(L1,' '),'value: ',3(ES12.5,' '),'\n')") &
         self%Has_displacementBC, self%boundaryDisplacement
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         displacement lower / upper bound: ',3(ES12.5,' '),' / ',3(ES12.5,' '),'\n')") &
         self%displacementLowerBound, self%displacementUpperBound
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
      write (IOBuffer, "('         damage BC: ',L1,' value: ',ES12.5,'\n')") &
         self%Has_damageBC, self%boundaryDamage
      PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))
   end subroutine MEF90DefMechVertexSetOptionsView

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
      type(MEF90DefMechGlobalOptions_Type)                   :: globalOptions
      type(MEF90DefMechCellSetOptions_Type)                  :: cellSetOptions
      type(MEF90DefMechFaceSetOptions_Type)                  :: faceSetOptions
      PetscBool                                              :: printHelp

      !!
      !! Problem-wide options
      !!
      PetscCall(MEF90DefMechGlobalOptionsSetFromOptions(self%comm, trim(self%prefix), globalOptions, ierr))

      !!
      !! The per-set options are read where they are needed, but they have to be visited once here so
      !! that they are registered with the options database, and so that the displacement bounds of all
      !! sets are known before the solvers are set up
      !!
      self%hasDisplacementBounds = PETSC_FALSE

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90CellSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'cs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90DefMechCellSetOptionsSetFromOptions(self%comm, trim(setPrefix), self%dim, cellSetOptions, ierr))
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

      character(len=MEF90MXSTRLEN, kind=c_char)              :: viewerType
      type(tIS)                                              :: setIS
      PetscInt, dimension(:), pointer                        :: setID
      PetscInt                                               :: set
      character(len=MEF90MXSTRLEN)                           :: setPrefix
      type(MEF90DefMechGlobalOptions_Type)                   :: globalOptions
      type(MEF90DefMechCellSetOptions_Type)                  :: cellSetOptions
      type(MEF90DefMechFaceSetOptions_Type)                  :: faceSetOptions
      type(MEF90DefMechVertexSetOptions_Type)                :: vertexSetOptions

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType /= 'ascii') return

      !! The problem-wide and per-set options print themselves: each of the option types below has its
      !! own view, and MEF90Object%view makes sure that a given set is only ever printed once.
      PetscCall(MEF90DefMechGlobalOptionsSetFromOptions(self%comm, trim(self%prefix), globalOptions, ierr))
      call globalOptions%view(viewer, ierr)

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90CellSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'cs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90DefMechCellSetOptionsSetFromOptions(self%comm, trim(setPrefix), self%dim, cellSetOptions, ierr))
         call cellSetOptions%view(viewer, ierr)
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90FaceSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'fs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90DefMechFaceSetOptionsSetFromOptions(self%comm, trim(setPrefix), faceSetOptions, ierr))
         call faceSetOptions%view(viewer, ierr)
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))

      PetscCall(DMGetLabelIdIS(self%megaDM, MEF90VertexSetLabelName, setIS, ierr))
      PetscCall(MEF90ISAllGatherMerge(self%comm, setIS, ierr))
      PetscCall(ISGetIndices(setIS, setID, ierr))
      do set = 1, size(setID)
         write (setPrefix, "(A,'vs',I4.4,'_')") trim(self%prefix), setID(set)
         PetscCall(MEF90DefMechVertexSetOptionsSetFromOptions(self%comm, trim(setPrefix), vertexSetOptions, ierr))
         call vertexSetOptions%view(viewer, ierr)
      end do
      PetscCall(ISRestoreIndices(setIS, setID, ierr))
      PetscCall(ISDestroy(setIS, ierr))
   end subroutine MEF90DefMechView
end module m_MEF90_DefMech_class
