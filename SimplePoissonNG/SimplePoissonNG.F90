!!!
!!! Try: mpiexec -n 2 ./Darwin-intel11.1-mef90-g/SimplePoissonNG2D -vs0001_tempBC 1 -vs0001_temp 0 -vs0002_tempBC 1 -vs0002_temp 0 -cs0001_force 0. -cs0002_force 1. -temp_pc_type bjacobi -temp_ksp_monitor -temp_snes_type ksponly --prefix SquareNG2
!!! 
Program  SimplePoissonNG
#include "SimplePoisson.inc"
#include <finclude/petscdef.h>
#include <finclude/petscbagdef.h>
   Use m_MEF90
   Use m_PoissonGlobalProperties
   Use M_POISSONCELLSETPROPERTIES
   Use m_PoissonVertexSetProperties
   Use M_POISSONASSEMBLY
   Use m_Poisson

   Use petsc
   Implicit NONE   

   Type(PoissonGlobalProperties_Type),parameter    :: defaultGlobalProperties = PoissonGlobalProperties_Type( &
                                                         0,                   & ! verbose
                                                         Poisson_EXOSingle,   & ! FileFormat
                                                         Poisson_Replace,     & ! FileMode
                                                         Poisson_SteadyState, & ! EvolutionLaw
                                                         Poisson_MIL,         & ! LoadingType
                                                         1.0_Kr,              & ! TimeMin
                                                         1.0_kr,              & ! TimeMax
                                                         1,                   & ! numTimeStep
                                                         1,                   & ! tempOffset
                                                         1,                   & ! refOffset
                                                         1,                   & ! bcOffset
                                                         2,                   & ! fluxoffset
                                                         PETSC_FALSE,         & ! addNullSpace
                                                         0.0_Kr)                ! initialTemp
                                                         
   Type(PoissonCellSetProperties_Type)             :: defaultCellSetProperties   = PoissonCellSetProperties_Type(DEFAULT_ELEMENT_SHORTID,0.0_Kr,0.0_Kr,0.0_Kr)
   Type(PoissonVertexSetProperties_Type),parameter :: defaultVertexSetProperties = PoissonVertexSetProperties_Type(PETSC_TRUE,0)
   Type(PoissonGlobalProperties_Type),pointer      :: GlobalProperties 
   Character(len=MEF90_MXSTRLEN)                   :: prefix,filename
   Type(PetscViewer),target                        :: energyViewer,logViewer
   Type(PetscViewer),Dimension(:),Pointer          :: energyViewerCellSet
   Type(MEF90Ctx_Type)                             :: MEF90Ctx
   Type(DM),target                                        :: mesh,tmp_mesh
   PetscErrorCode                                  :: iErr
   Character(len=MEF90_MXSTRLEN)                   :: IOBuffer
   PetscBool                                       :: flg
   Type(SNES),target                               :: snesTemp
   Type(TS),target                                 :: tsTemp
   Type(KSP),target                                :: kspTemp
   Type(PC),target                                 :: pcTemp
   Type(Mat),target                                :: matTemp
   Type(Vec),target                                :: solTemp,resTemp,RHSTemp
   Type(Vec),target                                :: flux,refTemp,bcTemp
   Integer                                         :: exoIN=0,exoOUT=0
   PetscReal,Dimension(:),Pointer                  :: energy,work
   PetscInt,dimension(:),Pointer                   :: CellSetGlobalID
   PetscInt                                        :: set
   Type(MatNullSpace),target                       :: nspTemp
   SNESConvergedReason                             :: reasonTemp
   PetscInt                                        :: itsTemp
   Type(SectionReal),target                        :: secTemp
   Type(IS),target                                 :: CellSetGlobalIS
   Type(Element_Type),Dimension(:),pointer         :: ElemType
      
   PetscInt                                        :: point,numCell,numVertex,numCellSetGlobal
   PetscReal,Dimension(:,:),Pointer                :: Coord
   PetscReal,Dimension(:),Pointer                  :: Coord2
   PetscInt                                        :: i,j
   PetscInt                                        :: TimeStepNum=1
   PetscReal,Dimension(:),Pointer                  :: Time
   
   PetscReal                                       :: rtol,atol,dtol
   PetscReal                                       :: TSinitialTime,TSsolutionTime
   PetscReal                                       :: TSinitialStep
   PetscInt                                        :: TSmaxStep=5000
   PetscInt                                        :: TSmaxTime
   PetscInt                                        :: maxits
   
   Call MEF90_Initialize()
   Call m_Poisson_Initialize(ierr);CHKERRQ(ierr)

   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-prefix',prefix,flg,ierr);CHKERRQ(ierr)
   If (.NOT. flg) Then
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_NULL,"Missing file prefix\n",ierr);CHKERRQ(ierr)
   End If

   Call MEF90CtxPoissonGlobalPropertiesCreate(MEF90Ctx,defaultGlobalProperties,ierr)
   Call PetscBagGetDataPoissonGlobalProperties(MEF90Ctx%GlobalPropertiesBag,GlobalProperties,ierr);CHKERRQ(ierr)  
   Call SimplePoissonEXOOpenInputFile(prefix,exoIN,ierr);CHKERRQ(ierr)

!<<<<<<
   !!!
   !!! Read DMMesh from exoIN
   !!!
   If (MEF90_NumProcs == 1) Then
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoIN,mesh,ierr);CHKERRQ(ierr)
   Else
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoIN,Tmp_mesh,ierr);CHKERRQ(ierr)
   
      Call DMMeshDistribute(Tmp_mesh,PETSC_NULL_CHARACTER,mesh,ierr);CHKERRQ(ierr)
      Call DMDestroy(Tmp_mesh,ierr);CHKERRQ(ierr)
   End If
   !!! 
   !!! Create a Section consistent with the element choice
   !!! The section is named 'default' so that it can be picked as the default
   !!! layout for all DM vector creation routines
   !!!
   Call DMMeshGetStratumSize(mesh,"depth",0,numVertex,ierr);CHKERRQ(ierr)
   Call DMMeshGetStratumSize(mesh,"height",0,numCell,ierr);CHKERRQ(ierr)
   Call DMMeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)   
   Call ISGetLocalSize(CellSetGlobalIS,numCellSetGlobal,ierr);CHKERRQ(ierr)
   Call ISGetIndicesF90(CellSetGlobalIS,CellSetGlobalID,ierr);CHKERRQ(ierr)


!!! This may go into a separate function since I would need to know the 
!!! element type in order to allocate anything more complicated than a 
!!! vertex or cell based section
   Call DMMeshGetSectionReal(mesh,'default',SecTemp,ierr);CHKERRQ(ierr)
   Do point = numCell,numCell+numVertex-1
      Call SectionRealSetFiberDimension(SecTemp,point,1,ierr);CHKERRQ(ierr)
   End Do
   Call SectionRealAllocate(SecTemp,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(secTemp,ierr);CHKERRQ(ierr)
   !!! I can destroy the section because its layout remains somehow cached in the DMMesh
!
! THIS NEEDS TO GO BACK INTO MEF90
!>>>   


   Call SimplePoissonGetTime(time,exoIN,MEF90Ctx,ierr);CHKERRQ(ierr)
   !!! This also needs to go into MEF90
   !!! I need to think about ways to make the interface generic
   
   
   !!! 
   !!! Create SNES or TS
   !!!
   If (GlobalProperties%TimeEvolution == Poisson_SteadyState) Then
      Call SNESCreate(PETSC_COMM_WORLD,snesTemp,ierr);CHKERRQ(ierr)
      Call SNESSetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
      Call SNESSetOptionsPrefix(snesTemp,'temp_',ierr);CHKERRQ(ierr)
   Else
      Call TSCreate(PETSC_COMM_WORLD,tsTemp,ierr);CHKERRQ(ierr)
      Call TSSetDM(tsTemp,mesh,ierr);CHKERRQ(ierr)
      Call TSSetOptionsPrefix(tsTemp,'temp_',ierr);CHKERRQ(ierr)
      Call TSGetSNES(tsTemp,snesTemp,ierr);CHKERRQ(ierr)
   End If

   !!! 
   !!! Create PoissonCtx
   !!!
   Call MEF90CtxPoissonVertexSetPropertiesCreate(MEF90Ctx,snesTemp,defaultVertexSetProperties,ierr);CHKERRQ(ierr)
   Call EXOGetCellSetElementType_Scal(PETSC_COMM_WORLD,exoIN,MEF90_DIM,ElemType,ierr)
   Call MEF90CtxPoissonCellSetPropertiesCreate(MEF90Ctx,snesTemp,defaultCellSetProperties,ElemType,ierr);CHKERRQ(ierr)


   !!! Stop here if -help is passed
   !!!
   !Call PetscOptionsGetBool(PETSC_NULL_CHARACTER,"-help",flg,flg,ierr);CHKERRQ(ierr)
   !If (flg) Then
   !   Call MEF90_Finalize()
   !   STOP
   !End If
   

   !!! Prepare output file
   Call SimplePoissonPrepareOutputEXO(mesh,prefix,exoIN,exoOUT,MEF90Ctx,ierr);CHKERRQ(ierr)
   
   !!! Open energy and log files
   Write(IOBuffer,*) 'step  time          energy        work          total\n'
   filename = trim(prefix) // '.log'
   Call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,logViewer,ierr);CHKERRQ(ierr)
   filename = trim(prefix) // '.ener'
   Call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,energyViewer,ierr);CHKERRQ(ierr)
      Call PetscViewerASCIIPrintf(energyViewer,IOBuffer,ierr);CHKERRQ(ierr)
   Allocate(energyViewerCellSet(numCellSetGlobal))
   Do set = 1, numCellSetGlobal
      Write(filename,300) trim(prefix),CellSetGlobalID(set)
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,energyViewerCellSet(set),ierr);CHKERRQ(ierr)
      Call PetscViewerASCIIPrintf(energyViewerCellSet(set),IOBuffer,ierr);CHKERRQ(ierr)
   End Do
300 Format(A,'-blk_',I4.4,'.ener')

   !!!
   !!! Get Matrix for the Jacobian / SNES and unknown vector
   !!!
   Call DMMeshSetMaxDof(mesh,1,iErr); CHKERRQ(iErr) 
   !Call DMMeshCreateMatrix(mesh,secTemp,MATAIJ,matTemp,iErr);CHKERRQ(iErr)
   Call DMCreateMatrix(mesh,MATAIJ,matTemp,iErr);CHKERRQ(iErr)
   ! I can use the generic DMCreateMatrix because I decided to associate a single DM with each 
   ! dof layout (just like DMComplex or DMDA do)

   !!! There are two ways to get global and local vectors:
   !!! Through the generic DM interface, which will only pull the layout from the 'default' section
   Call DMGetGlobalVector(mesh,solTemp,ierr);CHKERRQ(ierr)
   !!! Through DMMeshCreateVector which lets use any layout
   !Call DMMeshCreateVector(mesh,secTemp,solTemp,ierr);CHKERRQ(ierr)

   Call VecDuplicate(solTemp,resTemp,ierr);CHKERRQ(ierr)
   Call VecDuplicate(solTemp,RHSTemp,ierr);CHKERRQ(ierr)
   
   If (GlobalProperties%LoadingType == Poisson_FILE) Then
      Call VecDuplicate(solTemp,flux,ierr);CHKERRQ(ierr)
      Call VecDuplicate(solTemp,refTemp,ierr);CHKERRQ(ierr)
      Call VecDuplicate(solTemp,bcTemp,ierr);CHKERRQ(ierr)
   End If
      
   !!! Adding a null space when some boundary conditions are prescribes breaks everything...
   !!! Need to add a flag and make adding the null space optional
   If (GlobalProperties%addNullSpace) Then
      Call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL_OBJECT,nspTemp,ierr);CHKERRQ(ierr)
      Call MatSetNullSpace(matTemp,nspTemp,ierr);CHKERRQ(ierr)
   End If

   !!! Not sure if this is still needed when using MatZeroRowsColumnsIS for BC handling
   Call MatSetOption(matTemp,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetFromOptions(matTemp,ierr);CHKERRQ(ierr)

   !!!
   !!! Set Jacobian and Function for the SNES / TS
   !!!
   If (GlobalProperties%TimeEvolution == Poisson_SteadyState) Then
      Call SNESSetFunction(snesTemp,resTemp,SimplePoissonOperator,MEF90Ctx,ierr);CHKERRQ(ierr)
      Call SNESSetJacobian(snesTemp,matTemp,matTemp,SimplePoissonBilinearForm,MEF90Ctx,ierr);CHKERRQ(ierr)
      Call SNESSetFromOptions(snesTemp,ierr);CHKERRQ(ierr)
   Else
      Call TSSetIFunction(tsTemp,resTemp,SimplePoissonTSIFunction,MEF90Ctx,ierr);CHKERRQ(ierr)
      Call TSSetIJacobian(tsTemp,matTemp,matTemp,SimplePoissonTSIJacobian,MEF90Ctx,ierr);CHKERRQ(ierr)
      If ((GlobalProperties%LoadingType == Poisson_MIL) .OR. (GlobalProperties%LoadingType == Poisson_CST)) Then
         Call TSSetRHSFunction(tsTemp,PETSC_NULL_OBJECT,SimplePoissonTSRHS_Cst,MEF90Ctx,ierr);CHKERRQ(ierr)
      Else
         !Call TSSetRHSFunction(tsTemp,PETSC_NULL_OBJECT,SimplePoissonTSRHS,MEF90Ctx,ierr);CHKERRQ(ierr)
      End If
      Call TSSetType(tsTemp,'rosw',ierr);CHKERRQ(ierr)
      Call TSRosWSetType(tsTemp,'ra3pw',ierr);CHKERRQ(ierr)
      TSinitialStep = (time(size(time))-time(0)) / (size(time) + 0.0_Kr) / 10.0_Kr
      TSinitialTime = time(1)
      Call TSSetInitialTimeStep(tsTemp,TSinitialTime,TSinitialStep,ierr);CHKERRQ(ierr)
      Call TSSetProblemType(tsTemp,TS_LINEAR,ierr);CHKERRQ(ierr)
      Call VecSet(solTemp,GlobalProperties%initialTemp,ierr);CHKERRQ(ierr)
      Call TSSetSolution(tsTemp,solTemp,ierr);CHKERRQ(ierr)
      Call TSSetFromOptions(tsTemp,ierr);CHKERRQ(ierr)
   End If
   
   !!! 
   !!! Set some KSP options
   !!!
   Call SNESGetKSP(snesTemp,kspTemp,ierr);CHKERRQ(ierr)
   Call KSPSetType(kspTemp,KSPCG,ierr);CHKERRQ(ierr)
   Call KSPSetInitialGuessNonzero(kspTemp,PETSC_TRUE,ierr);CHKERRQ(ierr)
   If (GlobalProperties%addNullSpace) Then
      Call KSPSetNullSpace(kspTemp,nspTemp,ierr);CHKERRQ(ierr)
   End If
   Call KSPGetTolerances(kspTemp,rtol,atol,dtol,maxits,ierr);CHKERRQ(ierr)
   rtol = 1.0D-8
   Call KSPSetTolerances(kspTemp,rtol,atol,dtol,maxits,ierr);CHKERRQ(ierr)
   Call KSPSetFromOptions(kspTemp,ierr);CHKERRQ(ierr)

   !!!
   !!! Set some PC options
   !!!
   Call KSPGetPC(kspTemp,pcTemp,ierr);CHKERRQ(ierr)
   Call PCSetFromOptions(pcTemp,ierr);CHKERRQ(ierr)
   Call DMMeshGetCoordinatesF90(mesh,Coord,ierr);CHKERRQ(ierr)
   Allocate(coord2(size(coord)))
   Do i = 1, size(coord,1)
      Do j = 1, MEF90_DIM
         coord2(MEF90_DIM * (i-1) + j) = coord(i,j)
      End Do
   End Do
   Call PCSetCoordinates(pcTemp,MEF90_DIM,size(coord2)/MEF90_DIM,coord,ierr)
   Call DMMeshRestoreCoordinatesF90(mesh,Coord,ierr);CHKERRQ(ierr)
   DeAllocate(coord2)

   !!!
   !!! Allocate arrays for the energies in each block
   !!!
   Allocate(energy(numCellSetGlobal),stat=ierr)
   Allocate(work(numCellSetGlobal),stat=ierr)
   
   !!!
   !!! MAIN LOOP: Solve Poisson Equation
   !!!
   Do TimeStepNum = 1, size(time)
      Write(IOBuffer,200) TimeStepNum,time(TimeStepNum)
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
200   Format('Solving time step ',I4,': t=',ES12.5,'\n')

      If ( GlobalProperties%TimeEvolution == Poisson_SteadyState) Then
         Select case(GlobalProperties%LoadingType)
         Case(Poisson_MIL)
         
            Call SimplePoissonFormInitialGuess_Cst(snesTemp,solTemp,time(timeStepNum),MEF90Ctx,ierr);CHKERRQ(ierr)
            Call SimplePoissonRHS_Cst(snesTemp,rhsTemp,time(TimeStepNum),MEF90Ctx,ierr)
         Case(Poisson_CST)
            Call SimplePoissonFormInitialGuess_Cst(snesTemp,solTemp,1.0_Kr,MEF90Ctx,ierr);CHKERRQ(ierr)
            Call SimplePoissonRHS_Cst(snesTemp,rhsTemp,1.0_Kr,MEF90Ctx,ierr)
         Case(Poisson_FILE)
            Call SimplePoissonLoadEXO(mesh,EXOin,bctemp,flux,reftemp,TimeStepNum,MEF90Ctx,ierr);CHKERRQ(ierr)
            Call SimplePoissonFormInitialGuess(snesTemp,solTemp,bcTemp,MEF90Ctx,ierr);CHKERRQ(ierr)
            Call SimplePoissonRHS(snesTemp,rhsTemp,flux,refTemp,MEF90Ctx,ierr)
         End Select
   
         !!! No need to solve over and over when doing CST loading
         If ((GlobalProperties%LoadingType /= Poisson_CST) .OR. (TimeStepNum == 1)) Then
            Call SNESSolve(snesTemp,rhsTemp,solTemp,ierr);CHKERRQ(ierr)
            !!! Check SNES / KSP convergence
            Call SNESGetConvergedReason(snesTemp,reasonTemp,ierr);CHKERRQ(ierr)
            Call SNESGetIterationNumber(snesTemp,itsTemp,ierr);CHKERRQ(ierr)
            Write(IOBuffer,110) itsTemp,reasonTemp
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   110      Format('SNESTemp converged in ',I4,' iterations. SNESConvergedReason is ', I4,'\n')
      
            !!! Compute energy and work
            !!!
            If ((GlobalProperties%LoadingType == Poisson_MIL) .OR. (GlobalProperties%LoadingType == Poisson_CST)) Then
               Call SimplePoissonEnergies_Cst(snesTemp,solTemp,time(timestepnum),MEF90Ctx,energy,work,ierr)
            Else
               Call SimplePoissonEnergies(snesTemp,solTemp,flux,MEF90Ctx,energy,work,ierr)
            End If
      
            !!! Print and save energy and work
            !!!
            Write(IOBuffer,*) '\n'
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
            Do set = 1,size(CellSetGlobalID)
               Write(IOBuffer,102) CellSetGlobalID(set),energy(set),work(set)
               Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
               Write(IOBuffer,202) timeStepNum,time(timeStepNum),energy(set),work(set),energy(set)-work(set)
               Call PetscViewerASCIIPrintf(energyViewerCellSet(set),IOBuffer,ierr);CHKERRQ(ierr)
            End Do
            Write(IOBuffer,103) sum(energy),sum(work)
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
            Write(IOBuffer,202) timeStepNum,time(timeStepNum),sum(energy),sum(work),sum(energy) - sum(work)
            Call PetscViewerASCIIPrintf(energyViewer,IOBuffer,ierr);CHKERRQ(ierr)
   102      Format('Cell set ',I4.4,' energy: ',ES12.5,' work: ',ES12.5,'\n')
   202      Format(I4,'  ',4(ES12.5,'  '),'\n')
   103      Format('=====================================================\n',         &
                 'Total:        energy: ',ES12.5,' work: ',ES12.5,'\n')
         End If   
         !!! Save results
         !!!
         Call SimplePoissonSaveEXO(mesh,EXOout,solTemp,TimeStepNum,time(TimeStepNum),sum(energy),sum(work),MEF90Ctx,ierr);CHKERRQ(ierr)
         Call PetscLogView(logviewer,ierr);CHKERRQ(ierr)
      Else
         Write(IOBuffer,209) time(timeStepNum)
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   209   Format("Timestepping to t=",ES12.5,"\n")
   
         Select case(GlobalProperties%LoadingType)
         Case(Poisson_MIL)
            Call SimplePoissonFormInitialGuess_Cst(snesTemp,solTemp,time(timeStepNum),MEF90Ctx,ierr);CHKERRQ(ierr)
         Case(Poisson_CST)
            Call SimplePoissonFormInitialGuess_Cst(snesTemp,solTemp,1.0_Kr,MEF90Ctx,ierr);CHKERRQ(ierr)
         Case(Poisson_FILE)
            Call SimplePoissonLoadEXO(mesh,EXOin,bctemp,flux,reftemp,TimeStepNum,MEF90Ctx,ierr);CHKERRQ(ierr)
            Call SimplePoissonFormInitialGuess(snesTemp,solTemp,bcTemp,MEF90Ctx,ierr);CHKERRQ(ierr)
         End Select
         
         Call TSSetDuration(tsTemp,TSmaxStep,time(timeStepNum),ierr);CHKERRQ(ierr)
         !Call TSView(tsTemp,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
         Call TSSolve(tsTemp,solTemp,TSsolutionTime,ierr);CHKERRQ(ierr)
         Call TSGetTimeStepNumber(tsTemp,itsTemp,ierr);CHKERRQ(ierr)
         Call TSGetConvergedReason(tsTemp,reasonTemp,ierr);CHKERRQ(ierr)
         Write(IOBuffer,210) itsTemp,reasonTemp,TSsolutionTime
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   210   Format('TSTemp converged in ',I4,' steps. TSConvergedReason is ', I4,' Analysis time is ',ES12.5,'\n')

         !!! Compute energy and work
         !!!
         If ((GlobalProperties%LoadingType == Poisson_MIL) .OR. (GlobalProperties%LoadingType == Poisson_CST)) Then
            Call SimplePoissonEnergies_Cst(snesTemp,solTemp,time(timestepnum),MEF90Ctx,energy,work,ierr)
         Else
            Call SimplePoissonEnergies(snesTemp,solTemp,flux,MEF90Ctx,energy,work,ierr)
         End If
   
         !!! Print and save energy and work
         !!!
         Write(IOBuffer,*) '\n'
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         Do set = 1,size(CellSetGlobalID)
            Write(IOBuffer,102) CellSetGlobalID(set),energy(set),work(set)
            Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
            Write(IOBuffer,202) timeStepNum,time(timeStepNum),energy(set),work(set),energy(set)-work(set)
            Call PetscViewerASCIIPrintf(energyViewerCellSet(set),IOBuffer,ierr);CHKERRQ(ierr)
         End Do
         Write(IOBuffer,103) sum(energy),sum(work)
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         Write(IOBuffer,202) timeStepNum,time(timeStepNum),sum(energy),sum(work),sum(energy) - sum(work)
         Call PetscViewerASCIIPrintf(energyViewer,IOBuffer,ierr);CHKERRQ(ierr)
         !!! Save results
         !!!
         Call SimplePoissonSaveEXO(mesh,EXOout,solTemp,TimeStepNum,TSsolutionTime,sum(energy),sum(work),MEF90Ctx,ierr);CHKERRQ(ierr)
         Call PetscLogView(logviewer,ierr);CHKERRQ(ierr)
      End If
   End Do
   !!!
   !!! Cleanup
   !!!
   Call EXCLOS(exoIN,ierr)
   Call EXCLOS(exoOUT,ierr)
   Call ISRestoreIndicesF90(CellSetGlobalIS,CellSetGlobalID,ierr);CHKERRQ(ierr)
   Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call PoissonCtxDestroy(MEF90Ctx,snesTemp,ierr);CHKERRQ(ierr)
   If ( GlobalProperties%TimeEvolution == Poisson_SteadyState) Then
      Call SNESDestroy(snesTemp,ierr);CHKERRQ(ierr)
   Else
      Call TSDestroy(tsTemp,ierr);CHKERRQ(ierr)
   End If
   Call VecDestroy(solTemp,ierr);CHKERRQ(ierr)
   Call VecDestroy(resTemp,ierr);CHKERRQ(ierr)
   Call VecDestroy(RHSTemp,ierr);CHKERRQ(ierr)
   If (GlobalProperties%LoadingType == Poisson_FILE) Then
      Call VecDestroy(flux,ierr);CHKERRQ(ierr)
      Call VecDestroy(refTemp,ierr);CHKERRQ(ierr)
      Call VecDestroy(bcTemp,ierr);CHKERRQ(ierr)
   End If
   !Call SectionRealDestroy(secTemp,ierr);CHKERRQ(ierr)
   Call DMDestroy(mesh,ierr);CHKERRQ(ierr);

   DeAllocate(work)
   DeAllocate(energy)
   
   Call PetscViewerDestroy(energyViewer,ierr);CHKERRQ(ierr)
   Do set = 1, numCellSetGlobal
      Call PetscViewerDestroy(energyViewerCellSet(set),ierr);CHKERRQ(ierr)
   End Do
   DeAllocate(energyViewerCellSet)
   Call PetscLogView(logviewer,ierr);CHKERRQ(ierr)
   Call PetscViewerDestroy(logViewer,ierr);CHKERRQ(ierr)

   If (GlobalProperties%verbose > 0) Then
      Call PetscOptionsView(PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
   EndIf
   Call MEF90_Finalize()


End Program  SimplePoissonNG
