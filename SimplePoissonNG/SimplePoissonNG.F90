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

   Type(PoissonGlobalProperties_Type),parameter    :: defaultGlobalProperties    = PoissonGlobalProperties_Type(0,PETSC_FALSE)
   Type(PoissonCellSetProperties_Type)             :: defaultCellSetProperties   = PoissonCellSetProperties_Type(DEFAULT_ELEMENT_SHORTID,0.0_Kr,0.0_Kr,0.0_Kr)
   Type(PoissonVertexSetProperties_Type),parameter :: defaultVertexSetProperties = PoissonVertexSetProperties_Type(PETSC_TRUE,0)
   Type(PoissonGlobalProperties_Type),pointer      :: GlobalProperties 
   Character(len=MEF90_MXSTRLEN)                   :: prefix,filename
   Type(MEF90Ctx_Type)                             :: MEF90Ctx
   Type(DM)                                        :: mesh,tmp_mesh
   PetscErrorCode                                  :: iErr
   PetscInt                                        :: exo_step,exo_field
   Character(len=MEF90_MXSTRLEN)                   :: IOBuffer
   PetscInt                                        :: numDim
   PetscBool                                       :: flg
   Type(SNES)                                      :: snesTemp
   Type(KSP)                                       :: kspTemp
   Type(PC)                                        :: pcTemp
   Type(Mat)                                       :: matTemp
   Type(Vec)                                       :: solTemp,resTemp,locTemp
   !Type(Vec)                                       :: ubTemp,lbTemp
   Integer                                         :: cpu_ws,io_ws,exoIN=0,exoOUT=0
   Real                                            :: exo_version
   Integer                                         :: exoerr
   MPI_Comm                                        :: IOComm
   PetscReal,Dimension(:),Pointer                  :: energy,work
   PetscInt,dimension(:),Pointer                   :: setID
   PetscInt                                        :: set
   Type(MatNullSpace)                              :: nspTemp
   SNESConvergedReason                             :: reasonTemp
   PetscInt                                        :: itsTemp
   Type(SectionReal)                               :: secTemp
   Type(VecScatter)                                :: ScatterSecToVec
   Type(IS)                                        :: CellSetGlobalIS
   Type(Element_Type),Dimension(:),pointer         :: ElemType
      
   PetscInt                                        :: point,numCell,numVertex
   PetscReal,Dimension(:,:),Pointer                :: Coord
   PetscReal,Dimension(:),Pointer                  :: Coord2
   PetscInt                                        :: i,j

   Call MEF90_Initialize()
   Call m_Poisson_Initialize(ierr);CHKERRQ(ierr)

   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-prefix',prefix,flg,ierr);CHKERRQ(ierr)
   If (.NOT. flg) Then
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_NULL,"Missing file prefix\n",ierr);CHKERRQ(ierr)
   End If

   Call MEF90CtxPoissonGlobalPropertiesCreate(MEF90Ctx,defaultGlobalProperties,ierr)
   Call PetscBagGetDataPoissonGlobalProperties(MEF90Ctx%GlobalPropertiesBag,GlobalProperties,ierr);CHKERRQ(ierr)      
   !!!
   !!! Open INPUT mesh
   !!!
   If (MEF90_MyRank == 0) Then
      cpu_ws = 8
      io_ws = 8
      filename = Trim(prefix)//'.gen'
      exoIN = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exo_version,exoerr)
      If (GlobalProperties%verbose > 0) Then
         Write(IOBuffer,99) exoERR,exoIN
         Call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr);CHKERRQ(ierr)
      End If
99 Format('EXOPEN status: ',I4,' exoIN: ',I4,'\n')         
      If (exoerr < 0) Then
         Write(IOBuffer,*) '\n\nError opening EXO file ',trim(filename),'\n\n'
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr);
         SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,IOBuffer,ierr);
      EndIf
   End If
   
   !!!
   !!! open OUTPUT mesh
   !!!
   If (GlobalProperties%splitIO) Then
      IOComm = PETSC_COMM_SELF
      Write(filename,100) trim(prefix),MEF90_MyRank
   Else
      IOComm = PETSC_COMM_WORLD
      Write(filename,101) trim(prefix)
   End If
100 Format(A,'-',I4.4,'.gen')
101 Format(A,'_out.gen')

   !!!
   !!! Read DMMesh from exoIN
   !!!
   If (MEF90_NumProcs == 1) Then
      Call DMmeshCreateExodusNG(PETSC_COMM_WORLD,exoIN,mesh,ierr);CHKERRQ(ierr)
   Else
      Call DMmeshCreateExodusNG(PETSC_COMM_WORLD,exoIN,Tmp_mesh,ierr);CHKERRQ(ierr)
   
      Call DMmeshDistribute(Tmp_mesh,PETSC_NULL_CHARACTER,mesh,ierr);CHKERRQ(ierr)
      Call DMDestroy(Tmp_mesh,ierr);CHKERRQ(ierr)
   End If

   !!! 
   !!! Create SNES and associate mesh
   !!!
   Call SNESCreate(PETSC_COMM_WORLD,snesTemp,ierr);CHKERRQ(ierr)
   Call SNESSetDM(snesTemp,mesh,ierr);CHKERRQ(ierr)
   Call SNESSetOptionsPrefix(snesTemp,'temp_',ierr);CHKERRQ(ierr)


   !!! 
   !!! Create a Section consistent with the element choice
   !!! The section is named 'default' so that it can be picked as the default
   !!! layout for all DM vector creation routines
   !!!
   Call DMmeshGetStratumSize(mesh,"depth",0,numVertex,ierr);CHKERRQ(ierr)
   Call DMmeshGetStratumSize(mesh,"height",0,numCell,ierr);CHKERRQ(ierr)
   Call DMMeshGetSectionReal(mesh,'default',SecTemp,ierr);CHKERRQ(ierr)
   Do point = numCell,numCell+numVertex-1
      Call SectionRealSetFiberDimension(SecTemp,point,1,ierr);CHKERRQ(ierr)
   End Do
   Call SectionRealAllocate(SecTemp,ierr);CHKERRQ(ierr)

   !!! 
   !!! Create PoissonCtx
   !!!
   Call MEF90CtxPoissonVertexSetPropertiesCreate(MEF90Ctx,snesTemp,defaultVertexSetProperties,ierr);CHKERRQ(ierr)
   Call EXOGetCellSetElementType_Scal(PETSC_COMM_WORLD,exoIN,MEF90_DIM,ElemType,ierr)
   Call MEF90CtxPoissonCellSetPropertiesCreate(MEF90Ctx,snesTemp,defaultCellSetProperties,ElemType,ierr);CHKERRQ(ierr)
   !!!
   !!! Get Matrix for the Jacobian / SNES and unknown vector
   !!!
   Call DMMeshSetMaxDof(mesh,1,iErr); CHKERRQ(iErr) 
   Call DMMeshCreateMatrix(mesh,secTemp,MATAIJ,matTemp,iErr);CHKERRQ(iErr)

   !!! There are two ways to get global and local vectors:
   !!! Through the generic DM interface, which will only pull the layout from the 'default' section
   !Call DMGetGlobalVector(mesh,solTemp,ierr);CHKERRQ(ierr)
   !!! Through DMMeshCreateVector which lets use any layout
   Call DMMeshCreateVector(mesh,secTemp,solTemp,ierr);CHKERRQ(ierr)

   Call VecDuplicate(solTemp,resTemp,ierr);CHKERRQ(ierr)
      
   !!! Adding a null space when some boundary conditions are prescribes breaks everything...
   !!! Need to add a flag and make adding the null space optional
   !Call MatNullSpaceCreate(PETSC_COMM_WORLD,PETSC_TRUE,0,PETSC_NULL_OBJECT,nspTemp,ierr);CHKERRQ(ierr)
   !Call MatSetNullSpace(matTemp,nspTemp,ierr);CHKERRQ(ierr)

   !!! Not sure if this is still needed when using MatZeroRowsColumnsIS for BC handling
   Call MatSetOption(matTemp,MAT_SPD,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_SYMMETRY_ETERNAL,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetOption(matTemp,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE,ierr);CHKERRQ(ierr)
   Call MatSetFromOptions(matTemp,ierr);CHKERRQ(ierr)

   !!!
   !!! Set Jacobian and Function for the SNES
   !!!
   Call SNESSetFunction(snesTemp,resTemp,SimplePoissonOperator,MEF90Ctx,ierr);CHKERRQ(ierr)
   Call SNESSetJacobian(snesTemp,matTemp,matTemp,SimplePoissonBilinearForm,MEF90Ctx,ierr);CHKERRQ(ierr)

   !!!
   !!! Testing SNESVI: does not look ready for prime time...
   !!!
   !Call VecDuplicate(solTemp,lbTemp,ierr);CHKERRQ(ierr)
   !Call VecSet(lbTemp,-1.0_Kr,ierr);CHKERRQ(ierr)
   !Call VecDuplicate(solTemp,ubTemp,ierr);CHKERRQ(ierr)
   !Call VecSet(ubTemp,0.2_Kr,ierr);CHKERRQ(ierr)
   !Call SNESVISetVariableBounds(snesTemp,lbTemp,ubTemp,ierr);CHKERRQ(ierr)

   Call SNESSetFromOptions(snesTemp,ierr);CHKERRQ(ierr)

   Call SNESGetKSP(snesTemp,kspTemp,ierr);CHKERRQ(ierr)
   Call KSPSetType(kspTemp,KSPCG,ierr);CHKERRQ(ierr)
   Call KSPSetInitialGuessNonzero(kspTemp,PETSC_TRUE,ierr);CHKERRQ(ierr)
!!!   Call KSPSetNullSpace(kspTemp,nspTemp,ierr);CHKERRQ(ierr)
   Call KSPSetFromOptions(kspTemp,ierr);CHKERRQ(ierr)

   Call KSPGetPC(kspTemp,pcTemp,ierr);CHKERRQ(ierr)
   Call PCSetFromOptions(pcTemp,ierr);CHKERRQ(ierr)
   !!! Setup GAMG here (coordinates,in particular)
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

   !!! Solve Poisson Equation
   Call SimplePoissonFormInitialGuess(snesTemp,solTemp,MEF90Ctx,ierr);CHKERRQ(ierr)
   Call SNESSolve(snesTemp,PETSC_NULL_OBJECT,solTemp,ierr);CHKERRQ(ierr)

   !!! Check SNES / KSP convergence
   Call SNESGetConvergedReason(snesTemp,reasonTemp,ierr);CHKERRQ(ierr)
   Call SNESGetIterationNumber(snesTemp,itsTemp,ierr);CHKERRQ(ierr)

   Write(IOBuffer,110) itsTemp,reasonTemp
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
110 Format('SNESTemp converged in in ',I4,' iterations. SNESConvergedReason is ', I4,'\n')
   
   !!! Compute energy and work
   !!! This is one of the few looks that need to be synchronized across processors!
   !!!
   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS,ierr);CHKERRQ(ierr)   
   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Allocate(energy(size(setID)),stat=ierr)
   Allocate(work(size(setID)),stat=ierr)

   Call SimplePoissonEnergies(snesTemp,solTemp,MEF90Ctx,energy,work,ierr)
   Write(IOBuffer,*) '\n'
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   Do set = 1,size(setID)
      Write(IOBuffer,102) setID(set),energy(set),work(set)
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
   End Do
   Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Call ISDestroy(CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Write(IOBuffer,103) sum(energy),sum(work)
   Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
102 Format('Cell set ',I4.4,' energy: ',ES12.5,' work: ',ES12.5,'\n')
103 Format('=====================================================\n',         &
           'Total:        energy: ',ES12.5,' work: ',ES12.5,'\n')

   !!! Save results
   !!!
   !!! format the output file
   !!!
   If (GlobalProperties%splitIO) Then
      cpu_ws = 8
      io_ws = 8
      exoOUT = EXCRE(trim(filename),EXCLOB,cpu_ws,io_ws,ierr)
      Call DMmeshViewExodusSplit(mesh,exoOUT,ierr)
      Call EXPVP (exoOUT,'g',3,ierr)
      Call EXPVAN(exoOUT,'g',3,(/'Energy         ','Ext Fluxes work','Total Energy   '/),ierr)
      Call EXPVP (exoOUT,'n',1,ierr)
      Call EXPVAN(exoOUT,'n',1,(/'U'/),ierr)
      !Call EXPVAN(exoOUT,'e',1,(/'F'/),ierr)
   Else
      If (MEF90_MyRank == 0) Then
         cpu_ws = 8
         io_ws = 8
         exoOUT = EXCRE(filename,EXCLOB,cpu_ws,io_ws,ierr)
         Call EXCOPY(exoIN,exoOUT,ierr)
         Call EXPVP (exoOUT,'g',3,ierr)
         Call EXPVAN(exoOUT,'g',3,(/'Elastic Energy ','Ext Fluxes work','Total Energy   '/),ierr)
         Call EXPVP (exoOUT,'n',1,ierr)
         Call EXPVAN(exoOUT,'n',1,(/'U'/),ierr)
         !Call EXPVAN(exoOUT,'e',1,(/'F'/),ierr)
      End If
   End If

   !!! 
   !!! Write in exo file
   !!!
   !!! Just as for global vectors, there are two ways to deal with local vectors
   !!! Obtaining them from the DM, and using DMLocalToGlobalBegin/End
   !Call DMGetLocalVector(mesh,locTemp,ierr);CHKERRQ(ierr)
   !Call DMGlobalToLocalBegin(mesh,solTemp,INSERT_VALUES,locTemp,ierr);CHKERRQ(ierr)
   !Call DMGlobalToLocalEnd(mesh,solTemp,INSERT_VALUES,locTemp,ierr);CHKERRQ(ierr)
   !!! Or by from a SectionReal, and using SectionRealToVec to copy from Vec to Section
   Call DMMeshCreateGlobalScatter(mesh,secTemp,ScatterSecToVec,ierr);CHKERRQ(ierr)
   Call SectionRealCreateLocalVector(secTemp,locTemp,ierr);CHKERRQ(ierr)   
   Call SectionRealToVec(secTemp,ScatterSecToVec,SCATTER_REVERSE,solTemp,ierr);CHKERRQ(ierr)
   Call VecScatterDestroy(ScatterSecToVec,ierr);CHKERRQ(ierr)

   Call VecViewExodusVertex(mesh,locTemp,IOComm,exoOUT,1,1,ierr)
   If ( (GlobalProperties%splitIO) .OR. (MEF90_MyRank == 0)) Then
      Call EXPGV(exoOUT,1,3,(/ sum(energy),sum(work),sum(energy)-sum(work) /),ierr)   
   End If
   Call DMrestoreLocalVector(mesh,locTemp,ierr);CHKERRQ(ierr)

   !!!
   !!! Cleanup
   !!!
   Call EXCLOS(exoIN,ierr)
   Call EXCLOS(exoOUT,ierr)
   Call PoissonCtxDestroy(MEF90Ctx,snesTemp,ierr);CHKERRQ(ierr)
   Call SNESDestroy(snesTemp,ierr);CHKERRQ(ierr)
   Call VecDestroy(solTemp,ierr);CHKERRQ(ierr)
   Call VecDestroy(resTemp,ierr);CHKERRQ(ierr)
   Call VecDestroy(locTemp,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(secTemp,ierr);CHKERRQ(ierr)
   Call DMDestroy(mesh,ierr);CHKERRQ(ierr);

   DeAllocate(work)
   DeAllocate(energy)
   
   If (GlobalProperties%verbose > 0) Then
      Call PetscOptionsView(PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
   EndIf
   Call MEF90_Finalize()


End Program  SimplePoissonNG
