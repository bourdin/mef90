#if defined PB_2D
#define M_POISSON m_Poisson2D
#define ELEMENT_SCAL Element2D_Scal
#define VECT Vect2D
#elif defined PB_3D
#define M_POISSON m_Poisson3D
#define ELEMENT_SCAL Element3D_Scal
#define VECT Vect3D
#endif

Module M_POISSON

#include "finclude/petscdef.h"
   Use m_MEF90
   Implicit NONE   

   Type LogInfo_Type
      PetscLogStage                                :: IO_Stage
      PetscLogStage                                :: Setup_Stage 
      PetscLogStage                                :: meshCreateExodus_Stage
      PetscLogStage                                :: meshDistribute_Stage
      PetscLogStage                                :: MatAssembly_Stage    
      PetscLogStage                                :: RHSAssembly_Stage
      PetscLogStage                                :: KSPSolve_Stage
      PetscLogStage                                :: PostProc_Stage
      
      PetscLogEvent                                :: MatAssemblyBlock_Event
      PetscLogEvent                                :: RHSAssemblyBlock_Event
      PetscLogEvent                                :: PostProc_Event
   End Type LogInfo_Type

   Type AppParam_Type
      PetscBool                                    :: Restart
      PetscBool                                    :: splitIO
      PetscInt                                     :: Verbose
      PetscInt                                     :: TestCase
      Character(len=MEF90_MXSTRLEN)                :: prefix
      Type(PetscViewer)                            :: LogViewer,MyLogViewer
   End Type AppParam_Type


   Type Poisson_AppCtx_Type
      Type(dm)                                     :: mesh
      Type(SNES)                                   :: snesU
      Type(EXO_Type)                               :: EXO
      Type(ELEMENT_SCAL),Dimension(:),Pointer      :: Elem
      Type(IS)                                     :: CellSetGlobalIS
      Type(IS)                                     :: VertexSetGlobalIS
      Type(Element_Type),Dimension(:),Pointer      :: ElementType
      !!! Since everything is done block after block, it is simpler to have one set of elements per block
      Type(SectionReal)                            :: GradU
      PetscReal                                    :: ElasticEnergy
      Type(Field)                                  :: U
      Type(Field)                                  :: F
      PetscReal                                    :: ExtForcesWork
      PetscReal                                    :: TotalEnergy
      Type(Flag)                                   :: BCFlag
      Type(Mat)                                    :: K
      Type(Field)                                  :: RHS

      Type(KSP)                                    :: KSPU
      Type(PC)                                     :: PCU
      ! Should be pulled from the SNES?
      Type(LogInfo_Type)                           :: LogInfo
      Type(AppParam_Type)                          :: AppParam
   !For TS
      Type(TS)                                     :: TS
      PetscInt                                     :: maxsteps,NumSteps
      PetscReal                                    :: maxtime
      PetscInt                                     :: VertVar_Temperature 
      PetscReal,Dimension(:),Pointer               :: Diff,B_Mensi
   End Type Poisson_AppCtx_Type
   
   
Contains
#undef __FUNCT__
#define __FUNCT__ "SimplePoissonInit"
   Subroutine SimplePoissonInit(AppCtx)
      Type(Poisson_AppCtx_Type)                    :: AppCtx
      
      PetscInt                                     :: ierr
      PetscInt                                     :: iBlk,iDoF      
      PetscBool                                    :: HasPrefix,Flag
      PetscInt,Dimension(:),Pointer                :: TmpFlag
      
      PetscReal,Dimension(:),Pointer               :: ValPtr
      PetscReal,Dimension(:),Pointer               :: Coord
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer,filename   
      Type(DM)                                     :: Tmp_mesh
      PetscReal                                    :: Val
      PetscInt,Dimension(:),Pointer                :: SizeScal
      PetscInt                                     :: numCell,numDim,c
      Type(IS)                                     :: setIS,vertexIS,cellIS
      PetscInt,Dimension(:),Pointer                :: setID,vertexID
      PetscInt                                     :: set,vertex
      Type(SectionReal)                            :: coordSection
      Integer                                      :: cpu_ws,io_ws,exoidIN
      Real                                         :: exo_version
      PetscInt                                     :: numSizes
      PetscInt,Dimension(:),Pointer                :: elemTypes
      PetscBool                                    :: flg

      Call MEF90_Initialize()
      AppCtx%AppParam%verbose = 0
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-verbose',AppCtx%AppParam%verbose,Flag,ierr);CHKERRQ(ierr)
      
      AppCtx%AppParam%Restart = PETSC_FALSE
      Call PetscOptionsGetBool(PETSC_NULL_CHARACTER,'-restart',AppCtx%AppParam%restart,Flag,ierr);CHKERRQ(ierr)
      
      AppCtx%AppParam%splitIO = PETSC_FALSE
      Call PetscOptionsGetBool(PETSC_NULL_CHARACTER,'-splitIO',AppCtx%AppParam%splitIO,Flag,ierr);CHKERRQ(ierr)
      
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-p',      AppCtx%AppParam%prefix,HasPrefix,ierr);CHKERRQ(ierr)
      If (.NOT. HasPrefix) Then
         Call PetscPrintf(PETSC_COMM_WORLD,"No mesh prefix given\n",ierr)
         Call MEF90_Finalize()
         STOP
      End If

      AppCtx%AppParam%TestCase = 1
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-test',      AppCtx%AppParam%TestCase,Flag,ierr);CHKERRQ(ierr)
      
      Call InitLog(AppCtx)
      If (AppCtx%AppParam%verbose > 1) Then
         Write(filename,101) Trim(AppCtx%AppParam%prefix),MEF90_MyRank
         Call PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,AppCtx%AppParam%MyLogViewer,ierr);CHKERRQ(ierr);  
         Write(IOBuffer,102) MEF90_MyRank,Trim(filename)
         Call PetscSynchronizedPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         Call PetscSynchronizedFlush (PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
   
         Write(filename,103) Trim(AppCtx%AppParam%prefix)
         Call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,AppCtx%AppParam%LogViewer,ierr);CHKERRQ(ierr);  
         Write(IOBuffer,104) Trim(filename)
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      End If
   
101 Format(A,'-',I4.4,'.log')
102 Format('Output from processor ',I4.4,' redirected to file ',A,'\n')
103 Format(A,'.log')
104 Format('Collective output redirected to file ',A,'\n')


      !!! Read mesh from <prefix>.gen
      If (MEF90_MyRank == 0) Then
         cpu_ws = 8
         io_ws = 8
         filename = Trim(AppCtx%AppParam%prefix)//'.gen'
         exoidIN = EXOPEN(filename,EXREAD,cpu_ws,io_ws,exo_version,ierr)
      End If
      
      !!! Read the mesh from <prefix.gen> and partition it
      If (MEF90_NumProcs == 1) Then
         Call DMmeshCreateExodusNG(PETSC_COMM_WORLD,exoidIN,AppCtx%mesh,ierr);CHKERRQ(ierr)
      Else
         Call DMmeshCreateExodusNG(PETSC_COMM_WORLD,exoidIN,Tmp_mesh,ierr);CHKERRQ(ierr)
      
         Call DMmeshDistribute(Tmp_mesh,PETSC_NULL_CHARACTER,AppCtx%mesh,ierr);CHKERRQ(ierr)
         Call DMDestroy(Tmp_mesh,ierr);CHKERRQ(ierr)
      End If
      
      !!! Prepare the output file(s)
      If (AppCtx%AppParam%splitIO) Then
         AppCtx%EXO%Comm = PETSC_COMM_SELF
         Write(AppCtx%EXO%filename,105) trim(AppCtx%AppParam%prefix),MEF90_MyRank
      Else
         AppCtx%EXO%Comm = PETSC_COMM_WORLD
         Write(AppCtx%EXO%filename,106) trim(AppCtx%AppParam%prefix)
      End If
      
      !!! Open the output file, save geometry and format the file if restart is false
      If (AppCtx%AppParam%restart) Then
         !!! Open existing result file. 
         If (AppCtx%AppParam%splitIO .OR. (MEF90_MyRank == 0)) Then
            AppCtx%EXO%exoid = EXOPEN(AppCtx%EXO%filename,EXWRIT,cpu_ws,io_ws,exo_version,ierr)
         End If
      Else
         If (AppCtx%AppParam%splitIO) Then
            cpu_ws = 8
            io_ws = 8
            AppCtx%EXO%exoid = EXCRE(trim(AppCtx%EXO%filename),EXCLOB,cpu_ws,io_ws,ierr)
            Call DMmeshViewExodusSplit(AppCtx%mesh,AppCtx%EXO%exoid,ierr)
         Else
            If (MEF90_MyRank == 0) Then
               cpu_ws = 8
               io_ws = 8
               AppCtx%EXO%exoid = EXCRE(AppCtx%EXO%filename,EXCLOB,cpu_ws,io_ws,ierr)
               Call EXCOPY(exoidIN,AppCtx%EXO%exoid,ierr)
            End If
         End If
         Call EXOFormat_SimplePoisson(AppCtx%EXO)
      End If
      !!! Close the geometry file
      If (MEF90_MyRank == 0) Then
         Call EXCLOS(exoidIN,ierr)
      End If
105 Format(A,'-',I4.4,'.gen')
106 Format(A,'_out.gen')
      !!!
      !!! Compute global IS for cell and vertex sets
      !!!
      Call DMmeshGetLabelIdIS(AppCtx%mesh,'Cell Sets',AppCtx%CellSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,AppCtx%CellSetGlobalIS)
      Call DMmeshGetLabelIdIS(AppCtx%mesh,'Vertex Sets',AppCtx%VertexSetGlobalIS,ierr);CHKERRQ(ierr)
      Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,AppCtx%VertexSetGlobalIS)
      
      !!!
      !!! Allocate Element array
      !!!
      Call DMmeshGetStratumSize(AppCtx%mesh,"height",0,numCell,ierr);CHKERRQ(ierr)
      Allocate(AppCtx%Elem(numCell))

      !!! 
      !!! Set element type to P1 Lagrange on all sets
      !!!
      Call DMmeshGetDimension(AppCtx%mesh,numDim,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(AppCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
      Allocate(AppCtx%ElementType(size(setID)))
      Do set = 1,size(setID)
         If (numDim == 2) Then
            AppCtx%ElementType(setID(set)) = MEF90_P1_Lagrange_2D_Scal
         Else
            AppCtx%ElementType(setID(set)) = MEF90_P1_Lagrange_3D_Scal
         End If
      End Do
      Allocate(elemTypes(size(setID)))
      elemTypes = -2
      numSizes = size(setID)
      Call PetscOptionsGetIntArray(PETSC_NULL_CHARACTER,"-elem_type",elemTypes,numSizes,flg,ierr);CHKERRQ(ierr)
      write(*,*) numSizes,elemTypes
      Do set = 1, numSizes
         Call Element_TypeFindByID(elemTypes(set),AppCtx%ElementType(set))
         Write(*,*) "Changed element type for set ", set, "to ",trim(AppCtx%ElementType(set)%name)
      End Do      
      DeAllocate(elemTypes)      
      Do set = 1,size(setID)
         Call DMmeshGetStratumIS(AppCtx%mesh,'Cell Sets',setID(set),cellIS,ierr); CHKERRQ(ierr)
         Call ElementInit(AppCtx%mesh,cellIS,AppCtx%Elem,2,AppCtx%ElementType(setID(set)))
         Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
      End Do    
      Call ISRestoreIndicesF90(AppCtx%CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   
      !!! Allocate the Section for U and F
      Allocate(SizeScal(1))
      SizeScal=1

      Call FieldCreateVertex(AppCtx%U,    'U',        AppCtx%mesh,SizeScal)
      Call FieldDestroy(AppCtx%U)
      Call FieldCreateVertex(AppCtx%U,    'U',        AppCtx%mesh,SizeScal)
      Call FieldCreateVertex(AppCtx%F,    'F',        AppCtx%mesh,SizeScal)
      Call FieldCreateVertex(AppCtx%RHS,  'RHS',      AppCtx%mesh,SizeScal)


      !!! Allocate and initialize the Section for the flag
      Call FlagCreateVertex(AppCtx%BCFlag,'BC',AppCtx%mesh,SizeScal)      
      Allocate(TmpFlag(1))
      TmpFlag = 1
      Call DMmeshGetLabelIdIS(AppCtx%mesh,'Vertex Sets',setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         Call DMmeshGetStratumIS(AppCtx%mesh,'Vertex Sets',setID(set),vertexIS,ierr); CHKERRQ(ierr)
         Call ISGetIndicesF90(vertexIS,vertexID,ierr);CHKERRQ(ierr)
         Do vertex = 1, size(vertexID)
            Call SectionIntUpdate(AppCtx%BCFlag%Sec,vertexID(vertex),TmpFlag,INSERT_VALUES,ierr);CHKERRQ(ierr)
         End Do
         Call ISRestoreIndicesF90(vertexIS,vertexID,ierr);CHKERRQ(ierr)
         Call ISDestroy(vertexIS,ierr);CHKERRQ(ierr)
      End Do
      DeAllocate(TmpFlag)
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)


      !!! Initialize the matrix and vector for the linear system
      Call DMmeshSetMaxDof(AppCtx%mesh,1,ierr);CHKERRQ(ierr) 
      Call DMmeshCreateMatrix(AppCtx%mesh,AppCtx%U%Sec,MATMPIAIJ,AppCtx%K,ierr);CHKERRQ(ierr)

      !!! Create the SNEs from which we will get the KSP and PC
      !Call SNESCreate(PETSC_COMM_WORLD,AppCtx%SNESU,ierr);CHKERRQ(ierr)
      !Call SNESSetJacobian(AppCtx%SNESU,AppCtx%K,AppCtx%K,PoissonMatAssembly,AppCtx,ierr);CHKERRQ(ierr)
      !Call SNESSetType(AppCtx%SNESU,SNESKSPONLY,ierr);CHKERRQ(ierr)
      !Call SNESGetKSP(AppCtx%SNESU,AppCtx%KSPU,ierr);CHKERRQ(ierr)      
      !Call SNESAppendOptionsPrefix(AppCtx%SNESU,"U_",ierr);CHKERRQ(ierr)
      !Call SNESSetFromOptions(AppCtx%SNESU,ierr);CHKERRQ(ierr)
            
      Call KSPCreate(PETSC_COMM_WORLD,AppCtx%KSPU,ierr);CHKERRQ(ierr)
      Call KSPAppendOptionsPrefix(AppCtx%KSPU,"U_",ierr);CHKERRQ(ierr)
      Call KSPSetOperators(AppCtx%KSPU,AppCtx%K,AppCtx%K,SAME_NONZERO_PATTERN,ierr);CHKERRQ(ierr)
      Call KSPSetInitialGuessNonzero(AppCtx%KSPU,PETSC_TRUE,ierr);CHKERRQ(ierr)
      Call KSPSetType(AppCtx%KSPU,KSPCG,ierr);CHKERRQ(ierr)

      Call KSPGetPC(AppCtx%KSPU,AppCtx%PCU,ierr);CHKERRQ(ierr)
      Call PCSetType(AppCtx%PCU,PCBJACOBI,ierr);CHKERRQ(ierr)
      Call KSPSetFromOptions(AppCtx%KSPU,ierr);CHKERRQ(ierr)
      
      
      
      !!! Read Force and BC from Data file or reformat it
      
      !!!
      !!! Create local vectors, then call 
      !!! VecLoadExodusVertex(dm,res,IOcomm,exoidout,step,offset)
      !!!
      If (AppCtx%AppParam%Restart) Then
         Call VecLoadExodusVertex(AppCtx%mesh,AppCtx%U%LocalVec,PETSC_COMM_WORLD,AppCtx%EXO%exoid,1,1,ierr);CHKERRQ(ierr)
         Call VecLoadExodusVertex(AppCtx%mesh,AppCtx%F%LocalVec,PETSC_COMM_WORLD,AppCtx%EXO%exoid,1,2,ierr);CHKERRQ(ierr)
      Else
         Select Case (AppCtx%AppParam%TestCase)
         Case(1)
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer,*) 'Setting U to 0 and F to 1\n'
               Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
            End If
            !!! U = 0,F=1
            Val = 1.0_Kr
            Call SectionRealSet(AppCtx%F,Val,ierr);CHKERRQ(ierr);
            Val = 0.0_Kr
            Call SectionRealSet(AppCtx%U,Val,ierr);CHKERRQ(ierr);
         Case(2)
            If (AppCtx%AppParam%verbose > 0) Then
               Write(IOBuffer,*) 'Solving Test Case 2: pure dirichlet problem,no force\n'
               Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
            End If

            Allocate(ValPtr(1))
            Call DMmeshGetDimension(AppCtx%mesh,numDim,ierr);CHKERRQ(ierr)
            Call DMmeshGetSectionReal(AppCtx%mesh,'coordinates',coordSection,ierr); CHKERRQ(ierr)
            Call DMmeshGetStratumIS(AppCtx%mesh,'height',0,vertexIS,ierr);CHKERRQ(ierr)
            Call ISGetIndicesF90(vertexIS,vertexID,ierr);CHKERRQ(ierr)
            Do vertex = 1, size(vertexID)
               Call SectionRealRestrict(coordSection,vertexID(vertex),coord,ierr);CHKERRQ(ierr)
               ValPtr = 0
               Do c = 1, numDim
                  ValPtr = ValPtr + coord(c)**2
               End Do
               Call SectionRealUpdate(AppCtx%U%Sec,vertexID(vertex),ValPtr,INSERT_VALUES,ierr);CHKERRQ(ierr)
               Call SectionRealRestore(coordSection,vertexID(vertex),coord,ierr);CHKERRQ(ierr)
            End Do
            Call ISRestoreIndicesF90(vertexIS,vertexID,ierr);CHKERRQ(ierr)
            Call ISDestroy(vertexIS,ierr);CHKERRQ(ierr)
            Call SectionRealDestroy(coordSection,ierr);CHKERRQ(ierr)
            DeAllocate(ValPtr)
         End Select            
      End If
      
      !!! Create the EXO case file
      Call Write_EXO_Case(AppCtx%AppParam%prefix,'%0.4d',MEF90_NumProcs)
   End Subroutine SimplePoissonInit


   
#undef __FUNCT__
#define __FUNCT__ "ExoFormat_SimplePoisson"
   Subroutine EXOFormat_SimplePoisson(EXO)
      Type(EXO_Type)                               :: EXO
      PetscInt                                     :: ierr
   
      EXO%num_nsproperties = 0
      EXO%num_ebproperties = 0
      EXO%num_globvariables = 3
      Call EXPVP (EXO%exoid,'g',3,ierr)
      Call EXPVAN(EXO%exoid,'g',3,(/'Elastic Energy ','Ext Forces work','Total Energy   '/),ierr)
      Call EXPVP (EXO%exoid,'n',2,ierr)
      Call EXPVAN(EXO%exoid,'n',2,(/'U','F'/),ierr)
      EXO%num_vertvariables = 2
#if defined PB_2D
      Call EXPVP (EXO%exoid,'e',2,ierr)
      Call EXPVAN(EXO%exoid,'e',2,(/'Grad U_X','Grad U_Y'/),ierr)
      EXO%num_cellvariables = 2
#elif defined PB_3D
      Call EXPVP (EXO%exoid,'e',3,ierr)
      Call EXPVAN(EXO%exoid,'e',3,(/'Grad U_X','Grad U_Y','Grad U_Z'/),ierr)
      EXO%num_cellvariables = 3
#endif
      Call EXPTIM(EXO%exoid,1,1.0_Kr,ierr)
   End Subroutine EXOFormat_SimplePoisson
   
#undef __FUNCT__
#define __FUNCT__ "InitLog"
   Subroutine InitLog(AppCtx)
      Type(Poisson_AppCtx_Type)                    :: AppCtx
      PetscInt                                     :: ierr
      
      Call PetscLogEventRegister('MatAssembly Block',0,AppCtx%LogInfo%MatAssemblyBlock_Event,ierr);CHKERRQ(ierr)
      Call PetscLogEventRegister('RHSAssembly Block',0,AppCtx%LogInfo%RHSAssemblyBlock_Event,ierr);CHKERRQ(ierr)
      Call PetscLogEventRegister('Post Processing',  0,AppCtx%LogInfo%PostProc_Event,        ierr);CHKERRQ(ierr)

      Call PetscLogStageRegister("meshCreateExodus",AppCtx%LogInfo%meshCreateExodus_Stage,ierr)
      Call PetscLogStageRegister("meshDistribute",  AppCtx%LogInfo%meshDistribute_Stage,  ierr)
      Call PetscLogStageRegister("IO Stage",        AppCtx%LogInfo%IO_Stage,              ierr)
      Call PetscLogStageRegister("Setup",           AppCtx%LogInfo%Setup_Stage,           ierr)
      Call PetscLogStageRegister("Mat Assembly",    AppCtx%LogInfo%MatAssembly_Stage,     ierr)
      Call PetscLogStageRegister("RHS Assembly",    AppCtx%LogInfo%RHSAssembly_Stage,     ierr)
      Call PetscLogStageRegister("KSP Solve",       AppCtx%LogInfo%KSPSolve_Stage,        ierr)
      Call PetscLogStageRegister("Post Proc",       AppCtx%LogInfo%PostProc_Stage,        ierr)
   End Subroutine InitLog
   
#undef __FUNCT__
#define __FUNCT__ "Solve"
   Subroutine Solve(AppCtx)
      Type(Poisson_AppCtx_Type)                       :: AppCtx
      
      PetscInt                                     :: ierr
      KSPConvergedReason                           :: KSPreason
      PetscInt                                     :: KSPNumIter
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
      
      Call SectionRealToVec(AppCtx%U%Sec,AppCtx%U%Scatter,SCATTER_FORWARD,AppCtx%U%Vec,ierr);CHKERRQ(ierr)

      Call DMmeshCreateVector(AppCtx%mesh,AppCtx%RHS,AppCtx%RHS%Vec,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(AppCtx%RHS%Sec,AppCtx%RHS%Scatter,SCATTER_FORWARD,AppCtx%RHS%Vec,ierr);CHKERRQ(ierr)

      Call KSPSolve(AppCtx%KSPU,AppCtx%RHS%Vec,AppCtx%U%Vec,ierr);CHKERRQ(ierr)
         
      Call KSPGetIterationNumber(AppCtx%KSPU,KSPNumIter,ierr);CHKERRQ(ierr)
      Call KSPGetConvergedReason(AppCtx%KSPU,KSPreason,ierr);CHKERRQ(ierr)
      Write(IOBuffer,100) KSPNumIter,KSPreason
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Call SectionRealToVec(AppCtx%U%Sec,AppCtx%U%Scatter,SCATTER_REVERSE,AppCtx%U%Vec,ierr);CHKERRQ(ierr)
 
100 Format('KSP converged in ',I5,' iterations. KSPConvergedReason is ',I2,'\n')
   End Subroutine Solve
   
 
#undef __FUNCT__
#define __FUNCT__ "PoissonMatAssembly"
!!! Make interface compatible with SNES
   Subroutine PoissonMatAssembly(PoissonSNES,X,K,Kpre,mStruct,AppCtx,ierr)   
      Type(SNES),intent(IN)                        :: PoissonSNES
      Type(Vec),intent(IN)                         :: X
      Type(Mat),intent(IN)                         :: K,Kpre
      MatStructure,intent(IN)                      :: mStruct
      Type(Poisson_AppCtx_Type),Intent(IN)         :: AppCtx
      PetscInt,Intent(OUT)                         :: ierr
      
      PetscInt                                     :: iBlk
      Type(IS)                                     :: setIS
      PetscInt,Dimension(:),Pointer                :: setID
      PetscInt                                     :: set
      
      Call MatInsertVertexBoundaryValues(K,AppCtx%U,AppCtx%BCFlag,AppCtx%mesh)
      Call MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyEnd  (K,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)

      Call DMmeshGetLabelIdIS(AppCtx%mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         !!! Test if block corresponds to a volume or surface element
         If (AppCtx%ElementType(set)%codim == 0) Then
            Call MatAssemblyBlock(K,setID(set),AppCtx)
         End If
      End Do
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      Call MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
      Call MatAssemblyEnd  (K,MAT_FINAL_ASSEMBLY,ierr);CHKERRQ(ierr)
   End Subroutine PoissonMatAssembly
      
   Subroutine MatAssemblyBlock(K,iBlk,AppCtx)
      Type(Mat),Intent(IN)                         :: K
      PetscInt,Intent(IN)                          :: iBlk
      Type(Poisson_AppCtx_Type),Intent(IN)         :: AppCtx
      
      Type(IS)                                     :: cellIS
      PetscInt,Dimension(:),Pointer                :: cellID
      PetscInt                                     :: cell
      PetscInt                                     :: numDof
      PetscInt                                     :: ierr,i
      PetscReal,Dimension(:,:),Pointer             :: MatElem
      PetscInt,Dimension(:),Pointer                :: BCFlag
      PetscInt                                     :: iDoF1,iDoF2,iGauss
      PetscLogDouble                               :: flops = 0
      PetscReal                                    :: lDiff
      PetscReal,Dimension(:),Pointer               :: T_Loc
      PetscReal                                    :: T_Elem
     
      !!! iBlk needs to be the GLOBAL set ID so that I can recover the proper material properties
      !!! get Diff as a PETScOption
      !!! 
      numDof = AppCtx%ElementType(iBlk)%numDof
      Allocate(MatElem(numDof,numDof))
      Allocate(BCFlag(numDof))
      Allocate(T_Loc(numDof))

      Call DMmeshGetStratumIS(AppCtx%mesh,'Cell Sets',iBlk,cellIS,ierr); CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         MatElem = 0.0_Kr
         BCFlag = 0
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec,AppCtx%mesh,cellID(cell),numDof,BCFlag,ierr);CHKERRQ(ierr)
         Do iGauss = 1,size(AppCtx%Elem(cellID(cell)+1)%Gauss_C)
            Select Case(AppCtx%AppParam%TestCase)
            Case(2)
               lDiff = AppCtx%Diff(iBlk) 
               !!! THIS ASSUMES THAT BLOCKS ARE NUMBERED SEQUENTIALLY
               !!! Use PetscOptions to def lDiff, instead
            Case(3)
               ! Mensi Law D(C) = A exp (B*C) 
               ! 1988 Mensi-Acker-Attolou Mater Struct
               ! C : water content, A and B material parameters
               Call SectionRealRestrictClosure(AppCtx%U%Sec,AppCtx%mesh,cellID(cell),numDoF,T_Loc,ierr);CHKERRQ(ierr)
               T_Elem = 0.0_Kr
               Do iDoF1 = 1,numDof
                  T_Elem = T_Elem + AppCtx%Elem(cellID(cell)+1)%BF(iDoF1,iGauss) * T_Loc(iDoF1)
               End DO
               lDiff = AppCtx%Diff(iBlk)*exp(AppCtx%B_Mensi(iBlk)*T_Elem) 
               !!! SAME AS ABOVE
            Case Default
               ldiff =1 
            End Select
            Do iDoF1 = 1,numDof
               If (BCFlag(iDoF1) == 0) Then
                  Do iDoF2 = 1,numDof
                    ! MatElem(iDoF1,iDoF1) = 1./2.
                     MatElem(iDoF2,iDoF1) = MatElem(iDoF2,iDoF1) + lDiff * AppCtx%Elem(cellID(cell)+1)%Gauss_C(iGauss) * &
                                           (AppCtx%Elem(cellID(cell)+1)%Grad_BF(iDoF1,iGauss) .DotP. AppCtx%Elem(cellID(cell)+1)%Grad_BF(iDoF2,iGauss) )
                  End Do
                  flops = flops + 3 * numDof
               End If
            End Do
         End Do
         Call DMmeshAssembleMatrix(K,AppCtx%mesh,AppCtx%U%Sec,cellID(cell),MatElem,ADD_VALUES,ierr);CHKERRQ(ierr)
      End Do
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
   
      Call PetscLogFlops(flops,ierr);CHKERRQ(ierr)
      DeAllocate(BCFlag)
      DeAllocate(MatElem)
   End Subroutine MatAssemblyBlock


#undef __FUNCT__
#define __FUNCT__ "RHSAssembly"
!!! Change interface.
   Subroutine RHSAssembly(SNESU,RHSVec,UVec,AppCtx)
      Type(SNES),intent(IN)                        :: SNESU
      Type(Vec),intent(IN)                         :: RHSVec
      !!! Technically, RHS is an INOUT arg, but the value of RHSVec
      !!! should not be modified by this function
      Type(Vec),intent(IN)                         :: UVec
      Type(Poisson_AppCtx_Type)                    :: AppCtx

      PetscInt                                     :: ierr
      PetscInt                                     :: set
      PetscInt,Dimension(:),Pointer                :: setID
      Type(IS)                                     :: setIS

      Call SectionRealZero(AppCtx%RHS%Sec,ierr);CHKERRQ(ierr)
      
      Call DMmeshGetLabelIdIS(AppCtx%mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         If (AppCtx%ElementType(set)%codim > 0) Then
            !!! Only surface forces
            Call RHSAssemblyBlock(AppCtx%RHS,setID(set),AppCtx)
         End If
      End Do ! set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)

      Call SectionRealComplete(AppCtx%RHS%Sec,ierr);CHKERRQ(ierr)

      Call FieldInsertVertexBoundaryValues(AppCtx%RHS,AppCtx%U,AppCtx%BCFlag,AppCtx%mesh)
      Call SectionRealToVec(AppCtx%RHS%Sec,AppCtx%RHS%Scatter,SCATTER_FORWARD,RHSVec,ierr);CHKERRQ(ierr)
      !Call VecCopy(AppCtx%RHS%Vec,RHSVec,ierr);CHKERRQ(ierr)
   End Subroutine RHSAssembly

#undef __FUNCT__
#define __FUNCT__ "RHSAssemblyBlock"
   Subroutine RHSAssemblyBlock(RHS,setID,AppCtx)
      Type(Field),Intent(IN)                       :: RHS
      PetscInt,Intent(IN)                          :: setID
      Type(Poisson_AppCtx_Type),Intent(IN)         :: AppCtx

      PetscInt                                     :: cell
      PetscInt,Dimension(:),Pointer                :: cellID
      Type(IS)                                     :: cellIS
      PetscInt                                     :: ierr
      PetscInt                                     :: numDof,iDoF
      PetscInt                                     :: iGauss
      PetscInt,Dimension(:),Pointer                :: BCFlag_Loc
      PetscReal,Dimension(:),Pointer               :: RHS_Loc,F_Loc
      PetscReal                                    :: F_Elem
      PetscLogDouble                               :: flops 
      
      Call PetscLogEventBegin(AppCtx%LogInfo%RHSAssemblyBlock_Event,ierr);CHKERRQ(ierr)
      flops = 0.0

      numDof = AppCtx%ElementType(setID)%numDof
      Allocate(F_Loc(numDof))
      Allocate(RHS_Loc(numDof))
      Allocate(BCFlag_Loc(numDof))

      Call DMmeshGetStratumIS(AppCtx%mesh,'Cell Sets',setID,cellIS,ierr); CHKERRQ(ierr)
      Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Do cell = 1,size(cellID)      
         RHS_Loc = 0.0_Kr
         Call SectionRealRestrictClosure(AppCtx%F%Sec,AppCtx%mesh,cellID(cell),numDof,F_Loc,ierr);CHKERRQ(ierr)
         Call SectionIntRestrictClosure(AppCtx%BCFlag%Sec,AppCtx%mesh,cellID(cell),numDof,BCFlag_Loc,ierr);CHKERRQ(ierr)
         Do iGauss = 1,Size(AppCtx%Elem(cellID(cell)+1)%Gauss_C)
            F_Elem = 0.0_Kr
            Do iDoF = 1,numDof
               F_Elem = F_Elem + AppCtx%Elem(cellID(cell)+1)%BF(iDoF,iGauss) * F_Loc(iDoF)
               flops = flops + 2.0
            End Do
            Do iDoF = 1,numDof
               If (BCFlag_Loc(iDoF) == 0) Then
                  RHS_Loc(iDoF) = RHS_Loc(iDoF) + AppCtx%Elem(cellID(cell)+1)%Gauss_C(iGauss) * ( F_Elem * AppCtx%Elem(cellID(cell)+1)%BF(iDoF,iGauss) )
                  flops = flops + 3.0
               End If
            End Do
         End Do
         Call SectionRealUpdateClosure(AppCtx%RHS%Sec,AppCtx%mesh,cellID(cell),RHS_Loc,ADD_VALUES,ierr);CHKERRQ(ierr)
      End Do ! cell
      Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)

      DeAllocate(BCFlag_Loc)
      DeAllocate(RHS_Loc)
      DeAllocate(F_Loc)
      Call PetscLogFlops(flops ,ierr);CHKERRQ(ierr)
      Call PetscLogEventEnd(AppCtx%LogInfo%RHSAssemblyBlock_Event,ierr);CHKERRQ(ierr)
   End Subroutine RHSAssemblyBlock

#undef __FUNCT__
#define __FUNCT__ "ComputeEnergy"
!!! Change interface to make compatible with TAO
   Subroutine ComputeEnergy(AppCtx)
      Type(Poisson_AppCtx_Type)                    :: AppCtx
      
      PetscInt                                     :: ierr
      PetscInt                                     :: NumDoF
      PetscReal,Dimension(:),Pointer               :: F,U
      PetscInt                                     :: set,cell
      PetscInt,Dimension(:),Pointer                :: setID,cellID
      Type(IS)                                     :: setIS,cellIS
      PetscInt                                     :: iDoF,iGauss
      Type(VECT)                                   :: Strain_Elem,Stress_Elem      
      PetscReal                                    :: F_Elem,U_Elem
      PetscReal                                    :: MyElasticEnergy,MyExtForcesWork
      PetscLogDouble                               :: flops

      MyElasticEnergy = 0.0_Kr
      MyExtForcesWork = 0.0_Kr
      flops = 0.0_Kr
      Call DMmeshGetLabelIdIS(AppCtx%mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         numDof = AppCtx%ElementType(setID(set))%numDof
         Call DMmeshGetStratumIS(AppCtx%mesh,'Cell Sets',setID(set),cellIS,ierr); CHKERRQ(ierr)
         Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         Do cell = 1,size(cellID)      
            Allocate(F(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%F%Sec,AppCtx%mesh,cellID(cell),NumDoF,F,ierr);CHKERRQ(ierr)
            Allocate(U(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%U%Sec,AppCtx%mesh,cellID(cell),NumDoF,U,ierr);CHKERRQ(ierr)
            Do iGauss = 1,Size(AppCtx%Elem(cellID(cell)+1)%BF,2)
               Strain_Elem = 0.0_Kr
               Stress_Elem = 0.0_Kr
               F_Elem      = 0.0_Kr
               U_Elem      = 0.0_Kr
               Do iDoF = 1,NumDoF
                  Stress_Elem = Stress_Elem + AppCtx%Elem(cellID(cell)+1)%Grad_BF(iDoF,iGauss) * U(iDoF)
                  Strain_Elem = Strain_Elem + AppCtx%Elem(cellID(cell)+1)%Grad_BF(iDoF,iGauss) * U(iDoF)
                  F_Elem = F_Elem + AppCtx%Elem(cellID(cell)+1)%BF(iDoF,iGauss) * F(iDoF)
                  U_Elem = U_Elem + AppCtx%Elem(cellID(cell)+1)%BF(iDoF,iGauss) * U(iDoF)
                  flops = flops + 4
               End Do
               MyElasticEnergy = MyElasticEnergy + AppCtx%Elem(cellID(cell)+1)%Gauss_C(iGauss) * ( (Strain_Elem .DotP. Strain_Elem) * 0.5_Kr)
               MyExtForcesWork = MyExtForcesWork - AppCtx%Elem(cellID(cell)+1)%Gauss_C(iGauss) * F_Elem * U_Elem
               flops = flops + 5
            End Do
            DeAllocate(U)
            DeAllocate(F)
         End Do ! cell
         Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
      End Do ! set
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      Call MPI_AllReduce(MyElasticEnergy,AppCtx%ElasticEnergy,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
      Call MPI_AllReduce(MyExtForcesWork,AppCtx%ExtForcesWork,1,MPIU_SCALAR,MPI_SUM,PETSC_COMM_WORLD,ierr);CHKERRQ(ierr)
      AppCtx%TotalEnergy = AppCtx%ElasticEnergy + AppCtx%ExtForcesWork
      Call PetscLogFlops(flops,ierr)
   End Subroutine ComputeEnergy
 
#undef __FUNCT__
#define __FUNCT__ "ComputeGradU"
   Subroutine ComputeGradU(AppCtx)
      Type(Poisson_AppCtx_Type)                       :: AppCtx
      
      PetscInt                                     :: ierr
      Type(VECT)                                   :: Grad
      PetscReal                                    :: Vol
      PetscInt                                     :: numDof
      PetscReal,Dimension(:),Pointer               :: U
      PetscInt                                     :: numDim
      PetscInt                                     :: set,cell
      PetscInt,Dimension(:),Pointer                :: setID,cellID
      Type(IS)                                     :: setIS,cellIS
      PetscInt                                     :: iDoF,iGauss
      PetscReal,Dimension(:),Pointer               :: Grad_Ptr
      PetscLogDouble                               :: flops = 0

      Call DMmeshGetDimension(AppCtx%mesh,numDim,ierr);CHKERRQ(ierr)

      Call DMmeshGetCellSectionReal(AppCtx%mesh,'GradU',numDim,AppCtx%GradU,ierr);CHKERRQ(ierr)
      Allocate(Grad_Ptr(numDim))
      
      Call DMmeshGetLabelIdIS(AppCtx%mesh,'Cell Sets',setIS,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Do set = 1,size(setID)
         numDof = AppCtx%ElementType(setID(set))%numDof
         Call DMmeshGetStratumIS(AppCtx%mesh,'Cell Sets',setID(set),cellIS,ierr); CHKERRQ(ierr)
         Call ISGetIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         Do cell = 1,size(cellID)      
            Allocate(U(NumDoF))
            Call SectionRealRestrictClosure(AppCtx%U%Sec,AppCtx%mesh,cellID(cell),NumDoF,U,ierr);CHKERRQ(ierr)
            Grad = 0.0_Kr
            Vol  = 0.0_Kr
            Do iGauss = 1,size(AppCtx%Elem(cellID(cell)+1)%Gauss_C)
               Do iDoF = 1,NumDoF
                  Grad = Grad + (AppCtx%Elem(cellID(cell)+1)%Gauss_C(iGauss) * U(iDoF) ) * AppCtx%Elem(cellID(cell)+1)%Grad_BF(iDoF,iGauss) 
                  Vol = Vol + AppCtx%Elem(cellID(cell)+1)%Gauss_C(iGauss) * AppCtx%Elem(cellID(cell)+1)%BF(iDoF,iGauss)
                  flops = flops + 3
               End Do
            End Do
            Grad = Grad / Vol
#if defined PB_2D
            Grad_Ptr = (/ Grad%X,Grad%Y /)
#elif defined PB_3D
            Grad_Ptr = (/ Grad%X,Grad%Y,Grad%Z /)
#endif
            Call SectionRealUpdateClosure(AppCtx%GradU,AppCtx%mesh,cellID(cell),Grad_Ptr,INSERT_VALUES,ierr)
            DeAllocate(U)
         End Do ! cell
         Call ISRestoreIndicesF90(cellIS,cellID,ierr);CHKERRQ(ierr)
         Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
      End Do ! set
      Call ISRestoreIndicesF90(setIS,setID,ierr);CHKERRQ(ierr)
      Call ISDestroy(setIS,ierr);CHKERRQ(ierr)
      DeAllocate(Grad_Ptr)
      Call PetscLogFlops(flops,ierr)
   End Subroutine ComputeGradU
 
#undef __FUNCT__
#define __FUNCT__ "SimplePoissonFinalize"
!!! Check for missing destroys / deallocate
   Subroutine SimplePoissonFinalize(AppCtx)   
      Type(Poisson_AppCtx_Type)                    :: AppCtx

      PetscInt                                     :: ierr
      Character(len=MEF90_MXSTRLEN)                :: filename
      Type(PetscViewer)                            :: LogViewer

      Call SectionRealDestroy(AppCtx%GradU,ierr);CHKERRQ(ierr)
      Call FieldDestroy(AppCtx%U)
      Call FieldDestroy(AppCtx%F)
      Call FieldDestroy(AppCtx%RHS)

      Call SectionIntDestroy(AppCtx%BCFlag,ierr);CHKERRQ(ierr)
      Call MatDestroy(AppCtx%K,ierr);CHKERRQ(ierr)
      Call KSPDestroy(AppCtx%KSPU,ierr);CHKERRQ(ierr)
      Call DMDestroy(AppCtx%mesh,ierr);CHKERRQ(ierr)
      
      If (AppCtx%AppParam%verbose > 1) Then
         Call PetscViewerFlush(AppCtx%AppParam%MyLogViewer,ierr);CHKERRQ(ierr)
         Call PetscViewerFlush(AppCtx%AppParam%LogViewer,ierr);CHKERRQ(ierr)
         Call PetscViewerDestroy(AppCtx%AppParam%MyLogViewer,ierr);CHKERRQ(ierr)
         Call PetscViewerDestroy(AppCtx%AppParam%LogViewer,ierr);CHKERRQ(ierr)
      End If
      Write(filename,103) Trim(AppCtx%AppParam%prefix)
103 Format(A,'-logsummary.txt')
      Call PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,LogViewer,ierr);CHKERRQ(ierr)
      Call PetscViewerASCIIPrintf(LogViewer,filename,ierr);CHKERRQ(ierr)
      Call PetscLogView(LogViewer,ierr);CHKERRQ(ierr)
      Call PetscViewerDestroy(LogViewer,ierr);CHKERRQ(ierr)
!       Call PetscLogPrintSummary(PETSC_COMM_WORLD,filename,ierr);CHKERRQ(ierr)
      Call MEF90_Finalize()
   End Subroutine SimplePoissonFinalize

End Module M_POISSON
