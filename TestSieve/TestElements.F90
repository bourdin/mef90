Program TestElements
#include 'finclude/petscdef.h'
   Use m_MEF90
   Use petsc
   Implicit NONE   

   Type(DM)                                     :: mesh
   Type(EXO_Type)                               :: EXO,MyEXO
   Type(Element2D_Scal),Dimension(:),Pointer    :: Elem2D_Scal
   Type(Element3D_Scal),Dimension(:),Pointer    :: Elem3D_Scal
   PetscInt,Dimension(:),Pointer                :: ElemTypeID
   Type(Element_Type),Dimension(:),Pointer      :: ElemType
   
   PetscBool                                    :: HasPrefix = PETSC_FALSE
   PetscBool                                    :: verbose = PETSC_FALSE
   PetscBool                                    :: splitIO = PETSC_FALSE
   PetscBool                                    :: flg
   PetscErrorCode                               :: iErr
   Character(len = 256)                         :: filename
   Character(len = 256)                         :: prefix
   Character(len = MEF90_MXSTRLEN)              :: IOBuffer
   Type(Field)                                  :: U,V
   Type(SectionReal)                            :: coordSection,GradU
   
   PetscReal,Dimension(:),Pointer               :: valU,valV
   
   Type(IS)                                     :: setIS,cellIS,CellSetGlobalIS
   PetscInt,Dimension(:),Pointer                :: setID,cellID
   PetscInt,Dimension(:),Pointer                :: sizeScal
   PetscInt                                     :: set,cell,vertex
   PetscInt                                     :: numDim,numCell,numVertex,numSizes
   PetscReal,Dimension(:),Pointer               :: coord
   
   PetscInt                                     :: conesize
   Type(DM)                                     :: tmpDM
   Integer                                      :: cpu_ws = 8
   Integer                                      :: io_ws = 8
   PetscReal                                    :: vers
   Integer                                      :: exoidIN
   Integer                                      :: offset = 1
   PetscReal                                    :: LpNorm
   PetscInt                                     :: GaussOrder = 2
     
   Call MEF90_Initialize()
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER,'-verbose',verbose,iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-p',prefix,HasPrefix,iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD,'No input file prefix given\n',iErr)
      Call MEF90_Finalize()
      STOP
   End If
   Call PetscOptionsGetBool(PETSC_NULL_CHARACTER,'-splitIO',splitIO,flg,ierr);CHKERRQ(ierr)
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-gaussorder',GaussOrder,HasPrefix,iErr)    

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   If (MEF90_Myrank == 0) Then
      exoidIN = EXOPEN(EXO%filename,EXREAD,cpu_ws,io_ws,vers,ierr)
   End If
   If (MEF90_numprocs == 1) Then
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoidIN,mesh,ierr);CHKERRQ(ierr)
   Else
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoidIN,tmpDM,ierr);CHKERRQ(ierr)
      Call DMMeshDistribute(tmpDM,PETSC_NULL_CHARACTER,mesh,ierr);CHKERRQ(iErr)
      Call DMDestroy(tmpDM,ierr);CHKERRQ(iErr)
   End If
   Call DMmeshGetDimension(mesh,numDim,ierr);CHKERRQ(ierr)
   
   If (splitIO) Then
      EXO%Comm = PETSC_COMM_SELF
      Write(EXO%filename,105) trim(prefix),MEF90_MyRank
      cpu_ws = 8
      io_ws = 8
      EXO%exoid = EXCRE(trim(EXO%filename),EXCLOB,cpu_ws,io_ws,ierr)
      Call DMmeshViewExodusSplit(mesh,EXO%exoid,ierr)
   Else
      EXO%Comm = PETSC_COMM_WORLD
      Write(EXO%filename,106) trim(prefix)
      cpu_ws = 8
      io_ws = 8
      EXO%exoid = EXCRE(EXO%filename,EXCLOB,cpu_ws,io_ws,ierr)
      If (MEF90_Myrank == 0) Then
         Call EXCOPY(exoidIN,EXO%exoid,ierr)
         Call EXOFormat_Scal(EXO,numDim)
      End If
   End If
   If (MEF90_Myrank == 0) Then
      Call EXCLOS(exoidIN,ierr)
   End If
   Call EXOFormat_Scal(EXO,numDim)
105 Format(A,'-',I4.4,'.gen')
106 Format(A,'_out.gen')

   Call DMMeshSetMaxDof(mesh,numDim,iErr); CHKERRQ(iErr) 
   !!!
   !!! Compute global IS for cell and vertex sets
   !!!
   Call DMmeshGetLabelIdIS(mesh,'Cell Sets',CellSetGlobalIS,ierr);CHKERRQ(ierr)
   Call MEF90_ISAllGatherMerge(PETSC_COMM_WORLD,CellSetGlobalIS)
      
   !!!
   !!! Allocate Element array and set element types
   !!!
   Call DMmeshGetStratumSize(mesh,'height',0,numCell,ierr);CHKERRQ(ierr)

   Call ISGetIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   Allocate(ElemTypeID(size(setID)))
   Allocate(ElemType(size(setID)))

   If (numDim == 2) Then
      Allocate(Elem2D_Scal(numCell))
      ElemTypeID = MEF90_P1_Lagrange_2D_Scal%shortID
   Else
      Allocate(Elem3D_Scal(numCell))
      ElemTypeID = MEF90_P1_Lagrange_3D_Scal%shortID
   End If

   numSizes = size(setID)
   Call PetscOptionsGetIntArray(PETSC_NULL_CHARACTER,'-elem_type',ElemTypeID,numSizes,flg,ierr);CHKERRQ(ierr)
   If ((numSizes /= 0) .AND. (numSizes /= size(setID))) Then
      Write(*,*) '[ERROR], was expecting ',size(setID), 'element types, but got ', numSizes
      Call MEF90_Finalize()
      STOP
   End If
   Do set = 1,size(setID)
      Call Element_TypeFindByID(elemTypeID(set),ElemType(set))
      Write(*,*) 'Set: ', set, 'element type', trim(ElemType(set)%name)
      Call DMmeshGetStratumIS(mesh,'Cell Sets',setID(set),cellIS,ierr); CHKERRQ(ierr)
      If (numDim == 2) Then
         Call ElementInit(mesh,cellIS,Elem2D_Scal,GaussOrder,ElemType(setID(set)))
      Else
         Call ElementInit(mesh,cellIS,Elem3D_Scal,GaussOrder,ElemType(setID(set)))
      End If
      Call ISDestroy(cellIS,ierr);CHKERRQ(ierr)
   End Do    
   DeAllocate(ElemTypeID)
   
   !!! Allocate the Section for U and GradU
   Allocate(SizeScal(1))
   SizeScal = 1
   Call FieldCreateVertex(U,'U',mesh,SizeScal)
   Call FieldCreateVertex(V,'V',mesh,SizeScal)
   DeAllocate(sizeScal)
   Call DMMeshGetVertexSectionReal(mesh,'GradU',numDim,gradU,ierr);CHKERRQ(ierr)
   Allocate(ValU(1))
   Allocate(ValV(1))

   !!! Initialize the section for U   
   Call DMMeshGetSectionReal(mesh,'coordinates',coordSection,ierr); CHKERRQ(ierr)
   Call DMmeshGetStratumSize(mesh,'depth',0,numVertex,ierr);CHKERRQ(ierr)

   !!! U = 1, V=1
   Do vertex = 1, numVertex
      Call SectionRealRestrict(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
      valU = 1
      ValV = 1
      Call SectionRealUpdate(U%Sec,vertex+numCell-1,valU,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealUpdate(V%Sec,vertex+numCell-1,valV,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealRestore(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
   End Do
   Do set = 1,size(setID)
      If (numDim == 2) Then
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem2D_Scal,LpNorm,ierr);CHKERRQ(ierr)
      Else
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem3D_Scal,LpNorm,ierr);CHKERRQ(ierr)
      End If
         Write(*,*) 'Set ', setID(set), '<1,1>_2 = ', LpNorm
   End Do
   
   Do vertex = 1, numVertex
      Call SectionRealRestrict(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
      valU = 1
      ValV = coord(1)
      Call SectionRealUpdate(U%Sec,vertex+numCell-1,valU,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealUpdate(V%Sec,vertex+numCell-1,valV,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealRestore(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
   End Do
   Do set = 1,size(setID)
      If (numDim == 2) Then
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem2D_Scal,LpNorm,ierr);CHKERRQ(ierr)
      Else
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem3D_Scal,LpNorm,ierr);CHKERRQ(ierr)
      End If
         Write(*,*) 'Set ', setID(set), '<1,X>_2 = ', LpNorm
   End Do
   
   Do vertex = 1, numVertex
      Call SectionRealRestrict(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
      valU = 1
      ValV = coord(2)
      Call SectionRealUpdate(U%Sec,vertex+numCell-1,valU,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealUpdate(V%Sec,vertex+numCell-1,valV,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealRestore(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
   End Do
   Do set = 1,size(setID)
      If (numDim == 2) Then
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem2D_Scal,LpNorm,ierr);CHKERRQ(ierr)
      Else
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem3D_Scal,LpNorm,ierr);CHKERRQ(ierr)
      End If
         Write(*,*) 'Set ', setID(set), '<1,Y>_2 = ', LpNorm
   End Do
   
   If (numDim == 3) Then
      Do vertex = 1, numVertex
         Call SectionRealRestrict(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
         valU = 1
         ValV = coord(3)
         Call SectionRealUpdate(U%Sec,vertex+numCell-1,valU,INSERT_VALUES,iErr); CHKERRQ(iErr)
         Call SectionRealUpdate(V%Sec,vertex+numCell-1,valV,INSERT_VALUES,iErr); CHKERRQ(iErr)
         Call SectionRealRestore(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
      End Do
      Do set = 1,size(setID)
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem3D_Scal,LpNorm,ierr);CHKERRQ(ierr)
         Write(*,*) 'Set ', setID(set), '<1,Z>_2 = ', LpNorm
      End Do
   End If
      
   Do vertex = 1, numVertex
      Call SectionRealRestrict(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
      valU = coord(1)
      ValV = coord(2)
      Call SectionRealUpdate(U%Sec,vertex+numCell-1,valU,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealUpdate(V%Sec,vertex+numCell-1,valV,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealRestore(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
   End Do
   If (numDim == 2) Then
      Do set = 1,size(setID)
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,U%Sec,Elem2D_Scal,LpNorm,ierr);CHKERRQ(ierr)
         Write(*,*) 'Set ', setID(set), '<X,X>_2 = ', LpNorm
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem2D_Scal,LpNorm,ierr);CHKERRQ(ierr)
         Write(*,*) 'Set ', setID(set), '<X,Y>_2 = ', LpNorm
         Call SectionRealL2DotProduct(mesh,setID(set),V%Sec,V%Sec,Elem2D_Scal,LpNorm,ierr);CHKERRQ(ierr)
         Write(*,*) 'Set ', setID(set), '<Y,Y>_2 = ', LpNorm
      End Do
   Else
      Do set = 1,size(setID)
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,U%Sec,Elem3D_Scal,LpNorm,ierr);CHKERRQ(ierr)
         Write(*,*) 'Set ', setID(set), '<X,X>_2 = ', LpNorm
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem3D_Scal,LpNorm,ierr);CHKERRQ(ierr)
         Write(*,*) 'Set ', setID(set), '<X,Y>_2 = ', LpNorm
         Call SectionRealL2DotProduct(mesh,setID(set),V%Sec,V%Sec,Elem3D_Scal,LpNorm,ierr);CHKERRQ(ierr)
         Write(*,*) 'Set ', setID(set), '<Y,Y>_2 = ', LpNorm
      End Do
   
   End If

   Do vertex = 1, numVertex
      Call SectionRealRestrict(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
      valU = 1
      ValV = coord(1)**2
      Call SectionRealUpdate(U%Sec,vertex+numCell-1,valU,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealUpdate(V%Sec,vertex+numCell-1,valV,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealRestore(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
   End Do
   If (numDim == 2) Then
      Do set = 1,size(setID)
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem2D_Scal,LpNorm,ierr);CHKERRQ(ierr)
         Write(*,*) 'Set ', setID(set), '<1,X^2>_2 = ', LpNorm
      End Do
   Else
      Do set = 1,size(setID)
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem3D_Scal,LpNorm,ierr);CHKERRQ(ierr)
         Write(*,*) 'Set ', setID(set), '<1,X^2>_2 = ', LpNorm
      End Do   
   End If

   Do vertex = 1, numVertex
      Call SectionRealRestrict(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
      valU = 1
      ValV = coord(1)*coord(2)
      Call SectionRealUpdate(U%Sec,vertex+numCell-1,valU,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealUpdate(V%Sec,vertex+numCell-1,valV,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealRestore(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
   End Do
   If (numDim == 2) Then
      Do set = 1,size(setID)
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem2D_Scal,LpNorm,ierr);CHKERRQ(ierr)
         Write(*,*) 'Set ', setID(set), '<1,XY>_2 = ', LpNorm
      End Do
   End If
   
   Do vertex = 1, numVertex
      Call SectionRealRestrict(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
      valU = 1
      ValV = coord(2)**2
      Call SectionRealUpdate(U%Sec,vertex+numCell-1,valU,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealUpdate(V%Sec,vertex+numCell-1,valV,INSERT_VALUES,iErr); CHKERRQ(iErr)
      Call SectionRealRestore(coordSection,vertex+numCell-1,coord,iErr); CHKERRQ(iErr)
   End Do
   If (numDim == 2) Then
      Do set = 1,size(setID)
         Call SectionRealL2DotProduct(mesh,setID(set),U%Sec,V%Sec,Elem2D_Scal,LpNorm,ierr);CHKERRQ(ierr)
         Write(*,*) 'Set ', setID(set), '<1,Y^2>_2 = ', LpNorm
      End Do
   End If


   Call ISRestoreIndicesF90(CellSetGlobalIS,setID,ierr);CHKERRQ(ierr)
   DeAllocate(valU)
   DeAllocate(valV)
   DeAllocate(ElemType)
   Call ISDestroy(CellSetGlobalIS)
   Call FieldDestroy(U)
   Call FieldDestroy(V)
   Call SectionRealDestroy(GradU,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(coordSection,ierr);CHKERRQ(ierr)
   If (EXO%exoid /= 0) Then
      Call EXCLOS(EXO%exoid,ierr)
   End If
   Call MEF90_Finalize()

Contains
   Subroutine EXOFormat_Scal(EXO,numDim)
      Type(EXO_Type),intent(INOUT)                 :: EXO
      PetscInt,intent(IN)                          :: numDim
      PetscErrorCode                               :: ierr
   
      EXO%num_nsproperties = 0
      EXO%num_ebproperties = 0
      EXO%num_globvariables = 3
      EXO%num_vertvariables = 1
      EXO%num_cellvariables = numDim
      Call EXPVP (EXO%exoid,'g',EXO%num_globvariables,ierr)
      Call EXPVAN(EXO%exoid,'g',EXO%num_globvariables,(/'L1 norm','L2 norm','H1 semi-norm'/),ierr)
      Call EXPVP (EXO%exoid,'n',EXO%num_vertvariables,ierr)
      Call EXPVAN(EXO%exoid,'n',EXO%num_vertvariables,(/'U'/),ierr)

      Call EXPVP (EXO%exoid,'e',EXO%num_cellvariables,ierr)
      print*, numDim
      If (numDim == 2) Then
         Call EXPVAN(EXO%exoid,'e',EXO%num_cellvariables,(/'Grad U_X','Grad U_Y'/),ierr)
      Else
         Call EXPVAN(EXO%exoid,'e',EXO%num_cellvariables,(/'Grad U_X','Grad U_Y','Grad U_Z'/),ierr)
      End If
   End Subroutine EXOFormat_Scal
End Program TestElements
