Program TestDuplicate

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc
 
   Implicit NONE   

   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO, MyEXO
   Type(Element2D_Scal), Dimension(:), Pointer  :: Elem2DA
   
   PetscBool                                    :: HasPrefix
   PetscBool                                    :: verbose
   PetscErrorCode                               :: iErr
   Character(len=256)                           :: filename
   Character(len=256)                           :: prefix
   Type(SectionReal)                            :: U_Sec, V_Sec
   Type(VecScatter)                             :: scatter
   PetscInt                                     :: i
   PetscInt                                     :: dof
   PetscInt, Dimension(:), Pointer              :: iRows
   PetscReal                                    :: one = 1.0_Kr
   Character(len=256)                           :: IO_Buffer
   Type(PetscViewer)                            :: MeshViewer

   Integer                                      :: cpu_ws = 0
   Integer                                      :: io_ws = 0
   Type(DM)                                     :: tmpDM
   Integer                                      :: exoid
   PetscReal                                    :: vers
     
   Call MEF90_Initialize()
   dof=1
   Call PetscOptionsGetInt(PETSC_NULL_CHARACTER, '-dof', dof, HasPrefix, iErr)    
   Call PetscOptionsHasName(PETSC_NULL_CHARACTER, '-verbose', verbose, iErr)    
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   If (MEF90_Myrank == 0) Then
      exoid = EXOPEN(EXO%filename,EXREAD,cpu_ws,io_ws,vers,ierr)
   End If
   If (MEF90_numprocs == 1) Then
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoid,MeshTopology%mesh,ierr);CHKERRQ(ierr)
   Else
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,exoid,tmpDM,ierr);CHKERRQ(ierr)
      Call DMMeshDistribute(tmpDM,PETSC_NULL_CHARACTER,MeshTopology%mesh,ierr);CHKERRQ(iErr)
      Call DMDestroy(tmpDM,ierr);CHKERRQ(iErr)
   End If
   If (MEF90_Myrank == 0) Then
      Call EXCLOS(exoid,ierr)
   End If

   Call MeshTopologyGetInfo(MeshTopology)
   
   Call DMMeshSetMaxDof(MeshTopology%Mesh, dof, iErr); CHKERRQ(iErr) 
   Call DMMeshGetVertexSectionReal(MeshTopology%mesh, 'U', dof, U_Sec, iErr); CHKERRQ(iErr)
   Call DMMeshCreateGlobalScatter(MeshTopology%mesh, U_Sec, scatter, iErr); CHKERRQ(iErr)



   Call PetscPrintf(PETSC_COMM_WORLD, "\n\nSection U_Sec: \n", iErr); CHKERRQ(iErr)
   Call SectionRealView(U_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   Call SectionRealGetSize(U_Sec, i, iErr)
   Write(IO_Buffer, *) 'Size of U_Sec: ', i, '\n'
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, IO_Buffer, iErr); CHKERRQ(iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr); CHKERRQ(iErr)   
   Call SectionRealDuplicate(U_Sec, V_Sec, iErr); CHKERRQ(iErr)
   Call PetscPrintf(PETSC_COMM_WORLD, "\n\nSection V_Sec: \n", iErr); CHKERRQ(iErr)
   Call SectionRealView(V_Sec, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)

   Call MeshTopologyDestroy(MeshTopology)
   Call SectionRealDestroy(U_Sec, iErr); CHKERRQ(iErr)
   Call SectionRealDestroy(V_Sec, iErr); CHKERRQ(iErr)
   Call VecScatterDestroy(scatter, iErr); CHKERRQ(iErr)
   Call MEF90_Finalize()
End Program TestDuplicate
