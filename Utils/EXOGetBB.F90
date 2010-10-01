Program EXOGetBB

#include "finclude/petscdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
   Use petsc
   Use petscmesh

   Implicit NONE   

   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO, MyEXO
   Type(Vect2D), Dimension(:), Pointer          :: CoordVect2D
   Type(Vect3D), Dimension(:), Pointer          :: CoordVect3D
   
   PetscBool                                    :: HasPrefix
   PetscReal, Dimension(:,:), Pointer           :: CoordArray
   PetscErrorCode                               :: iErr
   Character(len=256)                           :: prefix
   
   Call MEF90_Initialize()
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, MeshTopology%mesh, ierr); CHKERRQ(iErr)

   Call MeshGetCoordinatesF90(MeshTopology%mesh, CoordArray, iErr); CHKERRQ(iErr)
   Write(*, *) 'Size(CoordArray)', Size(CoordArray,1), Size(CoordArray,2)
   Write(*,*) 'X range: ', MinVal(CoordArray(:,1)), MaxVal(CoordArray(:,1))
   Write(*,*) 'Y range: ', MinVal(CoordArray(:,2)), MaxVal(CoordArray(:,2))
   If (Size(CoordArray,2) == 3) Then
      Write(*,*) 'Z range: ', MinVal(CoordArray(:,3)), MaxVal(CoordArray(:,3))
   End If
   Call MeshRestoreCoordinatesF90(MeshTopology%mesh, CoordArray, iErr); CHKERRQ(iErr)

   Call MeshDestroy(MeshTopology%Mesh)
   Call MEF90_Finalize()
End Program EXOGetBB
