Program TestCoordinates

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscdmmeshdef.h"

   Use m_MEF90
   Use petsc
   Use petscvec
   Use petscmat
   Use petscdmmesh

   Implicit NONE   

   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO, MyEXO
   Type(Element2D_Scal), Dimension(:), Pointer  :: Elem2DA
   Type(Vect2D), DImension(:), Pointer          :: CoordVect2D
   
   PetscBool                                   :: HasPrefix
   PetscReal, Dimension(:,:), Pointer           :: array
   PetscErrorCode                               :: iErr
   PetscInt                                     :: NumVert
   Character(len=256)                           :: prefix
   Type(Mesh)                                   :: Tmp_Mesh
   Type(SectionReal)                            :: coordSection
   
   Call MEF90_Initialize()
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-p', prefix, HasPrefix, iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file prefix given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'


   Call MeshCreateExodus(PETSC_COMM_WORLD, EXO%filename, Tmp_mesh, ierr); CHKERRQ(iErr)
   Call MeshDistribute(Tmp_mesh, PETSC_NULL_CHARACTER, MeshTopology%mesh, ierr); CHKERRQ(iErr)
   Call MeshDestroy(Tmp_mesh, ierr); CHKERRQ(iErr)
   Call MeshTopologyGetInfo(MeshTopology, PETSC_COMM_WORLD)
   Write(*,*) 'OK MeshTopologyGetInfo'

   Call MeshGetCoordinatesF90(MeshTopology%mesh, array, iErr); CHKERRQ(iErr)
   Write(MEF90_MyRank + 100, *) 'Size(array)', Size(array,1), Size(Array,2)
   Write(MEF90_MyRank + 100, *) 'array(:,1)'
   Write(MEF90_MyRank + 100, *) array(:,1)
   Write(MEF90_MyRank + 100, *) 'array(:,2)'
   Write(MEF90_MyRank + 100, *) array(:,2)
   Call MeshRestoreCoordinatesF90(MeshTopology%mesh, array, iErr); CHKERRQ(iErr)

   Write(MEF90_MyRank + 100, *) 'Now creating the coordinates section'
   Call MeshGetSectionReal(MeshTopology%mesh, 'coordinates', coordSection, iErr); CHKERRQ(ierr)

   Call MeshGetCoordinatesF90(MeshTopology%mesh, array, iErr); CHKERRQ(iErr)
   Write(MEF90_MyRank + 100, *) 'Size(array)', Size(array,1), Size(Array,2)
   Write(MEF90_MyRank + 100, *) 'array(:,1)'
   Write(MEF90_MyRank + 100, *) array(:,1)
   Write(MEF90_MyRank + 100, *) 'array(:,2)'
   Write(MEF90_MyRank + 100, *) array(:,2)
   Call MeshRestoreCoordinatesF90(MeshTopology%mesh, array, iErr); CHKERRQ(iErr)


   Call MEF90_Finalize()
End Program TestCoordinates
