Program TestCoordinates

#include "finclude/petscdef.h"

   Use m_MEF90
   Use petsc

   Implicit NONE   

   Type(MeshTopology_Type)                      :: MeshTopology
   Type(EXO_Type)                               :: EXO,MyEXO
   Type(Element2D_Scal),Dimension(:),Pointer    :: Elem2DA
   Type(Vect2D),Dimension(:),Pointer            :: CoordVect2D
   Type(DM)                                     :: tmpDM
   
   PetscBool                                    :: HasPrefix
   PetscReal,Dimension(:,:),Pointer             :: array
   PetscErrorCode                               :: iErr
   Character(len=256)                           :: prefix
   Type(SectionReal)                            :: coordSection
   Integer                                      :: cpu_ws = 0
   Integer                                      :: io_ws = 0
   PetscReal                                    :: vers
   
   Call MEF90_Initialize()

   Call PetscOptionsGetString(PETSC_NULL_CHARACTER,'-p',prefix,HasPrefix,iErr)    
   If (.NOT. HasPrefix) Then
      Call PetscPrintf(PETSC_COMM_WORLD,"No input file prefix given\n",iErr)
      Call MEF90_Finalize()
      STOP
   End If

   EXO%Comm = PETSC_COMM_WORLD
   EXO%filename = Trim(prefix)//'.gen'
   If (MEF90_MyRank == 0) Then
      EXO%exoid = EXOPEN(EXO%filename,EXREAD,cpu_ws,io_ws,vers,ierr)
   End If
   If (MEF90_numprocs == 1) Then
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,EXO%exoid,MeshTopology%mesh,ierr);CHKERRQ(ierr)
   Else
      Call DMMeshCreateExodusNG(PETSC_COMM_WORLD,EXO%exoid,tmpDM,ierr);CHKERRQ(ierr)
      Call DMMeshDistribute(tmpDM,PETSC_NULL_CHARACTER,MeshTopology%mesh,ierr);CHKERRQ(iErr)
      Call DMDestroy(tmpDM,ierr);CHKERRQ(iErr)
   End If
   If (MEF90_MyRank == 0) Then
      Call EXCLOS(EXO%exoid,ierr)
   End If


   Call MeshTopologyGetInfo(MeshTopology)
   Write(*,*) 'OK MeshTopologyGetInfo'

   Call DMMeshGetCoordinatesF90(MeshTopology%mesh,array,iErr); CHKERRQ(iErr)
   Write(MEF90_MyRank + 100,*) 'Size(array)',Size(array,1),Size(Array,2)
   Write(MEF90_MyRank + 100,*) 'array(:,1)'
   Write(MEF90_MyRank + 100,*) array(:,1)
   Write(MEF90_MyRank + 100,*) 'array(:,2)'
   Write(MEF90_MyRank + 100,*) array(:,2)
   Call DMMeshRestoreCoordinatesF90(MeshTopology%mesh,array,iErr); CHKERRQ(iErr)

   Write(MEF90_MyRank + 100,*) 'Now creating the coordinates section'
   Call DMMeshGetSectionReal(MeshTopology%mesh,'coordinates',coordSection,iErr); CHKERRQ(ierr)
   Call SectionRealView(coordSection,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
   Call SectionRealDestroy(coordSection,ierr);CHKERRQ(ierr)

   Call DMMeshGetCoordinatesF90(MeshTopology%mesh,array,iErr); CHKERRQ(iErr)
   Write(MEF90_MyRank + 100,*) 'Size(array)',Size(array,1),Size(Array,2)
   Write(MEF90_MyRank + 100,*) 'array(:,1)'
   Write(MEF90_MyRank + 100,*) array(:,1)
   Write(MEF90_MyRank + 100,*) 'array(:,2)'
   Write(MEF90_MyRank + 100,*) array(:,2)
   Call DMMeshRestoreCoordinatesF90(MeshTopology%mesh,array,iErr); CHKERRQ(iErr)

   Call MEF90_Finalize()
End Program TestCoordinates
