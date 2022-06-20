Program TestMaterials
#include <petsc/finclude/petsc.h>
   Use m_MEF90
   Use m_MEF90_Materials
   Use petsc
   Implicit NONE   

   PetscInt                            :: i,numMat=1
   PetscBool                           :: flg
   PetscBag,Dimension(:),Pointer       :: matBag2D,matBag3D
   Type(MEF90MatProp2D_Type),pointer   :: matProp2D
   Type(MEF90MatProp3D_Type),pointer   :: matProp3D
   Type(MatS2D)                        :: E2D
   Type(MatS3D)                        :: E3D
   Character(len=256)                  :: IOBuffer
   Character(len=80)                   :: name,prefix
   PetscErrorCode                      :: ierr
   

   PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
   PetscCallA(MEF90Initialize(ierr))
   PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,'-n',numMat,flg,ierr))
   
   Allocate(matBag2D(numMat))
   Allocate(matBag3D(numMat))

   Do i = 1,numMat   
      Write(iobuffer,99) i
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))

      PetscCallA(PetscBagCreate(PETSC_COMM_WORLD,sizeofMEF90MatProp2D,matBag2D(i),ierr))
      write(name,200) i
      write(prefix,201) i
      PetscCallA(PetscBagRegisterMEF90MatProp(matBag2D(i),name,prefix,MEF90Mathium2D,ierr))
      PetscCallA(PetscBagView(matBag2D(i),PETSC_VIEWER_STDOUT_WORLD,ierr))
      PetscCallA(PetscBagGetDataMEF90MatProp(matBag2D(i),matProp2D,ierr))
      Write(*,*) 'MatProp2D: ',matProp2D
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD,'\n',ierr))
      E2D = MEF90MatS2DIdentity
      Write(*,*) matprop2D%HookesLaw*E2D
   EndDo


   Do i = 1,numMat   
      Write(iobuffer,99) i
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))

      PetscCallA(PetscBagCreate(PETSC_COMM_WORLD,sizeofMEF90MatProp3D,matBag3D(i),ierr))
      write(name,300) i
      write(prefix,301) i

      PetscCallA(PetscBagRegisterMEF90MatProp(matBag3D(i),name,prefix,MEF90Mathium3D,ierr))
      PetscCallA(PetscBagView(matBag3D(i),PETSC_VIEWER_STDOUT_WORLD,ierr))
      PetscCallA(PetscBagGetDataMEF90MatProp(matBag3D(i),matProp3D,ierr))
      Write(*,*) 'MatProp3D: ',matProp3D

      matProp3D%Density = -123456.0
      PetscCallA(PetscBagSetFromOptions(matBag3D(i),ierr))
      ! PetscBagSetFromOptions resets only to CL options, not to default options
      ! This could probably be fixed by inserting the options in the registration routine, if desired
      ! Calling PetscBagSetFromOptions means the the help message will be re-displayed
      PetscCallA(PetscBagGetDataMEF90MatProp(matBag3D(i),matProp3D,ierr))
      PetscCallA(PetscPrintf(PETSC_COMM_WORLD,'\n',ierr))
   EndDo

99  format('registering material ',I2.2,'\n')
200 format('2D material ',I2.2)
201 format('mat2D',I2.2,'_')
300 format('3D material ',I2.2)
301 format('mat3D',I2.2,'_')
   DeAllocate(matBag2D)
   DeAllocate(matBag3D)
   PetscCallA(MEF90Finalize(ierr))
   PetscCallA(PetscFinalize(ierr))
End Program TestMaterials
