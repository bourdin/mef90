#include "SimplePoisson.inc"
Module m_PoissonGlobalProperties_Type
#include "finclude/petscdef.h"
   Use m_MEF90
   Use,Intrinsic :: iso_c_binding
   Implicit none
   Private  
   Public :: PoissonGlobalProperties_Type

   Type PoissonGlobalProperties_Type   
      PetscInt       :: verbose
      PetscEnum      :: FileFormat
      PetscEnum      :: FileMode
      PetscEnum      :: LoadingType
      PetscReal      :: TimeMin
      PetscReal      :: TimeMax
      PetscInt       :: numTimeStep
      PetscInt       :: tempoffset
      PetscInt       :: fluxoffset
   End Type PoissonGlobalProperties_Type   
End Module m_PoissonGlobalProperties_Type

Module m_PoissonGlobalProperties   
#include "finclude/petscdef.h"
#include "finclude/petscbagdef.h"
   Use m_MEF90
   Use m_PoissonGlobalProperties_Type
   Implicit None

   Private :: PetscBagGetData

   PetscSizeT,protected    :: sizeofPoissonGlobalProperties

   Enum,bind(c)
      enumerator  :: Poisson_EXOSingle = 0,    &
                     Poisson_EXOSplit 
   End Enum
   Character(len = MEF90_MXSTRLEN),dimension(5),protected   :: Poisson_FileFormatList
   
   Enum,bind(c)
      enumerator  :: Poisson_Replace = 0,    &
                     Poisson_Append 
   End Enum
   Character(len = MEF90_MXSTRLEN),dimension(5),protected   :: Poisson_FileModeList
   
   Enum,bind(c) 
      Enumerator  :: Poisson_CST = 0,        &   
                     Poisson_MIL,            &
                     Poisson_FILE
   End Enum
   Character(len = MEF90_MXSTRLEN),dimension(6),protected   :: Poisson_LoadingList

   Interface PetscBagGetData
      Subroutine PetscBagGetData(bag,data,ierr)
         Use m_PoissonGlobalProperties_Type
         PetscBag                                     :: bag
         Type(PoissonGlobalProperties_Type),pointer   :: data
         PetscErrorCode                               :: ierr
      End subroutine PetscBagGetData
   End interface
Contains

#undef __FUNCT__
#define __FUNCT__ "PoissonGlobalPropertiesInitialize"
!!!
!!!  PoissonGlobalPropertiesInitialize - setup all protected variables related to m_PoissonGlobalProperties:
!!!                                      namely: sizeofPoissonGlobalProperties
!!!
   Subroutine PoissonGlobalPropertiesInitialize(ierr)
      PetscErrorCode,intent(OUT)                   :: ierr

      Type(PoissonGlobalProperties_Type),Target    :: GlobalProperties
      character(len=1),pointer                     :: dummychar(:)
      PetscSizeT                                   :: sizeofchar
   
      Call PetscDataTypeGetSize(PETSC_CHAR,sizeofchar,ierr)
      sizeofPoissonGlobalProperties = size(transfer(GlobalProperties,dummychar))*sizeofchar
      
      Poisson_FileModeList(1) = 'Replace'
      Poisson_FileModeList(2) = 'Append'
      Poisson_FileModeList(3) = 'FileMode'
      Poisson_FileModeList(4) = '_Poisson_FileMode'
      Poisson_FileModeList(5) = ''
      
      Poisson_FileFormatList(1) = 'EXOSingle'
      Poisson_FileFormatList(2) = 'EXOSplit'
      Poisson_FileFormatList(3) = 'FileFormat'
      Poisson_FileFormatList(4) = '_Poisson_FileFormat'
      Poisson_FileFormatList(5) = '' 
      
      Poisson_LoadingList(1) = 'CST'
      Poisson_LoadingList(2) = 'MIL'
      Poisson_LoadingList(3) = 'FILE'
      Poisson_LoadingList(4) = 'Loading'
      Poisson_LoadingList(5) = '_Poisson_Loading'
      Poisson_LoadingList(6) = ''

   End Subroutine PoissonGlobalPropertiesInitialize

#undef __FUNCT__
#define __FUNCT__ "PetscBagGetDataPoissonGlobalProperties"
!!!
!!!  PetscBagGetDataPoissonGlobalProperties - Custom interface to PetscGetData
!!!
   Subroutine PetscBagGetDataPoissonGlobalProperties(bag,data,ierr)
      PetscBag                                        :: bag
      Type(PoissonGlobalProperties_Type),pointer      :: data
      PetscErrorCode                                  :: ierr
      
      Call PetscBagGetData(bag,data,ierr)
   End Subroutine PetscBagGetDataPoissonGlobalProperties

#undef __FUNCT__
#define __FUNCT__ "PetscBagRegisterPoissonGlobalProperties"
!!!
!!!  PetscBagRegisterPoissonGlobalProperties - Register a PetscBag associated with a PoissonGlobalProperties_Type
!!!                                            Typically, the bag is in a MEF90Ctx_Type
!!!
   Subroutine PetscBagRegisterPoissonGlobalProperties(bag,name,prefix,default,ierr)
      PetscBag                                        :: bag
      Character(len=*),intent(IN)                     :: prefix,name
      type(PoissonGlobalProperties_Type),intent(IN)   :: default
      PetscErrorCode,intent(OUT)                      :: ierr

      Type(PoissonGlobalProperties_Type),pointer      :: GlobalProperties
      PetscBool                                       :: flg

      Call PetscBagGetDataPoissonGlobalProperties(bag,GlobalProperties,ierr);CHKERRQ(ierr)
      Call PetscBagSetName(bag,trim(name),"Global properties object",ierr);CHKERRQ(ierr)
      Call PetscBagSetOptionsPrefix(bag,trim(prefix),ierr);CHKERRQ(ierr)
      Call PetscBagRegisterInt(bag,GlobalProperties%verbose,default%verbose,'verbose','verbosity level',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,GlobalProperties%FileFormat,Poisson_FileFormatList,default%FileFormat,'file_format','I/O file format.',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,GlobalProperties%FileMode,Poisson_FileModeList,default%FileMode,'file_mode','I/O file mode.',ierr);CHKERRQ(ierr)
      Call PetscBagRegisterEnum(bag,GlobalProperties%LoadingType,Poisson_LoadingList,default%LoadingType,'Loading','Loading Type',ierr);CHKERRQ(ierr)
      
      Select Case (GlobalProperties%LoadingType)
         Case (Poisson_MIL)
            Call PetscBagRegisterReal(bag,GlobalProperties%timemin,default%timemin,'mil_timemin','Starting time',ierr);CHKERRQ(ierr)
            Call PetscBagRegisterReal(bag,GlobalProperties%timemax,default%timemax,'mil_timemax','End Time',ierr);CHKERRQ(ierr)
            Call PetscBagRegisterInt(bag,GlobalProperties%numtimeStep,default%numTimeStep,'mil_numTimeStep','number of time step in MIL',ierr);CHKERRQ(ierr)
         Case (Poisson_File)
            Call PetscBagRegisterInt(bag,GlobalProperties%tempoffset,default%tempoffset,'file_tempoffset','Offset at which temperature is stored in file',ierr);CHKERRQ(ierr)
            Call PetscBagRegisterInt(bag,GlobalProperties%fluxoffset,default%fluxoffset,'file_fluxoffset','Offset at which flux is stored in file',ierr);CHKERRQ(ierr)
      End Select
      
      If ( (GlobalProperties%FileMode == Poisson_FILE) .AND. (GlobalProperties%FileMode == Poisson_Replace)) Then
         Call PetscPrintf(PETSC_COMM_WORLD,'[WARNING]: input values of temperature will be overwritten with computed ones in EXO file\n',ierr);CHKERRQ(ierr)
      End If
   End Subroutine PetscBagRegisterPoissonGlobalProperties

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxPoissonGlobalPropertiesCreate"
!!!
!!!  
!!!  MEF90CtxPoissonGlobalPropertiesCreate - Call PetscBagRegisterPoissonGlobalProperties for the global properties
!!!                                          of a MEF90Ctx_Type. Mostly for symmetry with CellSet and VertexSet functions
!!!  
!!!  (c) 2012 Blaise Bourdin bourdin@lsu.edu
!!!
   Subroutine MEF90CtxPoissonGlobalPropertiesCreate(MEF90Ctx,defaultGlobalProperties,ierr)
      Type(MEF90Ctx_Type),intent(OUT)                 :: MEF90Ctx
      Type(PoissonGlobalProperties_Type),Intent(IN)   :: defaultGlobalProperties
      PetscErrorCode,Intent(OUT)                      :: ierr
   
      Type(PoissonGlobalProperties_Type),pointer      :: GlobalProperties
      
      Call PetscBagCreate(PETSC_COMM_WORLD,sizeofPoissonGlobalProperties,MEF90Ctx%GlobalPropertiesBag,ierr)
      Call PetscBagRegisterPoissonGlobalProperties(MEF90Ctx%GlobalPropertiesBag,"SimplePoissonNG",PETSC_NULL_CHARACTER,defaultGlobalProperties,ierr);CHKERRQ(ierr)

      Call PetscBagGetDataPoissonGlobalProperties(MEF90Ctx%GlobalPropertiesBag,GlobalProperties,ierr);CHKERRQ(ierr)
      If (GlobalProperties%verbose > 0) Then
         Call PetscBagView(MEF90Ctx%GlobalPropertiesBag,PETSC_VIEWER_STDOUT_WORLD,ierr);CHKERRQ(ierr)
         Call PetscPrintf(PETSC_COMM_WORLD,"\n",ierr);CHKERRQ(ierr)
      End If
   End Subroutine MEF90CtxPoissonGlobalPropertiesCreate
End Module m_PoissonGlobalProperties