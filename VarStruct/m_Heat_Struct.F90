Module m_Heat_Struct
#include "finclude/petscdef.h"

   Use m_MEF90

   Implicit NONE

   Type MatHeat_Type
      PetscInt                                     :: Type_Law
      PetscReal                                    :: Diffusivity
      PetscReal                                    :: Diffusivity2
      PetscReal                                    :: Diffusivity3
   End Type MatHeat_Type
   
   Type HeatSchemeParam_Type
      PetscInt                                     :: HeatMaxIter
      Character(len=MEF90_MXSTRLEN)                :: HeatTSType
      Character(len=MEF90_MXSTRLEN)                :: TSTypeOpt1
   End Type HeatSchemeParam_Type

   Type Heat_TestCase_Type
      PetscInt                                  :: Index
      Character(len=MEF90_MXSTRLEN)             :: Description
   End Type

 Contains
   Subroutine HeatSchemeParam_View(dSchemeParam, viewer)
      Type(HeatSchemeParam_Type)                :: dSchemeParam
      Type(PetscViewer)                            :: viewer
      PetscInt                                     :: iErr
      Character(len=MEF90_MXSTRLEN)                :: IOBuffer
  
      Write(IOBuffer, "('-HeatMaxIter ', I1, A)")          dSchemeParam%HeatMaxIter, '\n'
      Call PetscViewerASCIIPrintf(viewer, IOBuffer, iErr); CHKERRQ(iErr)

   End Subroutine HeatSchemeParam_View
      
   Subroutine HeatSchemeParam_GetFromOptions(dSchemeParam)
      Type(HeatSchemeParam_Type)                   :: dSchemeParam
      PetscInt                                     :: iErr
      PetscBool                                    :: flag

      dSchemeParam%HeatMaxIter      =  1000
      dSchemeParam%HeatTSType       = 'rosw' 
      dSchemeParam%TSTypeOpt1       = 'ra3pw' 
      
      Call PetscOptionsGetInt(PETSC_NULL_CHARACTER,   '-heatmaxiter',          dSchemeParam%HeatMaxIter, flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-heattstype',          dSchemeParam%HeatTSType,   flag, iErr); CHKERRQ(iErr) 
      Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-tstypeopt1',          dSchemeParam%TSTypeOpt1,   flag, iErr); CHKERRQ(iErr) 

   End Subroutine HeatSchemeParam_GetFromOptions
   
   Subroutine MatHeat_Write(MeshTopology, MatProp, filename)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(MatHeat_Type), Dimension(:), Pointer    :: MatProp
      Character(len=*)                             :: filename
      PetscMPIInt                                  :: rank
      PetscInt                                     :: iBlk
      
      Open(File = filename, Unit = F_OUT, Status = 'Unknown')
      Rewind(F_OUT)
      Write(F_OUT, *) MeshTopology%Num_Elem_Blks_Global
      Do iBlk = 1, Size(MatProp)
      Write(F_OUT,120) iBlk, MatProp(iBlk)%Type_Law, MatProp(iBlk)%Diffusivity, MatProp(iBlk)%Diffusivity2
      End Do
      Close(F_OUT)
      
120   Format(I6, '      ', I6, 2(ES12.5,' '))   
   End Subroutine MatHeat_Write
   
   
   Subroutine MatHeat_Read(MeshTopology, MatProp, filename)
      Type(MeshTopology_Type)                      :: MeshTopology
      Type(MatHeat_Type), Dimension(:), Pointer    :: MatProp
      Character(len=*)                             :: filename

      PetscInt                                     :: iBlk, iErr
      
      PetscInt                                     :: NumBlks, IdxMin, IdxMax
      PetscInt                                     :: Idx, Type_Law 
      PetscReal                                    :: Diffusivity, Diffusivity2
   
      Open(File = filename, Unit = F_IN, Status = 'Unknown', Action = 'Read')
      Rewind(F_IN)
      Read(F_IN, *) NumBlks
      If (NumBlks /= MeshTopology%Num_Elem_Blks_Global) Then
         SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ, 'MatHeatRead: non matching blocks numbers', iErr)
      End If
      !!! Reading the file once first to get the right number of blocks
      IdxMin =  100000000
      IdxMax = -100000000
      Do iBlk = 1, NumBlks
         Read(F_IN, *) Idx
         IdxMin = Min(IdxMin, Idx)
         IdxMax = Max(IdxMax, Idx)
      End Do
      Allocate(MatProp(IdxMin:IdxMax))
      Rewind(F_IN)
      Read(F_IN, *) Idx
      Do iBlk = 1, NumBlks
      Read(F_IN, *) Idx, Type_Law,  Diffusivity, Diffusivity2
         MatProp(Idx)%Type_Law  = Type_Law
         MatProp(Idx)%Diffusivity = Diffusivity
         MatProp(Idx)%Diffusivity2 = Diffusivity2
      End Do
      Close(F_IN)
      Return
!120   Format(I6, '      ', 10(ES12.5,' '))   
!120   Format(*)
   End Subroutine MatHeat_Read

#undef __FUNCT__
#define __FUNCT__ "View_Available_Diffusion_Laws"
   Subroutine View_Available_Diffusion_Laws
      Type(Heat_TestCase_Type), Dimension(:) , Pointer                 :: DiffusionLaw
      Character(len=MEF90_MXSTRLEN)                               ::  IOBuffer
      PetscInt                                                    :: NumLaws   
      PetscInt                                                    :: i, iErr 

      NumLaws = 5
      Allocate(DiffusionLaw(NumLaws))
      Do i = 1, NumLaws 
             DiffusionLaw(i)%Index = i
      End Do
      DiffusionLaw(1)%Description = "Constant diffusion parameter D(T, v) = k"
      DiffusionLaw(2)%Description = "Increasing diffusion paramter, Mensi's law D(T, v) = k exp(B*T)"
      DiffusionLaw(3)%Description = "Decreasing diffusion parameter D(T, v) = k(A-T)**B"
      DiffusionLaw(4)%Description = "Diffusion depending on damage field D(T, v) = k/v"
      DiffusionLaw(5)%Description = "Diffusion depending on crack opening [NOT IMPLEMENTED]"

      Write(IOBuffer, *) '\nDiffusion laws:\n'
      Do i = 1, NumLaws
         Write(IOBuffer, "('   [',I2.2,'] ',A)"), DiffusionLaw(i)%Index, Trim(DiffusionLaw(i)%Description)//'\n'
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr);   CHKERRQ(iErr)
      End Do
   DeAllocate(DiffusionLaw)

   End Subroutine View_Available_Diffusion_Laws

End Module m_Heat_Struct
