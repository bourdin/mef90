#include "../MEF90/mef90.inc"
module m_MEF90_DefMechAT_class
#include "petsc/finclude/petsc.h"

   use petscsys
   use m_MEF90_Parameters
   use m_MEF90_LinAlg
   use m_MEF90_BaseClass
   use iso_c_binding
   implicit none(type)

   private
   public :: MEF90DefMechAT_Type
   public :: MEF90DefMechAT_setFromOptions

!!!
!!!
!!!  MEF90_DefMechAT_Type: The abstract class used to define a generalized Ambrosio-Tortorelli phase field model
!!!
!!!  (c) 2022 Blaise Bourdin bourdin@lsu.edu
!!!      2025 Blaise Bourdin bourdin@mcmaster.ca
!!!

   type, abstract, extends(MEF90Object) :: MEF90DefMechAT_Type
      PetscReal                                        :: fractureToughness = 0.0_Kr
      class(mef90Mat), allocatable                     :: toughnessAnisotropyMatrix
      PetscReal                                        :: internalLength = -1.0_Kr
      PetscBool                                        :: isElastic = PETSC_FALSE
      PetscReal                                        :: cw = -1.0_Kr
      PetscInt                                         :: aOrder = -1
      PetscInt                                         :: wOrder = -1
      character(len=MEF90MXSTRLEN)                     :: type = ''

      contains
      ! procedure, pass(self) :: setFromOptions => MEF90DefMechAT_setFromOptions
      procedure, pass(self) :: view => MEF90DefMechAT_View
      procedure(ATInterface), pass(self), deferred     :: a
      procedure(ATInterface), pass(self), deferred     :: Da
      procedure(ATInterface), pass(self), deferred     :: D2a
      procedure(ATInterface), pass(self), deferred     :: w
      procedure(ATInterface), pass(self), deferred     :: Dw
      procedure(ATInterface), pass(self), deferred     :: D2w
   end type MEF90DefMechAT_Type

   abstract interface
      PetscReal function ATInterface(self, alpha)
         use petscsys
         import :: MEF90DefMechAT_Type
         class(MEF90DefMechAT_Type), intent(IN)        :: self
         PetscReal                                     :: alpha
      end function ATInterface
   end interface

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechAT_view"
!!!
!!!
!!!  MEF90DefMechAT_view: the default viewer for a MEF90_DefMechAT_Type
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90DefMechAT_View(self,viewer,ierr)
      class(MEF90DefMechAT_Type), intent(in) :: self
      type(tPetscViewer), intent(in) :: viewer
      PetscErrorCode, intent(inout) :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char) :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char) :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType == 'ascii') then
         write(IOBuffer, "(A,': Options for MEF90DefMechAT_type\n')") trim(self%prefix) // "_damage"
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         type: ',A,'\n')") trim(self%type)
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         fracture toughness (Gc): ',ES12.5,' [N.m^(-1)]\n')") self%fractureToughness
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         select type (k => self%toughnessAnisotropyMatrix)
         type is (MatS2D)
            write(IOBuffer, "('         fracture toughness anisotropy matrix (): ', 2(ES12.5,', '), ES12.5, ' []\n')") k
         type is (MatS3D)
            write(IOBuffer, "('         fracture toughness anisotropy matrix (): ', 5(ES12.5,', '), ES12.5, ' []\n')") k
         class default
            write(IOBuffer, *) 'somehow fracture toughness anisotropy matrix is neither MatS2D not MatS3D. This is wrong\n'
         end select
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         internalLength (ell): ',ES12.5,' [m]\n')") self%internalLength
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         elastic (): ',L1,' [bool]\n')") self%isElastic
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
      end if
   end subroutine MEF90DefMechAT_View

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechAT_setFromOptions"
!!!
!!!
!!!  MEF90DefMechAT_setFromOptions: initializes a MEF90_DefMechAT_Type from options
!!!                                 subclasses may implement their own setFromOptions if they have additional 
!!!                                 parameters, but it is expected that these will call MEF90DefMechAT_setFromOptions
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90DefMechAT_setFromOptions(self,ierr)
      class(MEF90DefMechAT_Type), intent(inout) :: self
      PetscErrorCode,intent(inout) :: ierr

      PetscBool :: printHelp
      PetscReal, dimension(:), allocatable :: tmpArray
      PetscInt :: nOpt

      self%internalLength = 1.0e-2

      PetscCall(PetscOptionsBegin(self%comm, trim(self%prefix) // "damage_", "Options for MEF90DefMechAT_type", "mef90DefMech", ierr))
         PetscCall(PetscOptionsReal('-internalLength', 'internal length (\ell)', '[m]', self%internalLength, self%internalLength, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-fractureToughness', 'fracture toughness (G_c))', '[N.m^(-1)]', self%fractureToughness, self%fractureToughness, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-elastic', 'set is elastic', '', PETSC_FALSE, self%isElastic, PETSC_NULL_BOOL, ierr))
         select type (k => self%toughnessAnisotropyMatrix)
         type is (MatS2D)
            nOpt = 3
            allocate(tmpArray(nOpt))
            tmpArray = k !MEF90MatS2DIdentity
            PetscCallA(PetscOptionsRealArray('-toughnessAnisotropyMatrix', 'toughness anisotropy matrix', '[]', tmpArray, nOpt, PETSC_NULL_BOOL, ierr))
            k = tmpArray
            deallocate(tmpArray)
         type is (MatS3D)
            nOpt = 6
            allocate(tmpArray(nOpt))
            tmpArray = k !MEF90MatS2DIdentity
            PetscCallA(PetscOptionsRealArray('-toughnessAnisotropyMatrix', 'toughness anisotropy matrix', '[]', tmpArray, nOpt, PETSC_NULL_BOOL, ierr))
            k = tmpArray
            deallocate(tmpArray)
         end select
      PetscCall(PetscOptionsEnd(ierr))

      PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", printHelp, PETSC_NULL_BOOL, ierr))
      if (printHelp) then
         call self%view(PETSC_VIEWER_STDOUT_WORLD,ierr)
      end if
   end subroutine MEF90DefMechAT_setFromOptions
end module m_MEF90_DefMechAT_class

