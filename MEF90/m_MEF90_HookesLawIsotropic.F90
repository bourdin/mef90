module m_MEF90_HookesLawIsotropic2D
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_LinAlg
   use m_MEF90_HookeLaw_Class
   use iso_c_binding

   implicit none(type)
   private
   public :: MEF90HookesLawIsotropic2D_type

   type, extends(MEF90HookesLaw_Type) :: MEF90HookesLawIsotropic2D_type
      PetscReal :: YoungsModulus = 1.0_Kr
      PetscReal :: PoissonRatio = 0.3_Kr
      PetscReal :: lambda = 0.0_Kr 
      PetscReal :: mu = 0.0_Kr
      PetscReal :: BulkModulus = 0.0_Kr
      PetscBool :: isPlaneStress = PETSC_FALSE
   contains
         procedure, pass(self) :: setFromOptions => MEF90HookesLawIsotropic2D_setFromOptions
         procedure, pass(self) :: view => MEF90HookesLawIsotropic2D_view
   end type MEF90HookesLawIsotropic2D_type

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLawIsotropic2D_setFromOptions"
!!!
!!!
!!!  MEF90HookesLawIsotropic2D_setFromOptions: initializes a MEF90HookesLawIsotropic2D_type from options
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90HookesLawIsotropic2D_setFromOptions(self,ierr)
      class(MEF90HookesLawIsotropic2D_Type), intent(inout) :: self
      PetscErrorCode,intent(inout) :: ierr

      PetscBool :: printHelp
      PetscCall(PetscOptionsBegin(self%comm, trim(self%prefix) // "HookesLaw_", "Options for MEF90HookesLawIsotropic2D_Type", "mef90HookesLaw", ierr))
         PetscCall(PetscOptionsReal('-YoungsModulus', 'Young''s modulus (E)', '[Pa]', self%YoungsModulus, self%YoungsModulus, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-PoissonRatio', 'Poisson ratio (\nu))', '[]', self%PoissonRatio, self%PoissonRatio, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsBool('-planeStress', '2D plane stress', '', PETSC_FALSE, self%isPlaneStress, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))
      if (self%isPlaneStress) then
         self%lambda = self%YoungsModulus * self%PoissonRatio / (1.0_kr - self%PoissonRatio**2)
         self%mu = self%YoungsModulus / (1.0_kr + self%PoissonRatio)*.5_kr
         self%BulkModulus = self%lambda + 2.0_kr * self%mu / 3.0_kr
      else
         self%lambda = self%YoungsModulus * self%PoissonRatio / (1.0_kr + self%PoissonRatio) / (1.0_kr - 2.0_kr * self%PoissonRatio)
         self%mu = self%YoungsModulus / (1.0_kr + self%PoissonRatio)*.5_kr
         self%BulkModulus = self%lambda + self%mu
      end if

      PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", printHelp, PETSC_NULL_BOOL, ierr))
      if (printHelp) then
         call self%view(PETSC_VIEWER_STDOUT_WORLD,ierr)
      end if
   end subroutine MEF90HookesLawIsotropic2D_setFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLawIsotropic2D_view"
!!!
!!!
!!!  MEF90HookesLawIsotropic2D_view: the default viewer for a MEF90_DefMechAT_Type
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90HookesLawIsotropic2D_View(self,viewer,ierr)
      class(MEF90HookesLawIsotropic2D_Type), intent(in) :: self
      type(tPetscViewer), intent(in) :: viewer
      PetscErrorCode, intent(inout) :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char) :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char) :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType == 'ascii') then
         write(IOBuffer, "(A,': Options for MEF90HookesLawIsotropic2D_type\n')") trim(self%prefix)
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Youngs modulus (E): ',ES12.5,' [Pa]\n')") self%YoungsModulus
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Poisson ratio (nu): ',ES12.5,' []\n')") self%PoissonRatio
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Bulk modulus (\kappa): ',ES12.5,' [Pa]\n')") self%BulkModulus
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Lame coefficient (\lambda): ',ES12.5,' [Pa]\n')") self%lambda
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Shear modulus (\mu): ',ES12.5,' [Pa]\n')") self%mu
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         plane stress (): ',L1,' [bool]\n')") self%isPlaneStress
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
      end if
   end subroutine MEF90HookesLawIsotropic2D_View

end module m_MEF90_HookesLawIsotropic2D

module m_MEF90_HookesLawIsotropic3D
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_LinAlg
   use m_MEF90_HookeLaw_Class
   use iso_c_binding

   implicit none(type)
   private

   public :: MEF90HookesLawIsotropic3D_type

   type, extends(MEF90HookesLaw_Type) :: MEF90HookesLawIsotropic3D_type
      PetscReal :: YoungsModulus = 1.0_Kr
      PetscReal :: PoissonRatio = 0.3_Kr
      PetscReal :: lambda = 0.0_Kr 
      PetscReal :: mu = 0.0_Kr
      PetscReal :: BulkModulus = 0.0_Kr
   contains
         procedure, pass(self) :: setFromOptions => MEF90HookesLawIsotropic3D_setFromOptions
         procedure, pass(self) :: view => MEF90HookesLawIsotropic3D_view
   end type MEF90HookesLawIsotropic3D_type

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLawIsotropic3D_setFromOptions"
!!!
!!!
!!!  MEF90HookesLawIsotropic3D_setFromOptions: initializes a MEF90HookesLawIsotropic3D_type from options
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90HookesLawIsotropic3D_setFromOptions(self,ierr)
      class(MEF90HookesLawIsotropic3D_Type), intent(inout) :: self
      PetscErrorCode,intent(inout) :: ierr

      PetscBool :: printHelp
      PetscCall(PetscOptionsBegin(self%comm, trim(self%prefix) // "HookesLaw_", "Options for MEF90HookesLawIsotropic3D_Type", "mef90HookesLaw", ierr))
         PetscCall(PetscOptionsReal('-YoungsModulus', 'Young''s modulus (E)', '[Pa]', self%YoungsModulus, self%YoungsModulus, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-PoissonRatio', 'Poisson ratio (\nu))', '[]', self%PoissonRatio, self%PoissonRatio, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))
      self%lambda = self%YoungsModulus * self%PoissonRatio / (1.0_kr + self%PoissonRatio) / (1.0_kr - 2.0_kr * self%PoissonRatio)
      self%mu = self%YoungsModulus / (1.0_kr + self%PoissonRatio)*.5_kr
      self%BulkModulus = self%lambda + 2.0_kr * self%mu / 3.0_kr

      PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", printHelp, PETSC_NULL_BOOL, ierr))
      if (printHelp) then
         call self%view(PETSC_VIEWER_STDOUT_WORLD,ierr)
      end if
   end subroutine MEF90HookesLawIsotropic3D_setFromOptions

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLawIsotropic3D_view"
!!!
!!!
!!!  MEF90HookesLawIsotropic3D_view: the default viewer for a MEF90_DefMechAT_Type
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90HookesLawIsotropic3D_View(self,viewer,ierr)
      class(MEF90HookesLawIsotropic3D_Type), intent(in) :: self
      type(tPetscViewer), intent(in) :: viewer
      PetscErrorCode, intent(inout) :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char) :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char) :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType == 'ascii') then
         write(IOBuffer, "(A,': Options for MEF90HookesLawIsotropic3D_type\n')") trim(self%prefix)
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Youngs modulus (E): ',ES12.5,' [Pa]\n')") self%YoungsModulus
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Poisson ratio (nu): ',ES12.5,' []\n')") self%PoissonRatio
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Bulk modulus (\kappa): ',ES12.5,' [Pa]\n')") self%BulkModulus
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Lame coefficient (\lambda): ',ES12.5,' [Pa]\n')") self%lambda
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Shear modulus (\mu): ',ES12.5,' [Pa]\n')") self%mu
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
      end if
   end subroutine MEF90HookesLawIsotropic3D_View

end module m_MEF90_HookesLawIsotropic3D
