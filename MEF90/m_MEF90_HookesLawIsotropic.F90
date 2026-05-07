module m_MEF90_HookesLawIsotropic2D
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_LinAlg
   use m_MEF90_HookesLaw_Class
   use iso_c_binding

   implicit none(type)
   private
   public :: MEF90HookesLawIsotropic2D

   type, extends(MEF90HookesLaw) :: MEF90HookesLawIsotropic2D
      PetscReal :: YoungsModulus = 1.0_Kr
      PetscReal :: PoissonRatio = 0.3_Kr
      PetscReal :: lambda = 0.0_Kr 
      PetscReal :: mu = 0.0_Kr
      PetscReal :: BulkModulus = 0.0_Kr
      PetscBool :: isPlaneStress = PETSC_FALSE
   contains
         procedure :: setFromOptions => MEF90HookesLawIsotropic2D_setFromOptions
         procedure :: view => MEF90HookesLawIsotropic2D_view
         procedure :: mult => MEF90HookesLawIsotropic2D_mult
         procedure :: multmult => MEF90HookesLawIsotropic2D_multmult
   end type MEF90HookesLawIsotropic2D

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLawIsotropic2D_setFromOptions"
!!!
!!!
!!!  MEF90HookesLawIsotropic2D_setFromOptions: initializes a MEF90HookesLawIsotropic2D from options
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HookesLawIsotropic2D_setFromOptions(self, ierr)
      class(MEF90HookesLawIsotropic2D), intent(inout) :: self
      PetscErrorCode,intent(inout)                    :: ierr

      PetscBool                                       :: printHelp

      PetscCall(PetscOptionsBegin(self%comm, trim(self%prefix) // "HookesLaw_Isotropic_", "Options for MEF90HookesLawIsotropic2D_Type", "mef90HookesLaw", ierr))
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

subroutine MEF90HookesLawIsotropic2D_View(self, viewer, ierr)
      class(MEF90HookesLawIsotropic2D), intent(in) :: self
      type(tPetscViewer), intent(in)               :: viewer
      PetscErrorCode, intent(inout)                :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)    :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)    :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType == 'ascii') then
         write(IOBuffer, "(A,': Options for MEF90HookesLaw\n')") trim(self%prefix)//"HookesLaw"
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Type: Isotropic\n')")
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

   subroutine MEF90HookesLawIsotropic2D_mult(A, phi, Aphi, ierr)
      class(MEF90HookesLawIsotropic2D), intent(in) :: A
      class(mef90Mat), intent(in)                  :: phi
      class(mef90Mat), allocatable, intent(out)    :: Aphi
      character(len=MEF90MXSTRLEN, kind=c_char)    :: IOBuffer
      PetscErrorCode, intent(inout)                :: ierr

      select type (phi2D => phi)
         type is (MatS2D)
            Aphi = A%lambda * trace(phi2D) * MEF90MatS2DIdentity + 2.0_Kr * A%mu * phi2D
         class default
            write (IOBuffer, *) "Incompatible arguments in "//__FUNCT__//'\n'
            PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
      end select ! phi
   end subroutine MEF90HookesLawIsotropic2D_mult

   subroutine MEF90HookesLawIsotropic2D_multmult(A, phi, psi, Aphipsi, ierr)
      class(MEF90HookesLawIsotropic2D), intent(in) :: A
      class(mef90Mat), intent(in)                  :: phi, psi
      PetscReal, intent(out)                       :: Aphipsi
      character(len=MEF90MXSTRLEN, kind=c_char)    :: IOBuffer
      PetscErrorCode, intent(inout)                :: ierr

      !!! This is really absurd
      select type (phi2D => phi)
         type is (MatS2D)
         select type (psi2D => psi)
            type is (MatS2D)
               Aphipsi = A%lambda * trace(phi2D) * trace(psi2D) &
                         + 2.0_Kr * A%mu * (phi2D .dotP. psi2D)
            class default
               write (IOBuffer, *) "Incompatible arguments in "//__FUNCT__//'\n'
               PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
               SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
            end select ! psi
         class default
            write (IOBuffer, *) "Incompatible arguments in "//__FUNCT__//'\n'
            PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
      end select ! phi
   end subroutine MEF90HookesLawIsotropic2D_multmult
end module m_MEF90_HookesLawIsotropic2D

module m_MEF90_HookesLawIsotropic3D
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_LinAlg
   use m_MEF90_HookesLaw_Class
   use iso_c_binding

   implicit none(type)
   private

   public :: MEF90HookesLawIsotropic3D

   type, extends(MEF90HookesLaw) :: MEF90HookesLawIsotropic3D
      PetscReal :: YoungsModulus = 1.0_Kr
      PetscReal :: PoissonRatio = 0.3_Kr
      PetscReal :: lambda = 0.0_Kr 
      PetscReal :: mu = 0.0_Kr
      PetscReal :: BulkModulus = 0.0_Kr
   contains
         procedure, pass(self) :: setFromOptions => MEF90HookesLawIsotropic3D_setFromOptions
         procedure, pass(self) :: view => MEF90HookesLawIsotropic3D_view
         procedure             :: mult => MEF90HookesLawIsotropic3D_mult
         procedure             :: multmult => MEF90HookesLawIsotropic3D_multmult
   end type MEF90HookesLawIsotropic3D

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLawIsotropic3D_setFromOptions"
!!!
!!!
!!!  MEF90HookesLawIsotropic3D_setFromOptions: initializes a MEF90HookesLawIsotropic3D_type from options
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90HookesLawIsotropic3D_setFromOptions(self,ierr)
      class(MEF90HookesLawIsotropic3D), intent(inout) :: self
      PetscErrorCode,intent(inout)                    :: ierr

      PetscBool :: printHelp
      PetscCall(PetscOptionsBegin(self%comm, trim(self%prefix) // "HookesLaw_Isotropic_", "Options for MEF90HookesLawIsotropic3D_Type", "mef90HookesLaw", ierr))
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
      class(MEF90HookesLawIsotropic3D), intent(in) :: self
      type(tPetscViewer), intent(in)               :: viewer
      PetscErrorCode, intent(inout)                :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)    :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)    :: viewerType

      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType == 'ascii') then
         write(IOBuffer, "(A,': Options for MEF90HookesLaw\n')") trim(self%prefix)//"HookesLaw"
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         Type: Isotropic\n')")
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

   subroutine MEF90HookesLawIsotropic3D_mult(A, phi, Aphi, ierr)
      class(MEF90HookesLawIsotropic3D), intent(in) :: A
      class(mef90Mat), intent(in)                  :: phi
      class(mef90Mat), allocatable, intent(out)    :: Aphi
      character(len=MEF90MXSTRLEN, kind=c_char)    :: IOBuffer
      PetscErrorCode, intent(inout)                :: ierr

      select type (phi3D => phi)
         type is (MatS3D)
            Aphi = A%lambda * trace(phi3D) * MEF90MatS3DIdentity + 2.0_Kr * A%mu * phi3D
         class default
            write (IOBuffer, *) "Incompatible arguments in "//__FUNCT__//'\n'
            PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
      end select ! phi
   end subroutine MEF90HookesLawIsotropic3D_mult

   subroutine MEF90HookesLawIsotropic3D_multmult(A, phi, psi, Aphipsi, ierr)
      class(MEF90HookesLawIsotropic3D), intent(in) :: A
      class(mef90Mat), intent(in)                  :: phi, psi
      PetscReal, intent(out)                       :: Aphipsi
      character(len=MEF90MXSTRLEN, kind=c_char)    :: IOBuffer
      PetscErrorCode, intent(inout)                :: ierr

      !!! This is really absurd
      select type (phi3D => phi)
         type is (MatS3D)
         select type (psi3D => psi)
            type is (MatS3D)
               Aphipsi = A%lambda * trace(phi3D) * trace(psi3D) &
                         + 2.0_Kr * A%mu * (phi3D .dotP. psi3D)
            class default
               write (IOBuffer, *) "Incompatible arguments in "//__FUNCT__//'\n'
               PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
               SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
            end select ! psi
         class default
            write (IOBuffer, *) "Incompatible arguments in "//__FUNCT__//'\n'
            PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
      end select ! phi
   end subroutine MEF90HookesLawIsotropic3D_multmult
end module m_MEF90_HookesLawIsotropic3D
