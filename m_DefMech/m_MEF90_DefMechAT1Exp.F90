#include "../MEF90/mef90.inc"
module m_MEF90_DefMechAT1exp
#include "petsc/finclude/petsc.h"
   use petscsys
   use m_MEF90_Parameters
   use m_MEF90_LinAlg

   use m_MEF90_DefMechAT_class
   use iso_c_binding
   implicit none(type)
   private
   public :: MEF90DefMechAT1exp_Type

!!! AT1exp, a variant of AT1 model with an exponential stiffness interpolation
!!! function:
!!!
!!! a_b(s)  = 1 + (e^{-bs} - 1) / (1 - e^-b) if b /= 0
!!! a_0(s)  = 1-s
!!!
!!! a is convex if b > 1 and
!!! a_b'(0) = -b / (1-e^{-b}) < -2 if b < 1.5
!!!
   type, extends(MEF90DefMechAT_Type) :: MEF90DefMechAT1exp_Type
      PetscReal :: b = 1.0_kr
   contains
      procedure, pass(self) :: a => aAT1exp
      procedure, pass(self) :: Da => DaAT1exp
      procedure, pass(self) :: D2a => D2aAT1exp

      procedure, pass(self) :: w => wAT1exp
      procedure, pass(self) :: Dw => DwAT1exp
      procedure, pass(self) :: D2w => D2wAT1exp

      procedure, pass(self) :: setFromOptions => MEF90DefMechAT1exp_setFromOptions
      procedure, pass(self) :: view => MEF90DefMechAT1exp_view
   end type MEF90DefMechAT1exp_Type

contains

#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechAT1exp_view"
!!!
!!!
!!!  MEF90DefMechAT1exp_view: the default viewer for a MEF90_DefMechAT_Type
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90DefMechAT1exp_View(self,viewer,ierr)
      class(MEF90DefMechAT1exp_Type), intent(in) :: self
      type(tPetscViewer), intent(in)             :: viewer
      PetscErrorCode, intent(inout)              :: ierr

      character(len=MEF90MXSTRLEN, kind=c_char)  :: IOBuffer
      character(len=MEF90MXSTRLEN, kind=c_char)  :: viewerType

      PetscCall(MEF90DefMechAT_View(self, viewer, ierr))
      PetscCall(PetscViewerGetType(viewer, viewerType, ierr))
      if (viewerType == 'ascii') then
         write(IOBuffer, "(A,': Options for MEF90DefMechAT1exp_type\n')") trim(self%prefix) // "damage_AT1exp"
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         b : ',ES12.5,' [N.m^(-1)]\n')") self%b
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
         write(IOBuffer, "('         cw : ',ES12.5,' [N.m^(-1)]\n')") self%cw
         PetscCall(PetscViewerASCIIPrintf(viewer, IOBuffer, ierr))  
      end if
   end subroutine MEF90DefMechAT1exp_View


#undef __FUNCT__
#define __FUNCT__ "MEF90DefMechAT1exp_setFromOptions"
!!!
!!!
!!!  MEF90DefMechAT1exp_setFromOptions: initializes a MEF90_DefMechAT1exp_Type from options
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90DefMechAT1exp_setFromOptions(self,ierr)
      class(MEF90DefMechAT1exp_Type), intent(inout) :: self
      PetscErrorCode,intent(inout)                  :: ierr


      self%b = 1.0_kr
      self%cw = 0.597674895613893_Kr
      PetscCall(PetscOptionsBegin(self%comm, trim(self%prefix) // "damage_AT1exp_", "Options for MEF90DefMechATexp_type", "mef90DefMech", ierr))
         PetscCall(PetscOptionsReal('-b', 'b parameter', '[]', self%b, self%b, PETSC_NULL_BOOL, ierr))
         PetscCall(PetscOptionsReal('-cw', 'cw parameter', '[]', self%cw, self%cw, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))

      self%aorder = 3
      self%worder = 1
      self%type = 'MEF90DefMechAT1exp'
      
      PetscCall(MEF90DefMechAT_setFromOptions(self, ierr))
   end subroutine MEF90DefMechAT1exp_setFromOptions

#undef __FUNCT__
#define __FUNCT__ "aAT1exp"
!!!
!!!
!!!  aAT1exp: the "a" function of the standard AT1exp model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function aAT1exp(self, alpha)
      class(MEF90DefMechAT1exp_Type), intent(IN) :: self
      PetscReal                                  :: alpha

      if (self%b == 0.0_kr) then
         aAT1exp = 1.0_kr - alpha
      else
         aAT1exp = 1.0_kr + (exp(-self%b * alpha) - 1.0_kr) / (1.0_kr - exp(-self%b))
      end if
   end function aAT1exp

#undef __FUNCT__
#define __FUNCT__ "DaAT1exp"
!!!
!!!
!!!  DaAT1exp: the derivative of the "a" function of the standard AT1exp model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DaAT1exp(self, alpha)
      class(MEF90DefMechAT1exp_Type), intent(IN) :: self
      PetscReal                                  :: alpha

      if (self%b == 0.0_kr) then
         DaAT1exp = -1.0_kr
      else
         DaAT1exp = -self%b * exp(-self%b * alpha) / (1.0_kr - exp(-self%b))
      end if
   end function DaAT1exp

#undef __FUNCT__
#define __FUNCT__ "D2aAT1exp"
!!!
!!!
!!!  D2aAT1exp: the second derivative of the "a" function of the standard AT1exp model, i.e. a(\alpha) = (1-\alpha)^2
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2aAT1exp(self, alpha)
      class(MEF90DefMechAT1exp_Type), intent(IN) :: self
      PetscReal                                  :: alpha

      if (self%b == 0.0_kr) then
         D2aAT1exp = 0.0_kr
      else
         D2aAT1exp = self%b**2 * exp(-self%b * alpha) / (1.0_kr - exp(-self%b))
      end if
   end function D2aAT1exp

#undef __FUNCT__
#define __FUNCT__ "wAT1exp"
!!!
!!!
!!!  wAT1exp: the "w" function of the standard AT1exp model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function wAT1exp(self, alpha)
      class(MEF90DefMechAT1exp_Type), intent(IN) :: self
      PetscReal                                  :: alpha

      wAT1exp = alpha
   end function wAT1exp

#undef __FUNCT__
#define __FUNCT__ "DwAT1exp"
!!!
!!!
!!!  DwAT1exp: the derivative of the "w" function of the standard AT1exp model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function DwAT1exp(self, alpha)
      class(MEF90DefMechAT1exp_Type), intent(IN) :: self
      PetscReal                                  :: alpha

      DwAT1exp = 1.0_kr
   end function DwAT1exp

#undef __FUNCT__
#define __FUNCT__ "D2wAT1exp"
!!!
!!!
!!!  D2wAT1exp: the second derivative of the "w" function of the standard AT1exp model, i.e. w(\alpha) = \alpha
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!
   PetscReal function D2wAT1exp(self, alpha)
      class(MEF90DefMechAT1exp_Type), intent(IN) :: self
      PetscReal                                  :: alpha

      D2wAT1exp = 0.0_kr
   end function D2wAT1exp
end module m_MEF90_DefMechAT1exp
