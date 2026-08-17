module m_MEF90_baseClass
#include "petsc/finclude/petsc.h"
   use petscsys
   use m_MEF90_Parameters
   use, intrinsic :: iso_c_binding
   implicit none(type, external)

   private
   public :: MEF90Object

!!! author: Blaise Bourdin (2025, bourdin@mcmaster.ca)
!!!
!!!  MEF90_baseClass: The base class for all MEF90 classes
!!!
!!!
!!!  NOTE FOR CLASSES THAT ARE HANDED TO PETSc AS AN APPLICATION CONTEXT
!!!
!!!  Extending MEF90Object gives a type deferred bindings, hence type-bound procedures, and Fortran
!!!  forbids passing a derived type that has type-bound or final procedures to an assumed-type dummy
!!!  argument. All of the PETSc Fortran context arguments are assumed-type,
!!!
!!!     subroutine SNESSetFunction(snes, r, f, ctx, ierr)
!!!        type(*) :: ctx
!!!
!!!  so an extension of MEF90Object can no longer be passed to SNESSetFunction, SNESSetJacobian,
!!!  TSSetIFunction, TAOSetObjective, ... the way a plain derived type could before mef90 0.5.2. The
!!!  compiler rejects it with
!!!
!!!     Error: Actual argument to assumed-type dummy has type parameters or is of derived type with
!!!            type-bound or FINAL procedures
!!!
!!!  The way out, used by MEF90DefMech_Type and MEF90HeatXfer_Type, is to give the class a
!!!
!!!     type(c_ptr) :: PETScCtx = C_NULL_PTR
!!!
!!!  component, set once at creation with PETScCtx = c_loc(self), and to hand THAT to PETSc. The
!!!  callbacks then declare their context argument as type(c_ptr) and recover the object with
!!!  call c_f_pointer(PETScCtx, self).
!!!
!!!  Two properties of that component are load-bearing, and both are easy to break:
!!!
!!!    - It must OUTLIVE every callback. Fortran passes actual arguments by reference, so what PETSc
!!!      records is the address of the argument it was given, not the address that argument contains.
!!!      Passing a local of the routine that registers the callbacks, or an inline c_loc(self) whose
!!!      result lives in a compiler temporary, makes PETSc record a stack slot that is dead by the time
!!!      the solver calls back. That is a use-after-scope, and it only misbehaves once the frame has
!!!      been reused, so it can survive a short test and corrupt a long run. Giving the component the
!!!      lifetime of the object it points to is what makes this safe.
!!!    - It must be ONE PER INSTANCE, which is why a local with the save attribute will not do: a saved
!!!      local is shared by every call of its routine, so creating a second context would silently
!!!      repoint the first solver at it.
!!!
!!!  The cost of the detour is that PETScCtx is an untyped address: c_f_pointer will happily reinterpret
!!!  whatever it is given, so handing a callback the wrong context is not a compile error.
!!!

   type, abstract :: MEF90Object
      MPIU_Comm :: comm
      character(len=MEF90MXSTRLEN, kind=c_char) :: prefix
      character(len=MEF90MXSTRLEN, kind=c_char) :: name = "NULL"
   contains
      procedure, pass(self)                                               :: setFromOptions => MEF90Object_setFromOptions
      procedure(MEF90ObjectViewInterface), pass(self), deferred           :: view_internal
      procedure, pass(self)                                               :: view => MEF90Object_view
   end type MEF90Object

   abstract interface
      subroutine MEF90ObjectViewInterface(self, viewer, ierr)
         use petscsys
         import :: MEF90Object
         class(MEF90Object), intent(in) :: self
         type(tPetscViewer), intent(in) :: viewer
         PetscErrorCode, intent(inout)  :: ierr
      end subroutine MEF90ObjectViewInterface
   end interface
contains

#undef __FUNCT__
#define __FUNCT__ "MEF90Object_view"
!!! author: Blaise Bourdin (2026, bourdin@mcmaster.ca)
!!!
!!!  MEF90Object_view: Views a MEF90Object and flagged it as viewed so that it is not viewed multiple times when -verbose is set to a high value
!!!

subroutine MEF90Object_view(self, viewer, ierr)
      class(MEF90Object), intent(in) :: self
      type(tPetscViewer), intent(in) :: viewer
      PetscErrorCode, intent(inout)  :: ierr

      PetscBool                      :: helpPrinted, flg

      PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-" // trim(self%name) // "_viewed", helpPrinted, PETSC_NULL_BOOL, ierr))
      if (.not. helpPrinted) then
         call self%view_internal(viewer, ierr)
         PetscCall(PetscOptionsSetValue(PETSC_NULL_OPTIONS, "-" //trim(self%name) // "_viewed", "true", ierr))
         ! This is just to silence a warning about the option being set but not used. The option is used to prevent printing the same object multiple times when -verbose is set to a high value.
         PetscCall(PetscOptionsGetBool(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-" // trim(self%name) // "_viewed", flg, PETSC_NULL_BOOL, ierr))
      end if
   end subroutine MEF90Object_view

! contains
#undef __FUNCT__
#define __FUNCT__ "MEF90Object_setFromOptions"
!!! author: Blaise Bourdin (2026, bourdin@mcmaster.ca)
!!!
!!!  MEF90Object_setFromOptions: generic version of the setFromOptions method for all MEF90Objects, which just prints the object if -verbose is set to a high value
!!!  Can be overridden in classes if they have options to set.
!!!

   subroutine MEF90Object_setFromOptions(self, ierr)
      class(MEF90Object), intent(inout)    :: self
      PetscErrorCode, intent(inout)        :: ierr

      PetscViewer                          :: stdoutViewer
      PetscInt                             :: verbose = 0

      PetscCall(PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-verbose", verbose, PETSC_NULL_BOOL, ierr))
      if (verbose > 0) then
         PetscCall(PetscViewerASCIIGetStdout(self%comm, stdoutViewer, ierr))
         call self%view(stdoutViewer, ierr)
      end if
   end subroutine MEF90Object_setFromOptions
end module m_MEF90_baseClass
