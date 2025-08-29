
module m_MEF90_HookesLaw
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_HookesLaw_class
   use m_MEF90_HookesLawIsotropic2D
   use m_MEF90_HookesLawIsotropic3D
   use petscsys
   implicit none(type)

   private
   public :: MEF90GetHookesLaw
   public :: MEF90HookesLaw_Type
   public :: MEF90HookesLawIsotropic2D_Type
   public :: MEF90HookesLawIsotropic3D_Type


   enum, bind(c)
      enumerator :: MEF90HookesLaw_TypeIsotropic = 0, &
         MEF90HookesLaw_TypeFull
   end enum

   character(len=MEF90MXSTRLEN), dimension(5), protected   :: MEF90HookesLaw_TypeList = [ &
      'Isotropic           ', &
      'Full                ', &
      'MEF90HookesLaw_Type ', &
      '_MEF90HookesLaw_Type', &
      '                    '  &
   ]

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90GetHookesLaw"
!!!
!!!
!!!  MEF90GetHookesLaw: Returns a MEF90HookesLaw_Type from set options
!!!
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   subroutine MEF90GetHookesLaw(comm, prefix, dim, HookesLaw, ierr)
      MPI_Comm, intent(in) :: comm
      character(len = MEF90MXSTRLEN), intent(in) :: prefix
      PetscInt, intent(IN) :: dim
      class(MEF90HookesLaw_Type), allocatable, intent(out) :: HookesLaw
      PetscErrorCode, intent(inout) :: ierr

      PetscEnum :: HookesLawType

      HookesLawType = MEF90HookesLaw_TypeIsotropic
      PetscCall(PetscOptionsBegin(comm, trim(prefix), "Options for MEF90HookesLaw_type", "MEF90", ierr))
         PetscCall(PetscOptionsEnum("-HookesLaw_type", "Hookes law type", "MEF90", MEF90HookesLaw_TypeList, MEF90HookesLaw_TypeIsotropic, HookesLawType, PETSC_NULL_BOOL, ierr))
      PetscCall(PetscOptionsEnd(ierr))


      select case (HookesLawType)
         case (MEF90HookesLaw_TypeIsotropic)
            select case(dim)
               case(2)
                  HookesLaw = MEF90HookesLawIsotropic2D_type(comm = comm, prefix = prefix)
                  !!! I think that I need to use comm = comm, prefix = prefix because MEF90HookesLawIsotropic2D_type extends 
                  !!! MEF90Object_Type, so that the order of positional arguments is not clear.
               case(3)
                  HookesLaw = MEF90HookesLawIsotropic3D_type(comm = comm, prefix = prefix)
               case default
                  write(*,"('In ', A, ' dim = ', I2, 'and should be 2 or 3')") __FUNCT__, dim
                  STOP
            end select
      end select
   end subroutine MEF90GetHookesLaw
end module m_MEF90_HookesLaw
