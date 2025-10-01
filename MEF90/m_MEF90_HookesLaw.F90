
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
   public :: MEF90HookesLawSum


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

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLawSum"
!!!
!!!
!!!  MEF90HookesLawSum: adding 2 MEF90HookesLaw_Type objects
!!!
!!!  (c) 2025 Blaise Bourdin bourdin@mcmaster.ca
!!!
   function MEF90HookesLawSum(A, B)
      class(MEF90HookesLaw_Type), intent(in) :: A, B
      class(MEF90HookesLaw_Type), allocatable :: MEF90HookesLawSum

      type(MEF90HookesLawIsotropic2D_Type) :: MEF90HookesLawSumIso2D
      character(len=MEF90MXSTRLEN) :: IOBuffer
      PetscErrorCode  :: ierr

      select type (AA => A)
      type is (MEF90HookesLawIsotropic2D_Type)
         select type (BB => B)
         type is (MEF90HookesLawIsotropic2D_Type)
            if (AA%isPlaneStress .eqv. BB%isPlaneStress) then
               MEF90HookesLawSumIso2D = AA
               MEF90HookesLawSumIso2D%isPlaneStress = AA%isPlaneStress
               MEF90HookesLawSumIso2D%lambda = AA%lambda + BB% lambda
               MEF90HookesLawSumIso2D%mu = AA%mu + BB% mu
               MEF90HookesLawSumIso2D%BulkModulus = MEF90HookesLawSumIso2D%lambda + MEF90HookesLawSumIso2D%mu
               if (AA%isPlaneStress) then
                  MEF90HookesLawSumIso2D%PoissonRatio = MEF90HookesLawSumIso2D%lambda / (MEF90HookesLawSumIso2D%lambda + MEF90HookesLawSumIso2D%mu) * 0.5_kr
                  MEF90HookesLawSumIso2D%YoungsModulus = 2.0_kr * MEF90HookesLawSumIso2D%mu * (1.0_kr + MEF90HookesLawSumIso2D%PoissonRatio)
               else
                  MEF90HookesLawSumIso2D%PoissonRatio = MEF90HookesLawSumIso2D%lambda / (MEF90HookesLawSumIso2D%lambda + 2.0_kr * MEF90HookesLawSumIso2D%mu) * 0.5_kr
                  MEF90HookesLawSumIso2D%YoungsModulus = 2.0_kr * MEF90HookesLawSumIso2D%mu * (1.0_kr + MEF90HookesLawSumIso2D%PoissonRatio)
               end if
               MEF90HookesLawSum = MEF90HookesLawSumIso2D
            else
               write (IOBuffer, *) "Incompatible planar Hooke laws in "//__FUNCT__//'\n'
               PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
               SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
            end if
         class default
            write (IOBuffer, *) "Incompatible Hooke law types in "//__FUNCT__//'\n'
            PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
         end select
      class default
         write (IOBuffer, *) "Unimplemented Hooke law types in "//__FUNCT__//'\n'
         PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
      end select
   end function MEF90HookesLawSum


end module m_MEF90_HookesLaw
