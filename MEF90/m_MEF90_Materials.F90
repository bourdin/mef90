module m_MEF90_Materials_Types
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_LinAlg
   implicit none(type, external)

   type MEF90HookesLaw2D
      type(Tens4OS2D)    :: fullTensor
      type(Tens4OS3D)    :: fullTensorLocal, fullTensor3D
      PetscReal          :: lambda, mu, YoungsModulus, PoissonRatio, BulkModulus
      PetscEnum          :: type
      PetscBool          :: isPlaneStress
   end type MEF90HookesLaw2D

   type MEF90HookesLaw3D
      type(Tens4OS3D)    :: fullTensor, fullTensorLocal
      PetscReal          :: lambda, mu, YoungsModulus, PoissonRatio, BulkModulus
      PetscEnum          :: type
#if (PETSC_SIZEOF_INT == 4)
      ! With 4-byte integers, this declared type is 4-bytes shy of being
      ! 8-byte aligned, which can be problematic for arrays of this type.
      ! The MEF90HookesLaw2D has an extra PetscBool that keeps it 8-byte
      ! aligned, so we'll mimic that:
      PetscBool          :: padding = PETSC_FALSE
#endif
   end type MEF90HookesLaw3D

   type MEF90RotationMatrix3D
      ! The Bunge (passive) convention is used. The rotation matrix transforms a vector from the global frame to the local frame.
      ! - rotation of a vector:              [V_local]_i    = [R]_ij . [V_global]_j
      ! - rotation of a 2nd order tensor:    [M_local]_ij   = [R]_ik . [M_global]_kl . [R^T]_lj
      ! - rotation of a fourth order tensor: [C_local]_ijkl = [R]_ip . [R]_jq . [R]_kr . [R]_ls . [C_global]_pqrs
      type(MAT3D)        :: fullTensor
      PetscReal          :: phi1, Phi, phi2
      type(Vect3D)       :: V1, V2, V3
      PetscBool          :: fromEuler
   end type MEF90RotationMatrix3D

   enum, bind(c)
      enumerator :: MEF90HookesLawTypeFull = 0, &
         MEF90HookesLawTypeIsotropic
   end enum
end module m_MEF90_Materials_Types

module m_MEF90_Materials
#include "petsc/finclude/petsc.h"
   use m_MEF90_Parameters
   use m_MEF90_Utils
   use m_MEF90_LinAlg
   use m_MEF90_Elements
   use m_MEF90_Ctx
   use m_MEF90_DMPlex
   use m_MEF90_Materials_Types
   implicit none(type)

   interface operator(+)
      module procedure MEF90HookesLaw2DSum, MEF90HookesLaw3DSum
   end interface

   interface operator(-)
      module procedure MEF90HookesLaw2DDiff, MEF90HookesLaw3DDiff
   end interface

   interface operator(*)
      module procedure MEF90HookesLaw2DXMatS2D, MEF90HookesLaw3DXMatS3D, MEF90HookesLaw2DXMat2D, MEF90HookesLaw3DXMat3D, ScalarXMEF90HookesLaw2D, ScalarXMEF90HookesLaw3D
   end interface

   PetscSizeT, protected   :: sizeofMEF90HookesLaw2D
   PetscSizeT, protected   :: sizeofMEF90HookesLaw3D

   character(len=MEF90MXSTRLEN), dimension(5), protected   :: MEF90HookesLawTypeList

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90MaterialsInitialize_Private"
!!!
!!!
!!!  MEF90MaterialsInitialize_Private:
!!!
!!!  (c) 2013-2014 Blaise Bourdin bourdin@lsu.edu
!!!
   subroutine MEF90MaterialsInitialize_Private(ierr)
      PetscErrorCode, intent(OUT)          :: ierr

      type(MEF90HookesLaw2D), target       :: HookesLaw2D
      type(MEF90HookesLaw3D), target       :: HookesLaw3D
      character(len=1), pointer            :: dummychar(:)
      PetscSizeT                           :: sizeofchar

      PetscCall(PetscDataTypeGetSize(PETSC_CHAR, sizeofchar, ierr))
      sizeofMEF90HookesLaw2D = size(transfer(HookesLaw2D, dummychar)) * sizeofchar
      sizeofMEF90HookesLaw3D = size(transfer(HookesLaw3D, dummychar)) * sizeofchar

      MEF90HookesLawTypeList(1) = 'Full'
      MEF90HookesLawTypeList(2) = 'Isotropic'
      MEF90HookesLawTypeList(3) = 'MEF90HookesLawTypeList'
      MEF90HookesLawTypeList(4) = '_MEF90HookesLawTypeList'
      MEF90HookesLawTypeList(5) = ''
   end subroutine MEF90MaterialsInitialize_Private

!!! Subroutine generating various types of Hooke's laws
#undef __FUNCT__
#define __FUNCT__ "MEF90HookeLawIsoLambdaMu2D"
   subroutine MEF90HookeLawIsoLambdaMu2D(A, lambda, mu)
      type(Tens4OS2D), intent(OUT)         :: A
      PetscReal, intent(IN)                :: lambda, mu
      A = 0.0_kr

      A%XXXX = lambda + 2.0_kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_kr * mu
   end subroutine MEF90HookeLawIsoLambdaMu2D

#undef __FUNCT__
#define __FUNCT__ "MEF90HookeLawIsoEnu2DPlaneStress"
   subroutine MEF90HookeLawIsoEnu2DPlaneStress(A, E, nu)
      PetscReal, intent(IN)                :: E, nu
      type(Tens4OS2D), intent(OUT)         :: A

      PetscReal                            :: Lambda, mu

      lambda = E * nu / (1.0_kr - nu**2)
      mu = E / (1.0_kr + nu)*.5_kr
      A = 0.0_kr
      A%XXXX = lambda + 2.0_kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_kr * mu
   end subroutine MEF90HookeLawIsoEnu2DPlaneStress

#undef __FUNCT__
#define __FUNCT__ "MEF90HookeLawIsoEnu2DPlaneStrain"
   subroutine MEF90HookeLawIsoEnu2DPlaneStrain(A, E, nu)
      PetscReal, intent(IN)                :: E, nu
      type(Tens4OS2D), intent(OUT)         :: A

      PetscReal                            :: Lambda, mu

      lambda = E * nu / (1.0_kr + nu) / (1.0_kr - 2.0_kr * nu)
      mu = E / (1.0_kr + nu)*.5_kr
      A = 0.0_kr
      A%XXXX = lambda + 2.0_kr * mu
      A%XXYY = lambda
      A%XYXY = mu
      A%YYYY = lambda + 2.0_kr * mu
   end subroutine MEF90HookeLawIsoEnu2DPlaneStrain

#undef __FUNCT__
#define __FUNCT__ "MEF90HookeLawIsoLambdaMu3D"
   subroutine MEF90HookeLawIsoLambdaMu3D(A, lambda, mu)
      PetscReal, intent(IN)                :: lambda, mu
      type(Tens4OS3D), intent(OUT)         :: A

      A = 0.0_kr
      A%XXXX = lambda + mu * 2.0_kr
      A%XXYY = lambda
      A%XXZZ = lambda

      A%XYXY = mu

      A%XZXZ = mu

      A%YYYY = lambda + mu * 2.0_kr
      A%YYZZ = lambda

      A%YZYZ = mu

      A%ZZZZ = lambda + mu * 2.0_kr
   end subroutine MEF90HookeLawIsoLambdaMu3D

#undef __FUNCT__
#define __FUNCT__ "MEF90HookeLawIsoENu3D"
   subroutine MEF90HookeLawIsoENu3D(A, E, nu)
      PetscReal, intent(IN)               :: E, nu
      type(Tens4OS3D), intent(OUT)        :: A

      real(Kind=Kr)                       :: Lambda, mu

      lambda = E * nu / (1.0_kr + nu) / (1 - 2.0_kr * nu)
      mu = E / (1.0_kr + nu)*.5_kr

      A = 0.0_kr
      A%XXXX = lambda + mu * 2.0_kr
      A%XXYY = lambda
      A%XXZZ = lambda

      A%XYXY = mu

      A%XZXZ = mu

      A%YYYY = lambda + mu * 2.0_kr
      A%YYZZ = lambda

      A%YZYZ = mu

      A%ZZZZ = lambda + mu * 2.0_kr
   end subroutine MEF90HookeLawIsoENu3D

!!! Overloading linear algebra functions with Hookes Laws.
#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw2DSum"
!!!
!!!
!!!  MEF90HookesLaw2DSum:
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   function MEF90HookesLaw2DSum(A, B)
      type(MEF90HookesLaw2D), intent(IN)           :: A, B
      type(MEF90HookesLaw2D)                       :: MEF90HookesLaw2DSum

      character(len=MEF90MXSTRLEN)                 :: IOBuffer
      PetscErrorCode                               :: ierr

      if ((A%type == MEF90HookesLawTypeIsotropic) .and. (B%type == MEF90HookesLawTypeIsotropic)) then
         MEF90HookesLaw2DSum%type = MEF90HookesLawTypeIsotropic
         if (A%isPlaneStress .eqv. B%isPlaneStress) then
            MEF90HookesLaw2DSum%lambda = A%lambda + B%lambda
            MEF90HookesLaw2DSum%mu = A%mu + B%mu
            MEF90HookesLaw2DSum%isPlaneStress = A%isPlaneStress
         else
            write (IOBuffer, *) "Incompatible planar Hooke law type in "//__FUNCT__//'\n'
            PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
         end if
         if (A%isPlaneStress) then
            MEF90HookesLaw2DSum%PoissonRatio = MEF90HookesLaw2DSum%lambda / (MEF90HookesLaw2DSum%lambda + MEF90HookesLaw2DSum%mu) * 0.5_kr
            MEF90HookesLaw2DSum%YoungsModulus = 2.0_kr * MEF90HookesLaw2DSum%mu * (1.0_kr + MEF90HookesLaw2DSum%PoissonRatio)
            MEF90HookesLaw2DSum%BulkModulus = MEF90HookesLaw2DSum%lambda + MEF90HookesLaw2DSum%mu
         else
            MEF90HookesLaw2DSum%PoissonRatio = MEF90HookesLaw2DSum%lambda / (MEF90HookesLaw2DSum%lambda + 2.0_kr * MEF90HookesLaw2DSum%mu) * 0.5_kr
            MEF90HookesLaw2DSum%YoungsModulus = 2.0_kr * MEF90HookesLaw2DSum%mu * (1.0_kr + MEF90HookesLaw2DSum%PoissonRatio)
            MEF90HookesLaw2DSum%BulkModulus = MEF90HookesLaw2DSum%lambda + MEF90HookesLaw2DSum%mu
         end if
      else if ((A%type == MEF90HookesLawTypeFull) .and. (B%type == MEF90HookesLawTypeFull)) then
         MEF90HookesLaw2DSum%type = MEF90HookesLawTypeFull
         MEF90HookesLaw2DSum%fullTensor = A%fullTensor + B%fullTensor
      else
         write (IOBuffer, *) "Incompatible planar Hooke law type in "//__FUNCT__//'\n'
         PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
      end if
   end function MEF90HookesLaw2DSum

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw3DSum"
!!!
!!!
!!!  MEF90HookesLaw3DSum:
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   function MEF90HookesLaw3DSum(A, B)
      type(MEF90HookesLaw3D), intent(IN)           :: A, B
      type(MEF90HookesLaw3D)                       :: MEF90HookesLaw3DSum

      character(len=MEF90MXSTRLEN)                 :: IOBuffer
      PetscErrorCode                               :: ierr

      if ((A%type == MEF90HookesLawTypeIsotropic) .and. (B%type == MEF90HookesLawTypeIsotropic)) then
         MEF90HookesLaw3DSum%type = MEF90HookesLawTypeIsotropic
         MEF90HookesLaw3DSum%lambda = A%lambda + B%lambda
         MEF90HookesLaw3DSum%mu = A%mu + B%mu
         MEF90HookesLaw3DSum%PoissonRatio = MEF90HookesLaw3DSum%lambda / (MEF90HookesLaw3DSum%lambda + MEF90HookesLaw3DSum%mu) * 0.5_kr
         MEF90HookesLaw3DSum%YoungsModulus = MEF90HookesLaw3DSum%mu * (3.0_kr * MEF90HookesLaw3DSum%lambda + 2.0_kr * MEF90HookesLaw3DSum%mu) / (MEF90HookesLaw3DSum%lambda + MEF90HookesLaw3DSum%mu)
         MEF90HookesLaw3DSum%BulkModulus = MEF90HookesLaw3DSum%lambda + MEF90HookesLaw3DSum%mu * 2.0_kr / 3.0_kr
      else if ((A%type == MEF90HookesLawTypeFull) .and. (B%type == MEF90HookesLawTypeFull)) then
         MEF90HookesLaw3DSum%type = MEF90HookesLawTypeFull
         MEF90HookesLaw3DSum%fullTensor = A%fullTensor + B%fullTensor
      else
         write (IOBuffer, *) "Incompatible planar Hooke law type in "//__FUNCT__//'\n'
         PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
      end if
   end function MEF90HookesLaw3DSum

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw2DDiff"
!!!
!!!
!!!  MEF90HookesLaw2DDiff:
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   function MEF90HookesLaw2DDiff(A, B)
      type(MEF90HookesLaw2D), intent(IN)           :: A, B
      type(MEF90HookesLaw2D)                       :: MEF90HookesLaw2DDiff

      character(len=MEF90MXSTRLEN)                 :: IOBuffer
      PetscErrorCode                               :: ierr

      if ((A%type == MEF90HookesLawTypeIsotropic) .and. (B%type == MEF90HookesLawTypeIsotropic)) then
         MEF90HookesLaw2DDiff%type = MEF90HookesLawTypeIsotropic
         if (A%isPlaneStress .eqv. B%isPlaneStress) then
            MEF90HookesLaw2DDiff%lambda = A%lambda - B%lambda
            MEF90HookesLaw2DDiff%mu = A%mu - B%mu
            MEF90HookesLaw2DDiff%isPlaneStress = A%isPlaneStress
         else
            write (IOBuffer, *) "Incompatible planar Hooke law type in "//__FUNCT__//'\n'
            PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
            SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
         end if
         if (A%isPlaneStress) then
            MEF90HookesLaw2DDiff%PoissonRatio = MEF90HookesLaw2DDiff%lambda / (MEF90HookesLaw2DDiff%lambda + MEF90HookesLaw2DDiff%mu) * 0.5_kr
            MEF90HookesLaw2DDiff%YoungsModulus = 2.0_kr * MEF90HookesLaw2DDiff%mu * (1.0_kr + MEF90HookesLaw2DDiff%PoissonRatio)
            MEF90HookesLaw2DDiff%BulkModulus = MEF90HookesLaw2DDiff%lambda + MEF90HookesLaw2DDiff%mu
         else
            MEF90HookesLaw2DDiff%PoissonRatio = MEF90HookesLaw2DDiff%lambda / (MEF90HookesLaw2DDiff%lambda + 2.0_kr * MEF90HookesLaw2DDiff%mu) * 0.5_kr
            MEF90HookesLaw2DDiff%YoungsModulus = 2.0_kr * MEF90HookesLaw2DDiff%mu * (1.0_kr + MEF90HookesLaw2DDiff%PoissonRatio)
            MEF90HookesLaw2DDiff%BulkModulus = MEF90HookesLaw2DDiff%lambda + MEF90HookesLaw2DDiff%mu
         end if
      else if ((A%type == MEF90HookesLawTypeFull) .and. (B%type == MEF90HookesLawTypeFull)) then
         MEF90HookesLaw2DDiff%type = MEF90HookesLawTypeFull
         MEF90HookesLaw2DDiff%fullTensor = A%fullTensor - B%fullTensor
      else
         write (IOBuffer, *) "Incompatible planar Hooke law type in "//__FUNCT__//'\n'
         PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
      end if
   end function MEF90HookesLaw2DDiff

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw3DDiff"
!!!
!!!
!!!  MEF90HookesLaw3DDiff:
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   function MEF90HookesLaw3DDiff(A, B)
      type(MEF90HookesLaw3D), intent(IN)           :: A, B
      type(MEF90HookesLaw3D)                       :: MEF90HookesLaw3DDiff

      character(len=MEF90MXSTRLEN)                 :: IOBuffer
      PetscErrorCode                               :: ierr

      if ((A%type == MEF90HookesLawTypeIsotropic) .and. (B%type == MEF90HookesLawTypeIsotropic)) then
         MEF90HookesLaw3DDiff%type = MEF90HookesLawTypeIsotropic
         MEF90HookesLaw3DDiff%lambda = A%lambda - B%lambda
         MEF90HookesLaw3DDiff%mu = A%mu - B%mu
         MEF90HookesLaw3DDiff%PoissonRatio = MEF90HookesLaw3DDiff%lambda / (MEF90HookesLaw3DDiff%lambda + MEF90HookesLaw3DDiff%mu) * 0.5_kr
         MEF90HookesLaw3DDiff%YoungsModulus = MEF90HookesLaw3DDiff%mu * (3.0_kr * MEF90HookesLaw3DDiff%lambda + 2.0_kr * MEF90HookesLaw3DDiff%mu) / (MEF90HookesLaw3DDiff%lambda + MEF90HookesLaw3DDiff%mu)
         MEF90HookesLaw3DDiff%BulkModulus = MEF90HookesLaw3DDiff%lambda + MEF90HookesLaw3DDiff%mu * 2.0_kr / 3.0_kr
      else if ((A%type == MEF90HookesLawTypeFull) .and. (B%type == MEF90HookesLawTypeFull)) then
         MEF90HookesLaw3DDiff%type = MEF90HookesLawTypeFull
         MEF90HookesLaw3DDiff%fullTensor = A%fullTensor - B%fullTensor
      else
         write (IOBuffer, *) "Incompatible planar Hooke law type in "//__FUNCT__//'\n'
         PetscCall(PetscPrintf(PETSC_COMM_SELF, IOBuffer, ierr))
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_SUP, IOBuffer)
      end if
   end function MEF90HookesLaw3DDiff

#undef __FUNCT__
#define __FUNCT__ "ScalarXMEF90HookesLaw2D"
!!!
!!!
!!!  ScalarXMEF90HookesLaw2D:
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   function ScalarXMEF90HookesLaw2D(t, A)
      PetscReal, intent(IN)                        :: t
      type(MEF90HookesLaw2D), intent(IN)           :: A
      type(MEF90HookesLaw2D)                       :: ScalarXMEF90HookesLaw2D

      if (A%type == MEF90HookesLawTypeIsotropic) then
         ScalarXMEF90HookesLaw2D%type = MEF90HookesLawTypeIsotropic
         ScalarXMEF90HookesLaw2D%lambda = t * A%lambda
         ScalarXMEF90HookesLaw2D%mu = t * A%mu
         ScalarXMEF90HookesLaw2D%isPlaneStress = A%isPlaneStress
         if (A%isPlaneStress) then
            ScalarXMEF90HookesLaw2D%PoissonRatio = ScalarXMEF90HookesLaw2D%lambda / (ScalarXMEF90HookesLaw2D%lambda + ScalarXMEF90HookesLaw2D%mu) * 0.5_kr
            ScalarXMEF90HookesLaw2D%YoungsModulus = 2.0_kr * ScalarXMEF90HookesLaw2D%mu * (1.0_kr + ScalarXMEF90HookesLaw2D%PoissonRatio)
            ScalarXMEF90HookesLaw2D%BulkModulus = ScalarXMEF90HookesLaw2D%lambda + ScalarXMEF90HookesLaw2D%mu
         else
            ScalarXMEF90HookesLaw2D%PoissonRatio = ScalarXMEF90HookesLaw2D%lambda / (ScalarXMEF90HookesLaw2D%lambda + 2.0_kr * ScalarXMEF90HookesLaw2D%mu) * 0.5_kr
            ScalarXMEF90HookesLaw2D%YoungsModulus = 2.0_kr * ScalarXMEF90HookesLaw2D%mu * (1.0_kr + ScalarXMEF90HookesLaw2D%PoissonRatio)
            ScalarXMEF90HookesLaw2D%BulkModulus = ScalarXMEF90HookesLaw2D%lambda + ScalarXMEF90HookesLaw2D%mu
         end if
      else
         ScalarXMEF90HookesLaw2D%type = MEF90HookesLawTypeFull
         ScalarXMEF90HookesLaw2D%fullTensor = t * A%fullTensor
      end if
   end function ScalarXMEF90HookesLaw2D

#undef __FUNCT__
#define __FUNCT__ "ScalarXMEF90HookesLaw3D"
!!!
!!!
!!!  ScalarXMEF90HookesLaw3D:
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!

   function ScalarXMEF90HookesLaw3D(t, A)
      PetscReal, intent(IN)                        :: t
      type(MEF90HookesLaw3D), intent(IN)           :: A
      type(MEF90HookesLaw3D)                       :: ScalarXMEF90HookesLaw3D

      if (A%type == MEF90HookesLawTypeIsotropic) then
         ScalarXMEF90HookesLaw3D%type = MEF90HookesLawTypeIsotropic
         ScalarXMEF90HookesLaw3D%lambda = t * A%lambda
         ScalarXMEF90HookesLaw3D%mu = t * A%mu
         ScalarXMEF90HookesLaw3D%PoissonRatio = ScalarXMEF90HookesLaw3D%lambda / (ScalarXMEF90HookesLaw3D%lambda + ScalarXMEF90HookesLaw3D%mu) * 0.5_kr
         ScalarXMEF90HookesLaw3D%YoungsModulus = ScalarXMEF90HookesLaw3D%mu * (3.0_kr * ScalarXMEF90HookesLaw3D%lambda + 2.0_kr * ScalarXMEF90HookesLaw3D%mu) / (ScalarXMEF90HookesLaw3D%lambda + ScalarXMEF90HookesLaw3D%mu)
         ScalarXMEF90HookesLaw3D%BulkModulus = ScalarXMEF90HookesLaw3D%lambda + ScalarXMEF90HookesLaw3D%mu * 2.0_kr / 3.0_kr
      else
         ScalarXMEF90HookesLaw3D%type = MEF90HookesLawTypeFull
         ScalarXMEF90HookesLaw3D%fullTensor = t * A%fullTensor
      end if
   end function ScalarXMEF90HookesLaw3D

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw2DXMatS2D"
!!!
!!!
!!!  MEF90HookesLaw2DXMatS2D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!

   function MEF90HookesLaw2DXMatS2D(A, X)
      type(MEF90HookesLaw2D), intent(IN)           :: A
      type(MatS2D), intent(IN)                     :: X
      type(MatS2D)                                 :: MEF90HookesLaw2DXMatS2D

      real(Kind=Kr)                                :: C1, C2

      select case (A%type)
      case (MEF90HookesLawTypeIsotropic)
         C1 = A%lambda + 2.0_kr * A%mu
         C2 = 2.0_kr * A%mu
         MEF90HookesLaw2DXMatS2D%XX = C1 * X%XX + A%lambda * X%YY
         MEF90HookesLaw2DXMatS2D%YY = A%lambda * X%XX + C1 * X%YY
         MEF90HookesLaw2DXMatS2D%XY = C2 * X%XY
      case (MEF90HookesLawTypeFull)
         MEF90HookesLaw2DXMatS2D = A%fullTensor * X
      end select
   end function MEF90HookesLaw2DXMatS2D

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw2DXMatS3D"
!!!
!!!
!!!  MEF90HookesLaw2DXMatS3D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!

   function MEF90HookesLaw3DXMatS3D(A, X)
      type(MEF90HookesLaw3D), intent(IN)           :: A
      type(MatS3D), intent(IN)                     :: X
      type(MatS3D)                                 :: MEF90HookesLaw3DXMatS3D

      real(Kind=Kr)                                :: C1, C2

      select case (A%type)
      case (MEF90HookesLawTypeIsotropic)
         C1 = A%lambda + 2.0_kr * A%mu
         C2 = 2.0_kr * A%mu
         MEF90HookesLaw3DXMatS3D%XX = C1 * X%XX + A%lambda * X%YY + A%lambda * X%ZZ
         MEF90HookesLaw3DXMatS3D%YY = A%lambda * X%XX + C1 * X%YY + A%lambda * X%ZZ
         MEF90HookesLaw3DXMatS3D%ZZ = A%lambda * X%XX + A%lambda * X%YY + C1 * X%ZZ
         MEF90HookesLaw3DXMatS3D%YZ = C2 * X%YZ
         MEF90HookesLaw3DXMatS3D%XZ = C2 * X%XZ
         MEF90HookesLaw3DXMatS3D%XY = C2 * X%XY

      case (MEF90HookesLawTypeFull)
         MEF90HookesLaw3DXMatS3D = A%fullTensor * X
      end select
   end function MEF90HookesLaw3DXMatS3D

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw2DXMat2D"
!!!
!!!
!!!  MEF90HookesLaw2DXMat2D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!
   function MEF90HookesLaw2DXMat2D(A, X)
      type(MEF90HookesLaw2D), intent(IN)           :: A
      type(Mat2D), intent(IN)                      :: X
      type(Mat2D)                                  :: MEF90HookesLaw2DXMat2D

      real(Kind=Kr)                                :: C1, C2

      select case (A%type)
      case (MEF90HookesLawTypeIsotropic)
         C1 = A%lambda + 2.0_kr * A%mu
         C2 = 2.0_kr * A%mu
         MEF90HookesLaw2DXMat2D%XX = C1 * X%XX + A%lambda * X%YY
         MEF90HookesLaw2DXMat2D%XY = C2 * X%XY
         MEF90HookesLaw2DXMat2D%YY = A%lambda * X%XX + C1 * X%YY
         MEF90HookesLaw2DXMat2D%YX = C2 * X%YX
      case (MEF90HookesLawTypeFull)
         MEF90HookesLaw2DXMat2D = A%fullTensor * X
      end select
   end function MEF90HookesLaw2DXMat2D

#undef __FUNCT__
#define __FUNCT__ "MEF90HookesLaw2DXMat3D"
!!!
!!!
!!!  MEF90HookesLaw2DXMat3D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!
   function MEF90HookesLaw3DXMat3D(A, X)
      type(MEF90HookesLaw3D), intent(IN)           :: A
      type(Mat3D), intent(IN)                      :: X
      type(Mat3D)                                  :: MEF90HookesLaw3DXMat3D

      real(Kind=Kr)                                :: C1, C2

      select case (A%type)
      case (MEF90HookesLawTypeIsotropic)
         C1 = A%lambda + 2.0_kr * A%mu
         C2 = 2.0_kr * A%mu
         MEF90HookesLaw3DXMat3D%XX = C1 * X%XX + A%lambda * X%YY + A%lambda * X%ZZ
         MEF90HookesLaw3DXMat3D%XY = C2 * X%XY
         MEF90HookesLaw3DXMat3D%XZ = C2 * X%XZ

         MEF90HookesLaw3DXMat3D%YX = C2 * X%YX
         MEF90HookesLaw3DXMat3D%YY = A%lambda * X%XX + C1 * X%YY + A%lambda * X%ZZ
         MEF90HookesLaw3DXMat3D%YZ = C2 * X%YZ

         MEF90HookesLaw3DXMat3D%XZ = C2 * X%XZ
         MEF90HookesLaw3DXMat3D%YZ = C2 * X%YZ
         MEF90HookesLaw3DXMat3D%ZZ = A%lambda * X%XX + A%lambda * X%YY + C1 * X%ZZ
      case (MEF90HookesLawTypeFull)
         MEF90HookesLaw3DXMat3D = A%fullTensor * X
      end select
   end function MEF90HookesLaw3DXMat3D
end module m_MEF90_Materials
