module m_MEF90_LinAlg_class
   implicit none(type, external)
   private
   public :: mef90Vect
   public :: mef90Mat
   public :: mef90Tens4OS

!!!
!!!
!!!  m_MEF90_LinAlg_class: Starting a OO implementation of the basic classes in m_MEF90_LinAlg
!!!
!!!  (c) 2020 Blaise Bourdin bourdin@lsu.edu
!!!      2022 Blaise Bourdin bourdin@mcmaster.ca
!!!

   type :: mef90Vect
   end type mef90Vect

   type :: mef90Mat
   end type mef90Mat

   type :: mef90Tens4OS
   end type mef90Tens4OS
end module m_MEF90_LinAlg_class

module m_MEF90_LinAlg
#include "petsc/finclude/petsc.h"
   use petscsys
   use m_MEF90_LinAlg_class, only: mef90Vect, mef90Mat, mef90Tens4OS
   use m_MEF90_Parameters

   implicit none(type)
   private
   public :: mef90Vect, mef90Mat, mef90Tens4OS
   public :: Vect2D, Vect3D
   public :: Mat2D, MatS2D, Mat3D, MatS3D
   public :: Tens4OS2D, Tens4OS3D
   public :: MEF90Vect2De1, MEF90Vect2De2
   public :: MEF90Vect3De1, MEF90Vect3De2, MEF90Vect3De3
   public :: MEF90Mat2DIdentity, MEF90MatS2DIdentity
   public :: MEF90Mat3DIdentity, MEF90MatS3DIdentity
   public :: MEF90Tens4OS2DIdentity, MEF90Tens4OS3DIdentity
   public :: operator(+), operator(-), operator(*), operator(/)
   public :: operator(.DotP.), operator(.CrossP.)
   public :: operator(.TensP.), operator(.SymP.), operator(.oDot.)
   public :: Transpose, Invert, Trace, Det
   public :: assignment(=)
   public :: Symmetrize, DeviatoricPart, HydrostaticPart
   public :: MEF90MatRaRt, MEF90MatRtaR
   public :: Norm, simplexNormal, sqrt
   public :: SpectralDecomposition, Diagonalize, Moment, Tens4OSTransform

   external :: dsyevd
   external :: dgetrf
   external :: dgetri

   type, extends(mef90Vect) :: Vect2D
      PetscReal          :: X = 0.0_Kr
      PetscReal          :: Y = 0.0_Kr
   end type Vect2D

   type, extends(mef90Vect) :: Vect3D
      PetscReal          :: X = 0.0_Kr
      PetscReal          :: Y = 0.0_Kr
      PetscReal          :: Z = 0.0_Kr
   end type Vect3D

   type, extends(mef90Mat) :: Mat2D
      PetscReal          :: XX = 0.0_Kr
      PetscReal          :: XY = 0.0_Kr
      PetscReal          :: YX = 0.0_Kr
      PetscReal          :: YY = 0.0_Kr
   end type Mat2D

   type, extends(mef90Mat) :: MatS2D
      PetscReal          :: XX = 0.0_Kr
      PetscReal          :: YY = 0.0_Kr
      PetscReal          :: XY = 0.0_Kr
   end type MatS2D

   type, extends(mef90Mat) :: Mat3D
      PetscReal          :: XX = 0.0_Kr
      PetscReal          :: XY = 0.0_Kr
      PetscReal          :: XZ = 0.0_Kr
      PetscReal          :: YX = 0.0_Kr
      PetscReal          :: YY = 0.0_Kr
      PetscReal          :: YZ = 0.0_Kr
      PetscReal          :: ZX = 0.0_Kr
      PetscReal          :: ZY = 0.0_Kr
      PetscReal          :: ZZ = 0.0_Kr
   end type Mat3D

   type, extends(mef90Mat) :: MatS3D
      PetscReal          :: XX = 0.0_Kr
      PetscReal          :: YY = 0.0_Kr
      PetscReal          :: ZZ = 0.0_Kr
      PetscReal          :: YZ = 0.0_Kr
      PetscReal          :: XZ = 0.0_Kr
      PetscReal          :: XY = 0.0_Kr
   end type MatS3D

 !! After much hesitation,
 !! - the terms are numbered in alphabetical order (i.e. XXYX and not XYXX)
 !! - the terms are stored in alphabetical order
 !! 2014-07: Changed ordering to rows of the upper triangular part
 !!          and naming to be consistent with Voigt notations
   type, extends(mef90Tens4OS) :: Tens4OS2D
      PetscReal          :: XXXX, XXYY, XXXY
      PetscReal          ::       YYYY, YYXY
      PetscReal          ::             XYXY

   end type Tens4OS2D

   type, extends(mef90Tens4OS) :: Tens4OS3D
      PetscReal          :: XXXX, XXYY, XXZZ, XXYZ, XXXZ, XXXY
      PetscReal          ::       YYYY, YYZZ, YYYZ, YYXZ, YYXY
      PetscReal          ::             ZZZZ, ZZYZ, ZZXZ, ZZXY
      PetscReal          ::                   YZYZ, YZXZ, YZXY
      PetscReal          ::                         XZXZ, XZXY
      PetscReal          ::                               XYXY
   end type Tens4OS3D

   type(Vect2D), parameter    :: MEF90Vect2De1 = Vect2D(1.0_kr, 0.0_kr)
   type(Vect2D), parameter    :: MEF90Vect2De2 = Vect2D(0.0_kr, 1.0_kr)

   type(Vect3D), parameter    :: MEF90Vect3De1 = Vect3D(1.0_kr, 0.0_kr, 0.0_kr)
   type(Vect3D), parameter    :: MEF90Vect3De2 = Vect3D(0.0_kr, 1.0_kr, 0.0_kr)
   type(Vect3D), parameter    :: MEF90Vect3De3 = Vect3D(0.0_kr, 0.0_kr, 1.0_kr)

   type(Mat2D), parameter     :: MEF90Mat2DIdentity = Mat2D(1.0_kr, 0.0_kr, &
                                                            0.0_kr, 1.0_kr)
   type(MatS2D), parameter    :: MEF90MatS2DIdentity = MatS2D(1.0_kr, 1.0_kr, 0.0_kr)
   type(Mat3D), parameter     :: MEF90Mat3DIdentity = Mat3D(1.0_kr, 0.0_kr, 0.0_kr, &
                                                            0.0_kr, 1.0_kr, 0.0_kr, &
                                                            0.0_kr, 0.0_kr, 1.0_kr)
   type(MatS3D), parameter    :: MEF90MatS3DIdentity = MatS3D(1.0_kr, 1.0_kr, 1.0_kr, &
                                                              0.0_kr, 0.0_kr, 0.0_kr)

   type(Tens4OS2D), parameter :: MEF90Tens4OS2DIdentity = Tens4OS2D(1.0_kr, 0.0_kr, 0.0_kr, &
                                                                    1.0_kr, 0.0_kr, &
                                                                    0.5_kr)
   type(Tens4OS3D), parameter :: MEF90Tens4OS3DIdentity = Tens4OS3D(1.0_kr, 0.0_kr, 0.0_kr, 0.0_kr, 0.0_kr, 0.0_kr, &
                                                                    1.0_kr, 0.0_kr, 0.0_kr, 0.0_kr, 0.0_kr, &
                                                                    1.0_kr, 0.0_kr, 0.0_kr, 0.0_kr, &
                                                                    0.5_kr, 0.0_kr, 0.0_kr, &
                                                                    0.5_kr, 0.0_kr, &
                                                                    0.5_kr)

   interface operator(+)
      module procedure SumVect2D, SumVect3D, SumMat2D, SumMat3D, SumMatS2D, SumMatS3D, SumTens4OS2D, SumTens4OS3D
   end interface

   interface operator(-)
      module procedure DifVect2D, DifVect3D, DifMat2D, DifMat3D, DifMatS2D, DifMatS3D, DifTens4OS2D, DifTens4OS3D
   end interface

   interface operator(*)
      module procedure DbleXVect2D, Vect2DXDble, DbleXVect3D, Vect3DXDble, &
         DbleXMat2D, Mat2DXDble, DbleXMat3D, Mat3DXDble, &
         DbleXMatS2D, MatS2DXDble, DbleXMatS3D, MatS3DXDble, &
         MatXVect2D, MatXVect3D, MatXVect2DS, MatXVect3DS, &
         DbleXTens4OS2D, Tens4OS2DXDble, Tens4OS2DXMatS2D, Tens4OS2DXMat2D, &
         DbleXTens4OS3D, Tens4OS3DXDble, Tens4OS3DXMatS3D, Tens4OS3DXMat3D, &
         Mat2DXMat2D, Mat2DXMatS2D, MatS2DXMatS2D, Mat3DXMat3D, Mat3DXMatS3D, MatS3DXMatS3D, &
         DotP2D, DotP3D
   end interface

   interface operator(/)
      module procedure Vect2DQuot, Vect3DQuot, Mat2DQuot, Mat3DQuot, MatS2DQuot, MatS3DQuot, Tens4OS2DQuot, Tens4OS3DQuot
   end interface

   interface operator(.DotP.)
      module procedure DotP2D, DotP3D, ContP2D, ContP3D, ContP2DS, ContP3DS, &
         Mat2DDotMatS2D, MatS2DDotMat2D, Mat3DDotMatS3D, MatS3DDotMat3D
   end interface

   interface operator(.CrossP.)
      module procedure CrossP3D
   end interface

   interface Transpose
      module procedure Transpose2D, Transpose3D
   end interface

   interface Invert
      module procedure InvertMat2D, InvertMatS2D, InvertMat3D, InvertMatS3D, InvertTens4OS2D, InvertTens4OS3D
   end interface

   interface operator(.TensP.)
      module procedure TensPVect2D, TensPVect3D
   end interface

   interface operator(.SymP.)
      module procedure SymPVect2D, SymPVect3D, SymPMatS2D, SymPMatS3D
   end interface

   interface operator(.oDot.)
      module procedure oDotMatS2D, oDotMatS3D
   end interface

   interface Trace
      module procedure Trace2D, Trace3D, Trace2DS, Trace3DS
   end interface

   interface Det
      module procedure DetMat2D, DetMatS2D, DetMat3D, DetMatS3D
   end interface

   interface assignment(=)
      module procedure Vect2D_Get_Real, Vect3D_Get_Real, &
         Vect2D_Get_VectR, Vect3D_Get_VectR, &
         Vect2DEQ, Vect3DEQ, Mat2D_Get_Real, Mat3D_Get_Real, &
         Mat2DEQ, Mat3DEQ, MatS2D_Get_Real, MatS3D_Get_Real, &
         MatS2DEQ, MatS3DEQ, MatS2D_Get_VectR, MatS3D_Get_VectR, &
         VectR_Get_MatS2D, VectR_Get_MatS3D, &
         Tens4OS2D_Get_Real, Tens4OS3D_Get_Real, &
         Tens4OS2DToArray, ArrayToTens4OS2D, &
         Tens4OS3DToArray, ArrayToTens4OS3D, &
         Mat2DGetArray, Mat3DGetArray, &
         MatS2DGetArray, MatS3DGetArray, &
         ArrayGetMat2D, ArrayGetMat3D, &
         ArrayGetMatS2D, ArrayGetMatS3D, &
         MatS3DToMat3D, MatS2DToMat2D, &
         Mat3DToMatS3D, Mat2DToMatS2D, &
         Tens4OS2D2Array4, Tens4OS3D2Array4, &
         Array42Tens4OS2D, Array42Tens4OS3D
   end interface

   interface Symmetrize
      module procedure Symmetrize2D, Symmetrize3D
   end interface

   interface DeviatoricPart
      module procedure DeviatoricPart2D, DeviatoricPart2DS, DeviatoricPart3D, DeviatoricPart3DS
   end interface

   interface HydrostaticPart
      module procedure HydrostaticPart2D, HydrostaticPart2DS, HydrostaticPart3D, HydrostaticPart3DS
   end interface

   interface MEF90MatRaRt
      module procedure RaRtMat2D, RaRtMatS2D, RaRtMat3D, RaRtMatS3D
   end interface

   interface MEF90MatRtaR
      module procedure RtaRMat2D, RtaRMatS2D, RtaRMat3D, RtaRMatS3D
   end interface

   interface Norm
      module procedure Vect2DNorm, Vect3DNorm, Mat2DNorm, MatS2DNorm, Mat3DNorm, MatS3DNorm
      ! Tens4OS2DNorm,Tens4OS3DNorm
   end interface

   interface simplexNormal
      module procedure simplexNormal2D, simplexNormal3D
   end interface simplexNormal

   interface sqrt
      module procedure Tens4OS2DSquareRoot, Tens4OS3DSquareRoot
   end interface

   interface SpectralDecomposition
      module procedure MatS2DSpectralDecomposition, MatS3DSpectralDecomposition
   end interface

   interface Diagonalize
      !!! Diagonalize(A,P,D) returns P,D such that A = P D P^{-1}
      !!! The diagonal entries of D are sorted in increasing order.
      module procedure MatS2DEigenVectorValues, MatS3DEigenVectorValues
   end interface

   interface Moment
      module procedure Mat2DMoment, MatS2DMoment, Mat3DMoment, MatS3DMoment
   end interface

   interface Tens4OSTransform
      module procedure Tens4OS2DTransform, Tens4OS3DTransform
   end interface

   !Interface MatSymToMat
   !   Module Procedure MatS3DToMat3D,MatS2DToMat2D
   !End Interface

!!$  Type(Vect2D),Parameter       :: e1_2D = [1.0_Kr,0.0_Kr]
!!$  Type(Vect2D),Parameter       :: e2_2D = [0.0_Kr,1.0_Kr]
!!$
!!$  Type(Vect3D),Parameter       :: e1_3D = [1.0_Kr,0.0_Kr,0.0_Kr]
!!$  Type(Vect3D),Parameter       :: e2_3D = [0.0_Kr,1.0_Kr,0.0_Kr]
!!$  Type(Vect3D),Parameter       :: e3_3D = [0.0_Kr,0.0_Kr,1.0_Kr]

contains
   function SumVect2D(V1, V2)
      type(Vect2D), intent(IN)                     :: V1
      type(Vect2D), intent(IN)                     :: V2
      type(Vect2D)                                :: SumVect2D

      SumVect2D%X = V1%X + V2%X
      SumVect2D%Y = V1%Y + V2%Y
   end function SumVect2D

   function SumVect3D(V1, V2)
      type(Vect3D), intent(IN)                     :: V1, V2
      type(Vect3D)                                :: SumVect3D

      SumVect3D%X = V1%X + V2%X
      SumVect3D%Y = V1%Y + V2%Y
      SumVect3D%Z = V1%Z + V2%Z
   end function SumVect3D

   function SumMat2D(M1, M2)
      type(Mat2D), intent(IN)                      :: M1, M2
      type(Mat2D)                                 :: SumMat2D

      SumMat2D%XX = M1%XX + M2%XX
      SumMat2D%XY = M1%XY + M2%XY
      SumMat2D%YX = M1%YX + M2%YX
      SumMat2D%YY = M1%YY + M2%YY
   end function SumMat2D

   function SumMatS2D(M1, M2)
      type(MatS2D), intent(IN)                     :: M1, M2
      type(MatS2D)                                :: SumMatS2D

      SumMatS2D%XX = M1%XX + M2%XX
      SumMatS2D%XY = M1%XY + M2%XY
      SumMatS2D%YY = M1%YY + M2%YY
   end function SumMatS2D

   function SumMat3D(M1, M2)
      type(Mat3D), intent(IN)                      :: M1, M2
      type(Mat3D)                                 :: SumMat3D

      SumMat3D%XX = M1%XX + M2%XX
      SumMat3D%XY = M1%XY + M2%XY
      SumMat3D%XZ = M1%XZ + M2%XZ
      SumMat3D%YX = M1%YX + M2%YX
      SumMat3D%YY = M1%YY + M2%YY
      SumMat3D%YZ = M1%YZ + M2%YZ
      SumMat3D%ZX = M1%ZX + M2%ZX
      SumMat3D%ZY = M1%ZY + M2%ZY
      SumMat3D%ZZ = M1%ZZ + M2%ZZ
   end function SumMat3D

   function SumMatS3D(M1, M2)
      type(MatS3D), intent(IN)                     :: M1, M2
      type(MatS3D)                                :: SumMatS3D

      SumMatS3D%XX = M1%XX + M2%XX
      SumMatS3D%YY = M1%YY + M2%YY
      SumMatS3D%ZZ = M1%ZZ + M2%ZZ
      SumMatS3D%YZ = M1%YZ + M2%YZ
      SumMatS3D%XZ = M1%XZ + M2%XZ
      SumMatS3D%XY = M1%XY + M2%XY
   end function SumMatS3D

   function SumTens4OS2D(T1, T2)
      type(Tens4OS2D), intent(IN)                  :: T1, T2
      type(Tens4OS2D)                             :: SumTens4OS2D

      SumTens4OS2D%XXXX = T1%XXXX + T2%XXXX
      SumTens4OS2D%XXYY = T1%XXYY + T2%XXYY
      SumTens4OS2D%XXXY = T1%XXXY + T2%XXXY

      SumTens4OS2D%YYYY = T1%YYYY + T2%YYYY
      SumTens4OS2D%YYXY = T1%YYXY + T2%YYXY

      SumTens4OS2D%XYXY = T1%XYXY + T2%XYXY

   end function SumTens4OS2D

   function SumTens4OS3D(T1, T2)
      type(Tens4OS3D), intent(IN)                  :: T1, T2
      type(Tens4OS3D)                             :: SumTens4OS3D

      SumTens4OS3D%XXXX = T1%XXXX + T2%XXXX
      SumTens4OS3D%XXYY = T1%XXYY + T2%XXYY
      SumTens4OS3D%XXZZ = T1%XXZZ + T2%XXZZ
      SumTens4OS3D%XXYZ = T1%XXYZ + T2%XXYZ
      SumTens4OS3D%XXXZ = T1%XXXZ + T2%XXXZ
      SumTens4OS3D%XXXY = T1%XXXY + T2%XXXY

      SumTens4OS3D%YYYY = T1%YYYY + T2%YYYY
      SumTens4OS3D%YYZZ = T1%YYZZ + T2%YYZZ
      SumTens4OS3D%YYYZ = T1%YYYZ + T2%YYYZ
      SumTens4OS3D%YYXZ = T1%YYXZ + T2%YYXZ
      SumTens4OS3D%YYXY = T1%YYXY + T2%YYXY

      SumTens4OS3D%ZZZZ = T1%ZZZZ + T2%ZZZZ
      SumTens4OS3D%ZZYZ = T1%ZZYZ + T2%ZZYZ
      SumTens4OS3D%ZZXZ = T1%ZZXZ + T2%ZZXZ
      SumTens4OS3D%ZZXY = T1%ZZXY + T2%ZZXY

      SumTens4OS3D%YZYZ = T1%YZYZ + T2%YZYZ
      SumTens4OS3D%YZXZ = T1%YZXZ + T2%YZXZ
      SumTens4OS3D%YZXY = T1%YZXY + T2%YZXY

      SumTens4OS3D%XZXZ = T1%XZXZ + T2%XZXZ
      SumTens4OS3D%XZXY = T1%XZXY + T2%XZXY

      SumTens4OS3D%XYXY = T1%XYXY + T2%XYXY
   end function SumTens4OS3D

   ! Overloading "-"
   function DifVect2D(V1, V2)
      type(Vect2D), intent(IN)                     :: V1
      type(Vect2D), intent(IN)                     :: V2
      type(Vect2D)                                 :: DifVect2D

      DifVect2D%X = V1%X - V2%X
      DifVect2D%Y = V1%Y - V2%Y
   end function DifVect2D

   function DifVect3D(V1, V2)
      type(Vect3D), intent(IN)                     :: V1, V2
      type(Vect3D)                                 :: DifVect3D

      DifVect3D%X = V1%X - V2%X
      DifVect3D%Y = V1%Y - V2%Y
      DifVect3D%Z = V1%Z - V2%Z
   end function DifVect3D

   function DifMat2D(M1, M2)
      type(Mat2D), intent(IN)                      :: M1, M2
      type(Mat2D)                                  :: DifMat2D

      DifMat2D%XX = M1%XX - M2%XX
      DifMat2D%XY = M1%XY - M2%XY
      DifMat2D%YX = M1%YX - M2%YX
      DifMat2D%YY = M1%YY - M2%YY
   end function DifMat2D

   function DifMatS2D(M1, M2)
      type(MatS2D), intent(IN)                      :: M1, M2
      type(MatS2D)                                  :: DifMatS2D

      DifMatS2D%XX = M1%XX - M2%XX
      DifMatS2D%YY = M1%YY - M2%YY
      DifMatS2D%XY = M1%XY - M2%XY
   end function DifMatS2D

   function DifMat3D(M1, M2)
      type(Mat3D), intent(IN)                      :: M1, M2
      type(Mat3D)                                  :: DifMat3D

      DifMat3D%XX = M1%XX - M2%XX
      DifMat3D%XY = M1%XY - M2%XY
      DifMat3D%XZ = M1%XZ - M2%XZ
      DifMat3D%YX = M1%YX - M2%YX
      DifMat3D%YY = M1%YY - M2%YY
      DifMat3D%YZ = M1%YZ - M2%YZ
      DifMat3D%ZX = M1%ZX - M2%ZX
      DifMat3D%ZY = M1%ZY - M2%ZY
      DifMat3D%ZZ = M1%ZZ - M2%ZZ
   end function DifMat3D

   function DifMatS3D(M1, M2)
      type(MatS3D), intent(IN)                     :: M1, M2
      type(MatS3D)                                 :: DifMatS3D

      DifMatS3D%XX = M1%XX - M2%XX
      DifMatS3D%YY = M1%YY - M2%YY
      DifMatS3D%ZZ = M1%ZZ - M2%ZZ
      DifMatS3D%YZ = M1%YZ - M2%YZ
      DifMatS3D%XZ = M1%XZ - M2%XZ
      DifMatS3D%XY = M1%XY - M2%XY
   end function DifMatS3D

   function DifTens4OS2D(T1, T2)
      type(Tens4OS2D), intent(IN)                  :: T1, T2
      type(Tens4OS2D)                              :: DifTens4OS2D

      DifTens4OS2D%XXXX = T1%XXXX - T2%XXXX
      DifTens4OS2D%XXYY = T1%XXYY - T2%XXYY
      DifTens4OS2D%XXXY = T1%XXXY - T2%XXXY

      DifTens4OS2D%YYYY = T1%YYYY - T2%YYYY
      DifTens4OS2D%YYXY = T1%YYXY - T2%YYXY

      DifTens4OS2D%XYXY = T1%XYXY - T2%XYXY
   end function DifTens4OS2D

   function DifTens4OS3D(T1, T2)
      type(Tens4OS3D), intent(IN)                  :: T1, T2
      type(Tens4OS3D)                              :: DifTens4OS3D

      DifTens4OS3D%XXXX = T1%XXXX - T2%XXXX
      DifTens4OS3D%XXYY = T1%XXYY - T2%XXYY
      DifTens4OS3D%XXZZ = T1%XXZZ - T2%XXZZ
      DifTens4OS3D%XXYZ = T1%XXYZ - T2%XXYZ
      DifTens4OS3D%XXXZ = T1%XXXZ - T2%XXXZ
      DifTens4OS3D%XXXY = T1%XXXY - T2%XXXY

      DifTens4OS3D%YYYY = T1%YYYY - T2%YYYY
      DifTens4OS3D%YYZZ = T1%YYZZ - T2%YYZZ
      DifTens4OS3D%YYYZ = T1%YYYZ - T2%YYYZ
      DifTens4OS3D%YYXZ = T1%YYXZ - T2%YYXZ
      DifTens4OS3D%YYXY = T1%YYXY - T2%YYXY

      DifTens4OS3D%ZZZZ = T1%ZZZZ - T2%ZZZZ
      DifTens4OS3D%ZZYZ = T1%ZZYZ - T2%ZZYZ
      DifTens4OS3D%ZZXZ = T1%ZZXZ - T2%ZZXZ
      DifTens4OS3D%ZZXY = T1%ZZXY - T2%ZZXY

      DifTens4OS3D%YZYZ = T1%YZYZ - T2%YZYZ
      DifTens4OS3D%YZXZ = T1%YZXZ - T2%YZXZ
      DifTens4OS3D%YZXY = T1%YZXY - T2%YZXY

      DifTens4OS3D%XZXZ = T1%XZXZ - T2%XZXZ
      DifTens4OS3D%XZXY = T1%XZXY - T2%XZXY

      DifTens4OS3D%XYXY = T1%XYXY - T2%XYXY
   end function DifTens4OS3D

   ! Overloading "*"
   function DbleXVect2D(D1, V1)
      PetscReal, intent(IN)                        :: D1
      type(Vect2D), intent(IN)                     :: V1
      type(Vect2D)                                 :: DbleXVect2D

      DbleXVect2D%X = D1 * V1%X
      DbleXVect2D%Y = D1 * V1%Y
   end function DbleXVect2D

   function Vect2DXDble(V1, D1)
      PetscReal, intent(IN)                        :: D1
      type(Vect2D), intent(IN)                     :: V1
      type(Vect2D)                                 :: Vect2DXDble

      Vect2DXDble%X = D1 * V1%X
      Vect2DXDble%Y = D1 * V1%Y
   end function Vect2DXDble

   function DbleXVect3D(D1, V1)
      PetscReal, intent(IN)                        :: D1
      type(Vect3D), intent(IN)                     :: V1
      type(Vect3D)                                 :: DbleXVect3D

      DbleXVect3D%X = D1 * V1%X
      DbleXVect3D%Y = D1 * V1%Y
      DbleXVect3D%Z = D1 * V1%Z
   end function DbleXVect3D

   function Vect3DXDble(V1, D1)
      PetscReal, intent(IN)                        :: D1
      type(Vect3D), intent(IN)                     :: V1
      type(Vect3D)                                 :: Vect3DXDble

      Vect3DXDble%X = D1 * V1%X
      Vect3DXDble%Y = D1 * V1%Y
      Vect3DXDble%Z = D1 * V1%Z
   end function Vect3DXDble

   function DbleXMat2D(D1, M1)
      PetscReal, intent(IN)                        :: D1
      type(Mat2D), intent(IN)                      :: M1
      type(Mat2D)                                  :: DbleXMat2D

      DbleXMat2D%XX = D1 * M1%XX
      DbleXMat2D%XY = D1 * M1%XY
      DbleXMat2D%YX = D1 * M1%YX
      DbleXMat2D%YY = D1 * M1%YY
   end function DbleXMat2D

   function Mat2DXDble(M1, D1)
      PetscReal, intent(IN)                        :: D1
      type(Mat2D), intent(IN)                      :: M1
      type(Mat2D)                                  :: Mat2DXDble

      Mat2DXDble%XX = D1 * M1%XX
      Mat2DXDble%XY = D1 * M1%XY
      Mat2DXDble%YX = D1 * M1%YX
      Mat2DXDble%YY = D1 * M1%YY
   end function Mat2DXDble

   function DbleXMatS2D(D1, M1)
      PetscReal, intent(IN)                        :: D1
      type(MatS2D), intent(IN)                     :: M1
      type(MatS2D)                                 :: DbleXMatS2D

      DbleXMatS2D%XX = D1 * M1%XX
      DbleXMatS2D%YY = D1 * M1%YY
      DbleXMatS2D%XY = D1 * M1%XY
   end function DbleXMatS2D

   function MatS2DXDble(M1, D1)
      PetscReal, intent(IN)                        :: D1
      type(MatS2D), intent(IN)                     :: M1
      type(MatS2D)                                 :: MatS2DXDble

      MatS2DXDble%XX = D1 * M1%XX
      MatS2DXDble%YY = D1 * M1%YY
      MatS2DXDble%XY = D1 * M1%XY
   end function MatS2DXDble

   function DbleXMat3D(D1, M1)
      PetscReal, intent(IN)                        :: D1
      type(Mat3D), intent(IN)                      :: M1
      type(Mat3D)                                  :: DbleXMat3D

      DbleXMat3D%XX = D1 * M1%XX
      DbleXMat3D%XY = D1 * M1%XY
      DbleXMat3D%XZ = D1 * M1%XZ
      DbleXMat3D%YX = D1 * M1%YX
      DbleXMat3D%YY = D1 * M1%YY
      DbleXMat3D%YZ = D1 * M1%YZ
      DbleXMat3D%ZX = D1 * M1%ZX
      DbleXMat3D%ZY = D1 * M1%ZY
      DbleXMat3D%ZZ = D1 * M1%ZZ
   end function DbleXMat3D

   function Mat3DXDble(M1, D1)
      PetscReal, intent(IN)                        :: D1
      type(Mat3D), intent(IN)                      :: M1
      type(Mat3D)                                  :: Mat3DXDble

      Mat3DXDble%XX = D1 * M1%XX
      Mat3DXDble%XY = D1 * M1%XY
      Mat3DXDble%XZ = D1 * M1%XZ
      Mat3DXDble%YX = D1 * M1%YX
      Mat3DXDble%YY = D1 * M1%YY
      Mat3DXDble%YZ = D1 * M1%YZ
      Mat3DXDble%ZX = D1 * M1%ZX
      Mat3DXDble%ZY = D1 * M1%ZY
      Mat3DXDble%ZZ = D1 * M1%ZZ
   end function Mat3DXDble

   function DbleXMatS3D(D1, M1)
      PetscReal, intent(IN)                        :: D1
      type(MatS3D), intent(IN)                     :: M1
      type(MatS3D)                                 :: DbleXMatS3D

      DbleXMatS3D%XX = D1 * M1%XX
      DbleXMatS3D%YY = D1 * M1%YY
      DbleXMatS3D%ZZ = D1 * M1%ZZ
      DbleXMatS3D%YZ = D1 * M1%YZ
      DbleXMatS3D%XZ = D1 * M1%XZ
      DbleXMatS3D%XY = D1 * M1%XY
   end function DbleXMatS3D

   function MatS3DXDble(M1, D1)
      real(Kind=Kr), intent(IN)                    :: D1
      type(MatS3D), intent(IN)                     :: M1
      type(MatS3D)                                 :: MatS3DXDble

      MatS3DXDble%XX = D1 * M1%XX
      MatS3DXDble%YY = D1 * M1%YY
      MatS3DXDble%ZZ = D1 * M1%ZZ
      MatS3DXDble%YZ = D1 * M1%YZ
      MatS3DXDble%XZ = D1 * M1%XZ
      MatS3DXDble%XY = D1 * M1%XY
   end function MatS3DXDble

   function DbleXTens4OS2D(D1, T1)
      PetscReal, intent(IN)                        :: D1
      type(Tens4OS2D), intent(IN)                  :: T1
      type(Tens4OS2D)                              :: DbleXTens4OS2D

      DbleXTens4OS2D%XXXX = D1 * T1%XXXX
      DbleXTens4OS2D%XXXY = D1 * T1%XXXY
      DbleXTens4OS2D%XXYY = D1 * T1%XXYY

      DbleXTens4OS2D%XYXY = D1 * T1%XYXY
      DbleXTens4OS2D%YYXY = D1 * T1%YYXY

      DbleXTens4OS2D%YYYY = D1 * T1%YYYY
   end function DbleXTens4OS2D

   function DbleXTens4OS3D(D1, T1)
      PetscReal, intent(IN)                        :: D1
      type(Tens4OS3D), intent(IN)                  :: T1
      type(Tens4OS3D)                              :: DbleXTens4OS3D

      DbleXTens4OS3D%XXXX = D1 * T1%XXXX
      DbleXTens4OS3D%XXYY = D1 * T1%XXYY
      DbleXTens4OS3D%XXZZ = D1 * T1%XXZZ
      DbleXTens4OS3D%XXYZ = D1 * T1%XXYZ
      DbleXTens4OS3D%XXXZ = D1 * T1%XXXZ
      DbleXTens4OS3D%XXXY = D1 * T1%XXXY

      DbleXTens4OS3D%YYYY = D1 * T1%YYYY
      DbleXTens4OS3D%YYZZ = D1 * T1%YYZZ
      DbleXTens4OS3D%YYYZ = D1 * T1%YYYZ
      DbleXTens4OS3D%YYXZ = D1 * T1%YYXZ
      DbleXTens4OS3D%YYXY = D1 * T1%YYXY

      DbleXTens4OS3D%ZZZZ = D1 * T1%ZZZZ
      DbleXTens4OS3D%ZZYZ = D1 * T1%ZZYZ
      DbleXTens4OS3D%ZZXZ = D1 * T1%ZZXZ
      DbleXTens4OS3D%ZZXY = D1 * T1%ZZXY

      DbleXTens4OS3D%YZYZ = D1 * T1%YZYZ
      DbleXTens4OS3D%YZXZ = D1 * T1%YZXZ
      DbleXTens4OS3D%YZXY = D1 * T1%YZXY

      DbleXTens4OS3D%XZXZ = D1 * T1%XZXZ
      DbleXTens4OS3D%XZXY = D1 * T1%XZXY

      DbleXTens4OS3D%XYXY = D1 * T1%XYXY
   end function DbleXTens4OS3D

   function Tens4OS2DXDble(T1, D1)
      PetscReal, intent(IN)                        :: D1
      type(Tens4OS2D), intent(IN)                  :: T1
      type(Tens4OS2D)                              :: Tens4OS2DXDble

      Tens4OS2DXDble%XXXX = D1 * T1%XXXX
      Tens4OS2DXDble%XXXY = D1 * T1%XXXY
      Tens4OS2DXDble%XXYY = D1 * T1%XXYY

      Tens4OS2DXDble%XYXY = D1 * T1%XYXY
      Tens4OS2DXDble%YYXY = D1 * T1%YYXY

      Tens4OS2DXDble%YYYY = D1 * T1%YYYY
   end function Tens4OS2DXDble

   function Tens4OS3DXDble(T1, D1)
      type(Tens4OS3D), intent(IN)                  :: T1
      PetscReal, intent(IN)                        :: D1
      type(Tens4OS3D)                              :: Tens4OS3DXDble

      Tens4OS3DXDble%XXXX = D1 * T1%XXXX
      Tens4OS3DXDble%XXYY = D1 * T1%XXYY
      Tens4OS3DXDble%XXZZ = D1 * T1%XXZZ
      Tens4OS3DXDble%XXYZ = D1 * T1%XXYZ
      Tens4OS3DXDble%XXXZ = D1 * T1%XXXZ
      Tens4OS3DXDble%XXXY = D1 * T1%XXXY

      Tens4OS3DXDble%YYYY = D1 * T1%YYYY
      Tens4OS3DXDble%YYZZ = D1 * T1%YYZZ
      Tens4OS3DXDble%YYYZ = D1 * T1%YYYZ
      Tens4OS3DXDble%YYXZ = D1 * T1%YYXZ
      Tens4OS3DXDble%YYXY = D1 * T1%YYXY

      Tens4OS3DXDble%ZZZZ = D1 * T1%ZZZZ
      Tens4OS3DXDble%ZZYZ = D1 * T1%ZZYZ
      Tens4OS3DXDble%ZZXZ = D1 * T1%ZZXZ
      Tens4OS3DXDble%ZZXY = D1 * T1%ZZXY

      Tens4OS3DXDble%YZYZ = D1 * T1%YZYZ
      Tens4OS3DXDble%YZXZ = D1 * T1%YZXZ
      Tens4OS3DXDble%YZXY = D1 * T1%YZXY

      Tens4OS3DXDble%XZXZ = D1 * T1%XZXZ
      Tens4OS3DXDble%XZXY = D1 * T1%XZXY

      Tens4OS3DXDble%XYXY = D1 * T1%XYXY
   end function Tens4OS3DXDble

   function MatXVect2D(M1, V1)
      type(Mat2D), intent(IN)                      :: M1
      type(Vect2D), intent(IN)                     :: V1
      type(Vect2D)                                 :: MatXVect2D

      MatXVect2D%X = M1%XX * V1%X + M1%XY * V1%Y
      MatXVect2D%Y = M1%YX * V1%X + M1%YY * V1%Y
   end function MatXVect2D

   function MatXVect2DS(M1, V1)
      type(MatS2D), intent(IN)                     :: M1
      type(Vect2D), intent(IN)                     :: V1
      type(Vect2D)                                 :: MatXVect2DS

      MatXVect2DS%X = M1%XX * V1%X + M1%XY * V1%Y
      MatXVect2DS%Y = M1%XY * V1%X + M1%YY * V1%Y
   end function MatXVect2DS

   function MatXVect3D(M1, V1)
      type(Mat3D), intent(IN)                      :: M1
      type(Vect3D), intent(IN)                     :: V1
      type(Vect3D)                                 :: MatXVect3D

      MatXVect3D%X = M1%XX * V1%X + M1%XY * V1%Y + M1%XZ * V1%Z
      MatXVect3D%Y = M1%YX * V1%X + M1%YY * V1%Y + M1%YZ * V1%Z
      MatXVect3D%Z = M1%ZX * V1%X + M1%ZY * V1%Y + M1%ZZ * V1%Z
   end function MatXVect3D

   function MatXVect3DS(M1, V1)
      type(MatS3D), intent(IN)                     :: M1
      type(Vect3D), intent(IN)                     :: V1
      type(Vect3D)                                 :: MatXVect3DS

      MatXVect3DS%X = M1%XX * V1%X + M1%XY * V1%Y + M1%XZ * V1%Z
      MatXVect3DS%Y = M1%XY * V1%X + M1%YY * V1%Y + M1%YZ * V1%Z
      MatXVect3DS%Z = M1%XZ * V1%X + M1%YZ * V1%Y + M1%ZZ * V1%Z
   end function MatXVect3DS

   function Tens4OS2DXMatS2D(T1, M1)
      type(Tens4OS2D), intent(IN)                  :: T1
      type(MatS2D), intent(IN)                     :: M1
      type(MatS2D)                                 :: Tens4OS2DXMatS2D

      Tens4OS2DXMatS2D%XX = T1%XXXX * M1%XX + T1%XXYY * M1%YY + T1%XXXY * M1%XY * 2.0_kr
      Tens4OS2DXMatS2D%YY = T1%XXYY * M1%XX + T1%YYYY * M1%YY + T1%YYXY * M1%XY * 2.0_kr
      Tens4OS2DXMatS2D%XY = T1%XXXY * M1%XX + T1%YYXY * M1%YY + T1%XYXY * M1%XY * 2.0_kr
   end function Tens4OS2DXMatS2D

   function Tens4OS3DXMatS3D(T1, M1)
      type(Tens4OS3D), intent(IN)                  :: T1
      type(MatS3D), intent(IN)                     :: M1
      type(MatS3D)                                 :: Tens4OS3DXMatS3D

      Tens4OS3DXMatS3D%XX = T1%XXXX * M1%XX + T1%XXYY * M1%YY + T1%XXZZ * M1%ZZ &
                              + (T1%XXYZ * M1%YZ + T1%XXXZ * M1%XZ + T1%XXXY * M1%XY) * 2.0_kr
      Tens4OS3DXMatS3D%YY = T1%XXYY * M1%XX + T1%YYYY * M1%YY + T1%YYZZ * M1%ZZ &
                              + (T1%YYYZ * M1%YZ + T1%YYXZ * M1%XZ + T1%YYXY * M1%XY) * 2.0_kr
      Tens4OS3DXMatS3D%ZZ = T1%XXZZ * M1%XX + T1%YYZZ * M1%YY + T1%ZZZZ * M1%ZZ &
                              + (T1%ZZYZ * M1%YZ + T1%ZZXZ * M1%XZ + T1%ZZXY * M1%XY) * 2.0_kr
      Tens4OS3DXMatS3D%YZ = T1%XXYZ * M1%XX + T1%YYYZ * M1%YY + T1%ZZYZ * M1%ZZ &
                              + (T1%YZYZ * M1%YZ + T1%YZXZ * M1%XZ + T1%YZXY * M1%XY) * 2.0_kr
      Tens4OS3DXMatS3D%XZ = T1%XXXZ * M1%XX + T1%YYXZ * M1%YY + T1%ZZXZ * M1%ZZ &
                              + (T1%YZXZ * M1%YZ + T1%XZXZ * M1%XZ + T1%XZXY * M1%XY) * 2.0_kr
      Tens4OS3DXMatS3D%XY = T1%XXXY * M1%XX + T1%YYXY * M1%YY + T1%ZZXY * M1%ZZ &
                              + (T1%YZXY * M1%YZ + T1%XZXY * M1%XZ + T1%XYXY * M1%XY) * 2.0_kr
   end function Tens4OS3DXMatS3D

#undef __FUNCT__
#define __FUNCT__ "Tens4OS2DXMat2D"
!!!
!!!
!!!  Tens4OS2DXMat2D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!
   function Tens4OS2DXMat2D(T1, M1)
      type(Tens4OS2D), intent(IN)                  :: T1
      type(Mat2D), intent(IN)                      :: M1
      type(MatS2D)                                 :: Tens4OS2DXMat2D

      Tens4OS2DXMat2D = T1 * symmetrize(M1)
   end function Tens4OS2DXMat2D

#undef __FUNCT__
#define __FUNCT__ "Tens4OS3DXMat3D"
!!!
!!!
!!!  Tens4OS3DXMat3D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!
   function Tens4OS3DXMat3D(T1, M1)
      type(Tens4OS3D), intent(IN)                  :: T1
      type(Mat3D), intent(IN)                      :: M1
      type(MatS3D)                                 :: Tens4OS3DXMat3D

      Tens4OS3DXMat3D = T1 * symmetrize(M1)
   end function Tens4OS3DXMat3D

   function Mat2DXMat2D(M1, M2)
      type(Mat2D), intent(IN)                      :: M1, M2
      type(Mat2D)                                 :: Mat2DXMat2D

      Mat2DXMat2D%XX = M1%XX * M2%XX + M1%XY * M2%YX
      Mat2DXMat2D%XY = M1%XX * M2%XY + M1%XY * M2%YY
      Mat2DXMat2D%YX = M1%YX * M2%XX + M1%YY * M2%YX
      Mat2DXMat2D%YY = M1%YX * M2%XY + M1%YY * M2%YY
   end function Mat2DXMat2D

   function Mat2DXMatS2D(M1, M2)
      type(Mat2D), intent(IN)                      :: M1
      type(MatS2D), intent(IN)                     :: M2
      type(Mat2D)                                  :: Mat2DXMatS2D

      Mat2DXMatS2D%XX = M1%XX * M2%XX + M1%XY * M2%XY
      Mat2DXMatS2D%XY = M1%XX * M2%XY + M1%XY * M2%YY
      Mat2DXMatS2D%YX = M1%YX * M2%XX + M1%YY * M2%XY
      Mat2DXMatS2D%YY = M1%YX * M2%XY + M1%YY * M2%YY
   end function Mat2DXMatS2D

   function MatS2DXMatS2D(M1, M2)
      type(MatS2D), intent(IN)                     :: M1, M2
      type(Mat2D)                                  :: MatS2DXMatS2D

      MatS2DXMatS2D%XX = M1%XX * M2%XX + M1%XY * M2%XY
      MatS2DXMatS2D%XY = M1%XX * M2%XY + M1%XY * M2%YY
      MatS2DXMatS2D%YX = M1%XY * M2%XX + M1%YY * M2%XY
      MatS2DXMatS2D%YY = M1%XY * M2%XY + M1%YY * M2%YY
   end function MatS2DXMatS2D

   function Mat3DXMat3D(M1, M2)
      type(Mat3D), intent(IN)                      :: M1, M2
      type(Mat3D)                                  :: Mat3DXMat3D

      Mat3DXMat3D%XX = M1%XX * M2%XX + M1%XY * M2%YX + M1%XZ * M2%ZX
      Mat3DXMat3D%XY = M1%XX * M2%XY + M1%XY * M2%YY + M1%XZ * M2%ZY
      Mat3DXMat3D%XZ = M1%XX * M2%XZ + M1%XY * M2%YZ + M1%XZ * M2%ZZ
      Mat3DXMat3D%YX = M1%YX * M2%XX + M1%YY * M2%YX + M1%YZ * M2%ZX
      Mat3DXMat3D%YY = M1%YX * M2%XY + M1%YY * M2%YY + M1%YZ * M2%ZY
      Mat3DXMat3D%YZ = M1%YX * M2%XZ + M1%YY * M2%YZ + M1%YZ * M2%ZZ
      Mat3DXMat3D%ZX = M1%ZX * M2%XX + M1%ZY * M2%YX + M1%ZZ * M2%ZX
      Mat3DXMat3D%ZY = M1%ZX * M2%XY + M1%ZY * M2%YY + M1%ZZ * M2%ZY
      Mat3DXMat3D%ZZ = M1%ZX * M2%XZ + M1%ZY * M2%YZ + M1%ZZ * M2%ZZ
   end function Mat3DXMat3D

   function Mat3DXMatS3D(M1, M2)
      type(Mat3D), intent(IN)                      :: M1
      type(MatS3D), intent(IN)                     :: M2
      type(Mat3D)                                  :: Mat3DXMatS3D

      Mat3DXMatS3D%XX = M1%XX * M2%XX + M1%XY * M2%XY + M1%XZ * M2%XZ
      Mat3DXMatS3D%XY = M1%XX * M2%XY + M1%XY * M2%YY + M1%XZ * M2%YZ
      Mat3DXMatS3D%XZ = M1%XX * M2%XZ + M1%XY * M2%YZ + M1%XZ * M2%ZZ
      Mat3DXMatS3D%YX = M1%YX * M2%XX + M1%YY * M2%XY + M1%YZ * M2%XZ
      Mat3DXMatS3D%YY = M1%YX * M2%XY + M1%YY * M2%YY + M1%YZ * M2%YZ
      Mat3DXMatS3D%YZ = M1%YX * M2%XZ + M1%YY * M2%YZ + M1%YZ * M2%ZZ
      Mat3DXMatS3D%ZX = M1%ZX * M2%XX + M1%ZY * M2%XY + M1%ZZ * M2%XZ
      Mat3DXMatS3D%ZY = M1%ZX * M2%XY + M1%ZY * M2%YY + M1%ZZ * M2%YZ
      Mat3DXMatS3D%ZZ = M1%ZX * M2%XZ + M1%ZY * M2%YZ + M1%ZZ * M2%ZZ
   end function Mat3DXMatS3D

   function MatS3DXMatS3D(M1, M2)
      type(MatS3D), intent(IN)                     :: M1, M2
      type(Mat3D)                                  :: MatS3DXMatS3D

      MatS3DXMatS3D%XX = M1%XX * M2%XX + M1%XY * M2%XY + M1%XZ * M2%XZ
      MatS3DXMatS3D%XY = M1%XX * M2%XY + M1%XY * M2%YY + M1%XZ * M2%YZ
      MatS3DXMatS3D%XZ = M1%XX * M2%XZ + M1%XY * M2%YZ + M1%XZ * M2%ZZ
      MatS3DXMatS3D%YX = M1%XY * M2%XX + M1%YY * M2%XY + M1%YZ * M2%XZ
      MatS3DXMatS3D%YY = M1%XY * M2%XY + M1%YY * M2%YY + M1%YZ * M2%YZ
      MatS3DXMatS3D%YZ = M1%XY * M2%XZ + M1%YY * M2%YZ + M1%YZ * M2%ZZ
      MatS3DXMatS3D%ZX = M1%XZ * M2%XX + M1%YZ * M2%XY + M1%ZZ * M2%XZ
      MatS3DXMatS3D%ZY = M1%XZ * M2%XY + M1%YZ * M2%YY + M1%ZZ * M2%YZ
      MatS3DXMatS3D%ZZ = M1%XZ * M2%XZ + M1%YZ * M2%YZ + M1%ZZ * M2%ZZ
   end function MatS3DXMatS3D

   ! Overloading "/"
   function Vect2DQuot(V1, D1)
      type(Vect2D), intent(IN)                     :: V1
      PetscReal, intent(IN)                        :: D1
      type(Vect2D)                                 :: Vect2DQuot

      Vect2DQuot%X = V1%X / D1
      Vect2DQuot%Y = V1%Y / D1
   end function Vect2DQuot

   function Vect3DQuot(V1, D1)
      type(Vect3D), intent(IN)                     :: V1
      PetscReal, intent(IN)                        :: D1
      type(Vect3D)                                 :: Vect3DQuot

      Vect3DQuot%X = V1%X / D1
      Vect3DQuot%Y = V1%Y / D1
      Vect3DQuot%Z = V1%Z / D1
   end function Vect3DQuot

   function Mat2DQuot(M1, D1)
      PetscReal, intent(IN)                        :: D1
      type(Mat2D), intent(IN)                      :: M1
      type(Mat2D)                                  :: Mat2DQuot

      Mat2DQuot%XX = M1%XX / D1
      Mat2DQuot%XY = M1%XY / D1
      Mat2DQuot%YX = M1%YX / D1
      Mat2DQuot%YY = M1%YY / D1
   end function Mat2DQuot

   function MatS2DQuot(M1, D1)
      PetscReal, intent(IN)                        :: D1
      type(MatS2D), intent(IN)                     :: M1
      type(MatS2D)                                 :: MatS2DQuot

      MatS2DQuot%XX = M1%XX / D1
      MatS2DQuot%YY = M1%YY / D1
      MatS2DQuot%XY = M1%XY / D1
   end function MatS2DQuot

   function Mat3DQuot(M1, D1)
      PetscReal, intent(IN)                        :: D1
      type(Mat3D), intent(IN)                      :: M1
      type(Mat3D)                                  :: Mat3DQuot

      Mat3DQuot%XX = M1%XX / D1
      Mat3DQuot%XY = M1%XY / D1
      Mat3DQuot%XZ = M1%XZ / D1
      Mat3DQuot%YX = M1%YX / D1
      Mat3DQuot%YY = M1%YY / D1
      Mat3DQuot%YZ = M1%YZ / D1
      Mat3DQuot%ZX = M1%ZX / D1
      Mat3DQuot%ZY = M1%ZY / D1
      Mat3DQuot%ZZ = M1%ZZ / D1
   end function Mat3DQuot

   function MatS3DQuot(M1, D1)
      PetscReal, intent(IN)                        :: D1
      type(MatS3D), intent(IN)                     :: M1
      type(MatS3D)                                 :: MatS3DQuot

      MatS3DQuot%XX = M1%XX / D1
      MatS3DQuot%YY = M1%YY / D1
      MatS3DQuot%ZZ = M1%ZZ / D1
      MatS3DQuot%YZ = M1%YZ / D1
      MatS3DQuot%XZ = M1%XZ / D1
      MatS3DQuot%XY = M1%XY / D1
   end function MatS3DQuot

   function Tens4OS2DQuot(T1, D1)
      type(Tens4OS2D), intent(IN)                  :: T1
      PetscReal, intent(IN)                        :: D1
      type(Tens4OS2D)                              :: Tens4OS2DQuot

      Tens4OS2DQuot%XXXX = T1%XXXX / D1
      Tens4OS2DQuot%XXXY = T1%XXXY / D1
      Tens4OS2DQuot%XXYY = T1%XXYY / D1

      Tens4OS2DQuot%XYXY = T1%XYXY / D1
      Tens4OS2DQuot%YYXY = T1%YYXY / D1

      Tens4OS2DQuot%YYYY = T1%YYYY / D1
   end function Tens4OS2DQuot

   function Tens4OS3DQuot(T1, D1)
      type(Tens4OS3D), intent(IN)                  :: T1
      PetscReal, intent(IN)                        :: D1
      type(Tens4OS3D)                              :: Tens4OS3DQuot

      Tens4OS3DQuot%XXXX = T1%XXXX / D1
      Tens4OS3DQuot%XXYY = T1%XXYY / D1
      Tens4OS3DQuot%XXZZ = T1%XXZZ / D1
      Tens4OS3DQuot%XXYZ = T1%XXYZ / D1
      Tens4OS3DQuot%XXXZ = T1%XXXZ / D1
      Tens4OS3DQuot%XXXY = T1%XXXY / D1

      Tens4OS3DQuot%YYYY = T1%YYYY / D1
      Tens4OS3DQuot%YYZZ = T1%YYZZ / D1
      Tens4OS3DQuot%YYYZ = T1%YYYZ / D1
      Tens4OS3DQuot%YYXZ = T1%YYXZ / D1
      Tens4OS3DQuot%YYXY = T1%YYXY / D1

      Tens4OS3DQuot%ZZZZ = T1%ZZZZ / D1
      Tens4OS3DQuot%ZZYZ = T1%ZZYZ / D1
      Tens4OS3DQuot%ZZXZ = T1%ZZXZ / D1
      Tens4OS3DQuot%ZZXY = T1%ZZXY / D1

      Tens4OS3DQuot%YZYZ = T1%YZYZ / D1
      Tens4OS3DQuot%YZXZ = T1%YZXZ / D1
      Tens4OS3DQuot%YZXY = T1%YZXY / D1

      Tens4OS3DQuot%XZXZ = T1%XZXZ / D1
      Tens4OS3DQuot%XZXY = T1%XZXY / D1

      Tens4OS3DQuot%XYXY = T1%XYXY / D1
   end function Tens4OS3DQuot

   ! dot product in 2D and 3D
   function DotP2D(V1, V2)
      type(Vect2D), intent(IN)                     :: V1, V2
      PetscReal                                    :: DotP2D

      DotP2D = V1%X * V2%X + V1%Y * V2%Y
   end function DotP2D

   function DotP3D(V1, V2)
      type(Vect3D), intent(IN)                     :: V1, V2
      PetscReal                                    :: DotP3D

      DotP3D = V1%X * V2%X + V1%Y * V2%Y + V1%Z * V2%Z
   end function DotP3D

   function ContP2D(M1, M2)
      ! tr(A^t x B)
      type(Mat2D), intent(IN)                      :: M1, M2
      PetscReal                                    :: ContP2D

      ContP2D = M1%XX * M2%XX + M1%XY * M2%XY + M1%YX * M2%YX + M1%YY * M2%YY
   end function ContP2D

   function ContP2DS(M1, M2)
      ! tr(A^t x B)
      type(MatS2D), intent(IN)                     :: M1, M2
      PetscReal                                    :: ContP2DS

      ContP2DS = M1%XX * M2%XX + M1%YY * M2%YY + 2.0_kr * M1%XY * M2%XY
   end function ContP2DS

#undef __FUNCT__
#define __FUNCT__ "Mat2DDotMatS2D"
!!!
!!!
!!!  Mat2DDotMatS2D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!
   function Mat2DDotMatS2D(M1, M2)
      type(Mat2D), intent(IN)                      :: M1
      type(MatS2D), intent(IN)                     :: M2
      PetscReal                                    :: Mat2DDotMatS2D

      Mat2DDotMatS2D = symmetrize(M1) .DotP.M2
   end function Mat2DDotMatS2D

#undef __FUNCT__
#define __FUNCT__ "MatS2DDotMat2D"
!!!
!!!
!!!  MatS2DDotMat2D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!
   function MatS2DDotMat2D(M1, M2)
      type(MatS2D), intent(IN)                     :: M1
      type(Mat2D), intent(IN)                      :: M2
      PetscReal                                    :: MatS2DDotMat2D

      MatS2DDotMat2D = M1.DotP.symmetrize(M2)
   end function MatS2DDotMat2D

#undef __FUNCT__
#define __FUNCT__ "Mat2DDotMatS3D"
!!!
!!!
!!!  Mat3DDotMatS3D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!
   function Mat3DDotMatS3D(M1, M2)
      type(Mat3D), intent(IN)                      :: M1
      type(MatS3D), intent(IN)                     :: M2
      PetscReal                                    :: Mat3DDotMatS3D

      Mat3DDotMatS3D = symmetrize(M1) .DotP.M2
   end function Mat3DDotMatS3D

#undef __FUNCT__
#define __FUNCT__ "MatS3DDotMat3D"
!!!
!!!
!!!  MatS3DDotMat3D:
!!!
!!!  (c) 2016 Blaise Bourdin bourdin@lsu.edu
!!!
   function MatS3DDotMat3D(M1, M2)
      type(MatS3D), intent(IN)                     :: M1
      type(Mat3D), intent(IN)                      :: M2
      PetscReal                                    :: MatS3DDotMat3D

      MatS3DDotMat3D = M1.DotP.symmetrize(M2)
   end function MatS3DDotMat3D

   function ContP3D(M1, M2)
      ! tr(A^t x B)
      type(Mat3D), intent(IN)                      :: M1, M2
      PetscReal                                   :: ContP3D

      ContP3D = M1%XX * M2%XX + M1%XY * M2%XY + M1%XZ * M2%XZ + &
           &    M1%YX * M2%YX + M1%YY * M2%YY + M1%YZ * M2%YZ + &
           &    M1%ZX * M2%ZX + M1%ZY * M2%ZY + M1%ZZ * M2%ZZ
   end function ContP3D

   function ContP3DS(M1, M2)
      ! tr(A^t x B)
      type(MatS3D), intent(IN)                     :: M1, M2
      PetscReal                                    :: ContP3DS

      ContP3DS = M1%XX * M2%XX + M1%YY * M2%YY + M1%ZZ * M2%ZZ + 2.0_kr * M1%YZ * M2%YZ &
                 + 2.0_kr * M1%XZ * M2%XZ + 2.0_kr * M1%XY * M2%XY
   end function ContP3DS

   ! cross product 3D
   function CrossP3D(V1, V2)
      type(Vect3D), intent(IN)                     :: V1, V2
      type(Vect3D)                                 :: CrossP3D

      CrossP3D%X = V1%Y * V2%Z - V1%Z * V2%Y
      CrossP3D%Y = V1%Z * V2%X - V1%X * V2%Z
      CrossP3D%Z = V1%X * V2%Y - V1%Y * V2%X
   end function CrossP3D

   ! Transpose
   function Transpose2D(M1)
      type(Mat2D), intent(IN)                     :: M1
      type(Mat2D)                                 :: Transpose2D

      Transpose2D%XX = M1%XX
      Transpose2D%XY = M1%YX
      Transpose2D%YX = M1%XY
      Transpose2D%YY = M1%YY
   end function Transpose2D

   function Transpose3D(M1)
      type(Mat3D), intent(IN)                      :: M1
      type(Mat3D)                                  :: Transpose3D

      Transpose3D%XX = M1%XX
      Transpose3D%XY = M1%YX
      Transpose3D%XZ = M1%ZX
      Transpose3D%YX = M1%XY
      Transpose3D%YY = M1%YY
      Transpose3D%YZ = M1%ZY
      Transpose3D%ZX = M1%XZ
      Transpose3D%ZY = M1%YZ
      Transpose3D%ZZ = M1%ZZ
   end function Transpose3D

   ! Tensor product
   function TensPVect2D(V1, V2)
      type(Vect2D), intent(IN)                     :: V1
      type(Vect2D), intent(IN)                     :: V2
      type(Mat2D)                                  :: TensPVect2D

      TensPVect2D%XX = V1%X * V2%X
      TensPVect2D%XY = V1%X * V2%Y
      TensPVect2D%YX = V1%Y * V2%X
      TensPVect2D%YY = V1%Y * V2%Y
   end function TensPVect2D

   function TensPVect3D(V1, V2)
      type(Vect3D), intent(IN)                     :: V1, V2
      type(Mat3D)                                  :: TensPVect3D

      TensPVect3D%XX = V1%X * V2%X
      TensPVect3D%XY = V1%X * V2%Y
      TensPVect3D%XZ = V1%X * V2%Z
      TensPVect3D%YX = V1%Y * V2%X
      TensPVect3D%YY = V1%Y * V2%Y
      TensPVect3D%YZ = V1%Y * V2%Z
      TensPVect3D%ZX = V1%Z * V2%X
      TensPVect3D%ZY = V1%Z * V2%Y
      TensPVect3D%ZZ = V1%Z * V2%Z
   end function TensPVect3D

   ! Symmetrized product
   function SymPVect2D(V1, V2)
      type(Vect2D), intent(IN)                     :: V1, V2
      type(MatS2D)                                 :: SymPVect2D

      SymPVect2D = Symmetrize(V1.TensP.V2)
   end function SymPVect2D

   function SymPVect3D(V1, V2)
      type(Vect3D), intent(IN)                     :: V1, V2
      type(MatS3D)                                 :: SymPVect3D

      SymPVect3D = Symmetrize(V1.TensP.V2)
   end function SymPVect3D

   function SymPMatS2D(M1, M2)
      type(MatS2D), intent(IN)                     :: M1, M2
      type(Tens4OS2D)                              :: SymPMatS2D

      SymPMatS2D%XXXX = M1%XX * M2%XX
      SymPMatS2D%XXYY = (M1%XX * M2%YY + M1%YY * M2%XX) * 0.5_kr
      SymPMatS2D%XXXY = (M1%XX * M2%XY + M1%XY * M2%XX) * 0.5_kr

      SymPMatS2D%YYYY = M1%YY * M2%YY
      SymPMatS2D%YYXY = (M1%YY * M2%XY + M1%XY * M2%YY) * 0.5_kr

      SymPMatS2D%XYXY = M1%XY * M2%XY
   end function SymPMatS2D

   function SymPMatS3D(M1, M2)
      type(MatS3D), intent(IN)                     :: M1, M2
      type(Tens4OS3D)                              :: SymPMatS3D

      SymPMatS3D%XXXX = M1%XX * M2%XX
      SymPMatS3D%XXYY = (M1%XX * M2%YY + M1%YY * M2%XX) * 0.5_kr
      SymPMatS3D%XXZZ = (M1%XX * M2%ZZ + M1%ZZ * M2%XX) * 0.5_kr
      SymPMatS3D%XXYZ = (M1%XX * M2%YZ + M1%YZ * M2%XX) * 0.5_kr
      SymPMatS3D%XXXZ = (M1%XX * M2%XZ + M1%XZ * M2%XX) * 0.5_kr
      SymPMatS3D%XXXY = (M1%XX * M2%XY + M1%XY * M2%XX) * 0.5_kr

      SymPMatS3D%YYYY = M1%YY * M2%YY
      SymPMatS3D%YYZZ = (M1%YY * M2%ZZ + M1%ZZ * M2%YY) * 0.5_kr
      SymPMatS3D%YYYZ = (M1%YY * M2%YZ + M1%YZ * M2%YY) * 0.5_kr
      SymPMatS3D%YYXZ = (M1%YY * M2%XZ + M1%XZ * M2%YY) * 0.5_kr
      SymPMatS3D%YYXY = (M1%YY * M2%XY + M1%XY * M2%YY) * 0.5_kr

      SymPMatS3D%ZZZZ = M1%ZZ * M2%ZZ
      SymPMatS3D%ZZYZ = (M1%ZZ * M2%YZ + M1%YZ * M2%ZZ) * 0.5_kr
      SymPMatS3D%ZZXZ = (M1%ZZ * M2%XZ + M1%XZ * M2%ZZ) * 0.5_kr
      SymPMatS3D%ZZXY = (M1%ZZ * M2%XY + M1%XY * M2%ZZ) * 0.5_kr

      SymPMatS3D%YZYZ = M1%YZ * M2%YZ
      SymPMatS3D%YZXZ = (M1%YZ * M2%XZ + M1%XZ * M2%YZ) * 0.5_kr
      SymPMatS3D%YZXY = (M1%YZ * M2%XY + M1%XY * M2%YZ) * 0.5_kr

      SymPMatS3D%XZXZ = M1%XZ * M2%XZ
      SymPMatS3D%XZXY = (M1%XZ * M2%XY + M1%XY * M2%XZ) * 0.5_kr

      SymPMatS3D%XYXY = M1%XY * M2%XY
   end function SymPMatS3D

   function oDotMatS2D(M1, M2)
      type(MatS2D), intent(IN)                     :: M1, M2
      type(Tens4OS2D)                              :: oDotMatS2D

      oDotMatS2D%XXXX = M1%XX * M2%XX
      oDotMatS2D%XXYY = M1%XY * M2%XY
      oDotMatS2D%XXXY = (M1%XX * M2%XY + M1%XY * M2%XX) * 0.5_kr

      oDotMatS2D%YYYY = M1%YY * M2%YY
      oDotMatS2D%YYXY = (M1%XY * M2%YY + M1%YY * M2%XY) * 0.5_kr

      oDotMatS2D%XYXY = (M1%XX * M2%YY + M1%YY * M2%XX) * 0.5_kr
   end function oDotMatS2D

   function oDotMatS3D(M1, M2)
      type(MatS3D), intent(IN)                     :: M1, M2
      type(Tens4OS3D)                              :: oDotMatS3D

      oDotMatS3D%XXXX = M1%XX * M2%XX
      oDotMatS3D%XXYY = M1%XY * M2%XY
      oDotMatS3D%XXZZ = M1%XZ * M2%XZ
      oDotMatS3D%XXYZ = (M1%XY * M2%XZ + M1%XZ * M2%YZ) * 0.5_kr
      oDotMatS3D%XXXZ = (M1%XX * M2%XZ + M1%XZ * M2%XX) * 0.5_kr
      oDotMatS3D%XXXY = (M1%XX * M2%XY + M1%XY * M2%XX) * 0.5_kr

      oDotMatS3D%YYYY = M1%YY * M2%YY
      oDotMatS3D%YYZZ = M1%YZ * M2%YZ
      oDotMatS3D%YYYZ = (M1%YY * M2%YZ + M1%YZ * M2%YY) * 0.5_kr
      oDotMatS3D%YYXZ = (M1%XY * M2%YZ + M1%YZ * M2%XY) * 0.5_kr
      oDotMatS3D%YYXY = (M1%XY * M2%YY + M1%YY * M2%XY) * 0.5_kr

      oDotMatS3D%ZZZZ = M1%ZZ * M2%ZZ
      oDotMatS3D%ZZYZ = (M1%YZ * M2%ZZ + M1%ZZ * M2%YZ) * 0.5_kr
      oDotMatS3D%ZZXZ = (M1%XZ * M2%ZZ + M1%ZZ * M2%XZ) * 0.5_kr
      oDotMatS3D%ZZXY = (M1%XZ * M2%YZ + M1%YZ * M2%XZ) * 0.5_kr

      oDotMatS3D%YZYZ = (M1%YY * M2%ZZ + M1%ZZ * M2%YY) * 0.5_kr
      oDotMatS3D%YZXZ = (M1%XY * M2%ZZ + M1%ZZ * M2%XY) * 0.5_kr
      oDotMatS3D%YZXY = (M1%XY * M2%YZ + M1%YZ * M2%XY) * 0.5_kr

      oDotMatS3D%XZXZ = (M1%XX * M2%ZZ + M1%ZZ * M2%XX) * 0.5_kr
      oDotMatS3D%XZXY = (M1%XX * M2%YZ + M1%YZ * M2%XX) * 0.5_kr

      oDotMatS3D%XYXY = (M1%XX * M2%YY + M1%YY * M2%XX) * 0.5_kr
   end function oDotMatS3D

   function Trace2D(M1)
      type(Mat2D), intent(IN)                      :: M1
      PetscReal                                    :: Trace2D

      Trace2D = M1%XX + M1%YY
   end function Trace2D

   function Trace2DS(M1)
      type(MatS2D), intent(IN)                     :: M1
      PetscReal                                    :: Trace2DS

      Trace2DS = M1%XX + M1%YY
   end function Trace2DS

   function Trace3D(M1)
      type(Mat3D), intent(IN)                      :: M1
      PetscReal                                    :: Trace3D

      Trace3D = M1%XX + M1%YY + M1%ZZ
   end function Trace3D

   function Trace3DS(M1)
      type(MatS3D), intent(IN)                     :: M1
      PetscReal                                    :: Trace3DS

      Trace3DS = M1%XX + M1%YY + M1%ZZ
   end function Trace3DS

   subroutine Vect2D_Get_Real(V1, R1)
      type(Vect2D), intent(OUT)                    :: V1
      PetscReal, intent(IN)                        :: R1

      V1%X = R1
      V1%Y = R1
   end subroutine Vect2D_Get_Real

   subroutine Vect3D_Get_Real(V1, R1)
      type(Vect3D), intent(OUT)                    :: V1
      PetscReal, intent(IN)                        :: R1

      V1%X = R1
      V1%Y = R1
      V1%Z = R1
   end subroutine Vect3D_Get_Real

   subroutine Vect2D_Get_VectR(V1, R1)
      type(Vect2D), intent(OUT)               :: V1
      PetscReal, dimension(2), intent(IN)     :: R1

      V1%X = R1(1)
      V1%Y = R1(2)
   end subroutine Vect2D_Get_VectR

   subroutine Vect3D_Get_VectR(V1, R1)
      type(Vect3D), intent(OUT)               :: V1
      PetscReal, dimension(3), intent(IN)     :: R1

      V1%X = R1(1)
      V1%Y = R1(2)
      V1%Z = R1(3)
   end subroutine Vect3D_Get_VectR

   subroutine Vect2DEQ(V1, V2)
      type(Vect2D), intent(OUT)                    :: V1
      type(Vect2D), intent(IN)                     :: V2

      V1%X = V2%X
      V1%Y = V2%Y
   end subroutine Vect2DEQ

   subroutine Vect3DEQ(V1, V2)
      type(Vect3D), intent(OUT)                    :: V1
      type(Vect3D), intent(IN)                     :: V2

      V1%X = V2%X
      V1%Y = V2%Y
      V1%Z = V2%Z
   end subroutine Vect3DEQ

   subroutine Mat2D_Get_Real(M1, R1)
      type(Mat2D), intent(OUT)                     :: M1
      PetscReal, intent(IN)                        :: R1

      M1%XX = R1; M1%XY = R1
      M1%YX = R1; M1%YY = R1
   end subroutine Mat2D_Get_Real

   subroutine Mat3D_Get_Real(M1, R1)
      type(Mat3D), intent(OUT)                     :: M1
      PetscReal, intent(IN)                        :: R1

      M1%XX = R1; M1%XY = R1; M1%XZ = R1
      M1%YX = R1; M1%YY = R1; M1%YZ = R1
      M1%ZX = R1; M1%ZY = R1; M1%ZZ = R1
   end subroutine Mat3D_Get_Real

   subroutine Mat2DEQ(M1, M2)
      type(Mat2D), intent(OUT)                     :: M1
      type(Mat2D), intent(IN)                      :: M2

      M1%XX = M2%XX; M1%XY = M2%XY
      M1%YX = M2%YX; M1%YY = M2%YY
   end subroutine Mat2DEQ

   subroutine Mat3DEQ(M1, M2)
      type(Mat3D), intent(OUT)                     :: M1
      type(Mat3D), intent(IN)                      :: M2

      M1%XX = M2%XX; M1%XY = M2%XY; M1%XZ = M2%XZ
      M1%YX = M2%YX; M1%YY = M2%YY; M1%YZ = M2%YZ
      M1%ZX = M2%ZX; M1%ZY = M2%ZY; M1%ZZ = M2%ZZ
   end subroutine Mat3DEQ

   subroutine MatS2D_Get_Real(M1, R1)
      type(MatS2D), intent(OUT)                    :: M1
      PetscReal, intent(IN)                        :: R1

      M1%XX = R1; M1%YY = R1; M1%XY = R1
   end subroutine MatS2D_Get_Real

   subroutine MatS3D_Get_Real(M1, R1)
      type(MatS3D), intent(OUT)                    :: M1
      PetscReal, intent(IN)                        :: R1

      M1%XX = R1; M1%YY = R1; M1%ZZ = R1
      M1%YZ = R1; M1%XZ = R1; M1%XY = R1
   end subroutine MatS3D_Get_Real

   subroutine MatS2DEQ(M1, M2)
      type(MatS2D), intent(OUT)                    :: M1
      type(MatS2D), intent(IN)                     :: M2

      M1%XX = M2%XX; M1%YY = M2%YY; M1%XY = M2%XY
   end subroutine MatS2DEQ

   subroutine MatS3DEQ(M1, M2)
      type(MatS3D), intent(OUT)                    :: M1
      type(MatS3D), intent(IN)                     :: M2

      M1%XX = M2%XX; M1%YY = M2%YY; M1%ZZ = M2%ZZ
      M1%YZ = M2%YZ; M1%XZ = M2%XZ; M1%XY = M2%XY
   end subroutine MatS3DEQ

   subroutine MatS2D_Get_VectR(M1, R1)
      type(MatS2D), intent(OUT)                    :: M1
      PetscReal, dimension(3), intent(IN)          :: R1

      M1%XX = R1(1)
      M1%YY = R1(2)
      M1%XY = R1(3)
   end subroutine MatS2D_Get_VectR

   subroutine MatS3D_Get_VectR(M1, R1)
      type(MatS3D), intent(OUT)                    :: M1
      PetscReal, dimension(6), intent(IN)          :: R1

      M1%XX = R1(1)
      M1%YY = R1(2)
      M1%ZZ = R1(3)
      M1%YZ = R1(4)
      M1%XZ = R1(5)
      M1%XY = R1(6)
   end subroutine MatS3D_Get_VectR

   subroutine VectR_Get_MatS2D(R1, M1)
      PetscReal, dimension(3), intent(OUT)          :: R1
      type(MatS2D), intent(IN)                      :: M1

      R1(1) = M1%XX
      R1(2) = M1%YY
      R1(3) = M1%XY
   end subroutine VectR_Get_MatS2D

   subroutine VectR_Get_MatS3D(R1, M1)
      PetscReal, dimension(6), intent(OUT)          :: R1
      type(MatS3D), intent(IN)                      :: M1

      R1(1) = M1%XX
      R1(2) = M1%YY
      R1(3) = M1%ZZ
      R1(4) = M1%YZ
      R1(5) = M1%XZ
      R1(6) = M1%XY
   end subroutine VectR_Get_MatS3D

   subroutine Mat2DGetArray(M, A)
      type(Mat2D), intent(OUT)                     :: M
      PetscReal, dimension(2, 2), intent(IN)       :: A

      M%XX = A(1, 1)
      M%XY = A(1, 2)
      M%YX = A(2, 1)
      M%YY = A(2, 2)
   end subroutine Mat2DGetArray

   subroutine MatS2DGetArray(M, A)
      type(MatS2D), intent(OUT)                    :: M
      PetscReal, dimension(2, 2), intent(IN)       :: A

      M%XX = A(1, 1)
      M%YY = A(2, 2)
      M%XY = A(1, 2)
   end subroutine MatS2DGetArray

   subroutine Mat3DGetArray(M, A)
      type(Mat3D), intent(OUT)                     :: M
      PetscReal, dimension(3, 3), intent(IN)       :: A

      M%XX = A(1, 1)
      M%XY = A(1, 2)
      M%XZ = A(1, 3)
      M%YX = A(2, 1)
      M%YY = A(2, 2)
      M%YZ = A(2, 3)
      M%ZX = A(3, 1)
      M%ZY = A(3, 2)
      M%ZZ = A(3, 3)
   end subroutine Mat3DGetArray

   subroutine MatS3DGetArray(M, A)
      type(MatS3D), intent(OUT)                    :: M
      PetscReal, dimension(3, 3), intent(IN)       :: A

      M%XX = A(1, 1)
      M%YY = A(2, 2)
      M%ZZ = A(3, 3)
      M%YZ = A(2, 3)
      M%XZ = A(1, 3)
      M%XY = A(1, 2)
   end subroutine MatS3DGetArray

   subroutine ArrayGetMat2D(A, M)
      PetscReal, dimension(2, 2), intent(OUT)        :: A
      type(Mat2D), intent(IN)                        :: M

      A(1, 1) = M%XX
      A(1, 2) = M%XY
      A(2, 1) = M%YX
      A(2, 2) = M%YY
   end subroutine ArrayGetMat2D

   subroutine ArrayGetMatS2D(A, M)
      PetscReal, dimension(2, 2), intent(OUT)        :: A
      type(MatS2D), intent(IN)                       :: M

      A(1, 1) = M%XX
      A(1, 2) = M%XY
      A(2, 1) = M%XY
      A(2, 2) = M%YY
   end subroutine ArrayGetMatS2D

   subroutine ArrayGetMat3D(A, M)
      PetscReal, dimension(3, 3), intent(OUT)        :: A
      type(Mat3D), intent(IN)                        :: M

      A(1, 1) = M%XX
      A(1, 2) = M%XY
      A(1, 3) = M%XZ
      A(2, 1) = M%YX
      A(2, 2) = M%YY
      A(2, 3) = M%YZ
      A(3, 1) = M%ZX
      A(3, 2) = M%ZY
      A(3, 3) = M%ZZ
   end subroutine ArrayGetMat3D

   subroutine ArrayGetMatS3D(A, M)
      PetscReal, dimension(3, 3), intent(OUT)        :: A
      type(MatS3D), intent(IN)                       :: M

      A(1, 1) = M%XX
      A(1, 2) = M%XY
      A(1, 3) = M%XZ
      A(2, 1) = M%XY
      A(2, 2) = M%YY
      A(2, 3) = M%YZ
      A(3, 1) = M%XZ
      A(3, 2) = M%YZ
      A(3, 3) = M%ZZ
   end subroutine ArrayGetMatS3D

   subroutine Tens4OS2DEQ(T1, T2)
      type(Tens4OS2D), intent(OUT)                 :: T1
      type(Tens4OS2D), intent(IN)                  :: T2

      T1%XXXX = T2%XXXX
      T1%XXXY = T2%XXXY
      T1%XXYY = T2%XXYY

      T1%XYXY = T2%XYXY
      T1%YYXY = T2%YYXY

      T1%YYYY = T2%YYYY
   end subroutine Tens4OS2DEQ

   subroutine Tens4OS3DEQ(T1, T2)
      type(Tens4OS3D), intent(OUT)                 :: T1
      type(Tens4OS3D), intent(IN)                  :: T2

      T1%XXXX = T2%XXXX
      T1%XXYY = T2%XXYY
      T1%XXZZ = T2%XXZZ
      T1%XXYZ = T2%XXYZ
      T1%XXXZ = T2%XXXZ
      T1%XXXY = T2%XXXY

      T1%YYYY = T2%YYYY
      T1%YYZZ = T2%YYZZ
      T1%YYYZ = T2%YYYZ
      T1%YYXZ = T2%YYXZ
      T1%YYXY = T2%YYXY

      T1%ZZZZ = T2%ZZZZ
      T1%ZZYZ = T2%ZZYZ
      T1%ZZXZ = T2%ZZXZ
      T1%ZZXY = T2%ZZXY

      T1%YZYZ = T2%YZYZ
      T1%YZXZ = T2%YZXZ
      T1%YZXY = T2%YZXY

      T1%XZXZ = T2%XZXZ
      T1%XZXY = T2%XZXY

      T1%XYXY = T2%XYXY
   end subroutine Tens4OS3DEQ

   subroutine Tens4OS2D_Get_Real(T1, D1)
      type(Tens4OS2D), intent(OUT)                 :: T1
      PetscReal, intent(IN)                        :: D1

      T1%XXXX = D1
      T1%XXXY = D1
      T1%XXYY = D1

      T1%XYXY = D1
      T1%YYXY = D1

      T1%YYYY = D1
   end subroutine Tens4OS2D_Get_Real

   subroutine Tens4OS3D_Get_Real(T1, D1)
      type(Tens4OS3D), intent(OUT)                 :: T1
      PetscReal, intent(IN)                        :: D1

      T1%XXXX = D1
      T1%XXYY = D1
      T1%XXZZ = D1
      T1%XXYZ = D1
      T1%XXXZ = D1
      T1%XXXY = D1

      T1%YYYY = D1
      T1%YYZZ = D1
      T1%YYYZ = D1
      T1%YYXZ = D1
      T1%YYXY = D1

      T1%ZZZZ = D1
      T1%ZZYZ = D1
      T1%ZZXZ = D1
      T1%ZZXY = D1

      T1%YZYZ = D1
      T1%YZXZ = D1
      T1%YZXY = D1

      T1%XZXZ = D1
      T1%XZXY = D1

      T1%XYXY = D1
   end subroutine Tens4OS3D_Get_Real

   subroutine Tens4OS2DToArray(A, T)
      PetscReal, dimension(3, 3), intent(OUT)        :: A
      type(Tens4OS2D), intent(IN)                    :: T

      A(1, 1) = T%XXXX
      A(1, 2) = T%XXYY
      A(1, 3) = T%XXXY * 2.0_kr

      A(2, 1) = A(1, 2)
      A(2, 2) = T%YYYY
      A(2, 3) = T%YYXY * 2.0_kr

      A(3, 1) = A(1, 3)
      A(3, 2) = A(2, 3)
      A(3, 3) = T%XYXY * 2.0_kr
   end subroutine Tens4OS2DToArray

   subroutine ArrayToTens4OS2D(T, A)
      type(Tens4OS2D), intent(OUT)                   :: T
      PetscReal, dimension(3, 3), intent(IN)         :: A

      T%XXXX = A(1, 1)
      T%XXYY = A(1, 2)
      T%XXXY = A(1, 3)*.5_kr

      T%YYYY = A(2, 2)
      T%YYXY = A(2, 3)*.5_kr

      T%XYXY = A(3, 3)*.5_kr
   end subroutine ArrayToTens4OS2D

   subroutine Tens4OS3DToArray(A, T)
      PetscReal, dimension(6, 6), intent(OUT)        :: A
      type(Tens4OS3D), intent(IN)                    :: T

      A(1, 1) = T%XXXX
      A(1, 2) = T%XXYY
      A(1, 3) = T%XXZZ
      A(1, 4) = T%XXYZ * 2.0_kr
      A(1, 5) = T%XXXZ * 2.0_kr
      A(1, 6) = T%XXXY * 2.0_kr

      A(2, 1) = A(1, 2)
      A(2, 2) = T%YYYY
      A(2, 3) = T%YYZZ
      A(2, 4) = T%YYYZ * 2.0_kr
      A(2, 5) = T%YYXZ * 2.0_kr
      A(2, 6) = T%YYXY * 2.0_kr

      A(3, 1) = A(1, 3)
      A(3, 2) = A(2, 3)
      A(3, 3) = T%ZZZZ
      A(3, 4) = T%ZZYZ * 2.0_kr
      A(3, 5) = T%ZZXZ * 2.0_kr
      A(3, 6) = T%ZZXY * 2.0_kr

      A(4, 1) = A(1, 4)
      A(4, 2) = A(2, 4)
      A(4, 3) = A(3, 4)
      A(4, 4) = T%YZYZ * 2.0_kr
      A(4, 5) = T%YZXZ * 2.0_kr
      A(4, 6) = T%YZXY * 2.0_kr

      A(5, 1) = A(1, 5)
      A(5, 2) = A(2, 5)
      A(5, 3) = A(3, 5)
      A(5, 4) = A(4, 5)
      A(5, 5) = T%XZXZ * 2.0_kr
      A(5, 6) = T%XZXY * 2.0_kr

      A(6, 1) = A(1, 6)
      A(6, 2) = A(2, 6)
      A(6, 3) = A(3, 6)
      A(6, 4) = A(4, 6)
      A(6, 5) = A(5, 6)
      A(6, 6) = T%XYXY * 2.0_kr
   end subroutine Tens4OS3DToArray

   subroutine ArrayToTens4OS3D(T, A)
      type(Tens4OS3D), intent(OUT)                   :: T
      PetscReal, dimension(6, 6), intent(IN)         :: A

      T%XXXX = A(1, 1)
      T%XXYY = A(1, 2)
      T%XXZZ = A(1, 3)
      T%XXYZ = A(1, 4)*.5_kr
      T%XXXZ = A(1, 5)*.5_kr
      T%XXXY = A(1, 6)*.5_kr

      T%YYYY = A(2, 2)
      T%YYZZ = A(2, 3)
      T%YYYZ = A(2, 4)*.5_kr
      T%YYXZ = A(2, 5)*.5_kr
      T%YYXY = A(2, 6)*.5_kr

      T%ZZZZ = A(3, 3)
      T%ZZYZ = A(3, 4)*.5_kr
      T%ZZXZ = A(3, 5)*.5_kr
      T%ZZXY = A(3, 6)*.5_kr

      T%YZYZ = A(4, 4)*.5_kr
      T%YZXZ = A(4, 5)*.5_kr
      T%YZXY = A(4, 6)*.5_kr

      T%XZXZ = A(5, 5)*.5_kr
      T%XZXY = A(5, 6)*.5_kr

      T%XYXY = A(6, 6)*.5_kr
   end subroutine ArrayToTens4OS3D

   !!! Overloading euclidian norm of derived types
   PetscReal function Vect2DNorm(V)
      type(Vect2D), intent(IN)                     :: V

      Vect2DNorm = sqrt(V%X**2 + V%Y**2)
   end function Vect2DNorm

   PetscReal function Vect3DNorm(V)
      type(Vect3D), intent(IN)                     :: V

      Vect3DNorm = sqrt(V%X**2 + V%Y**2 + V%Z**2)
   end function Vect3DNorm

   PetscReal function Mat2DNorm(M)
      type(Mat2D), intent(IN)                      :: M

      Mat2DNorm = sqrt(M%XX**2 + M%XY**2 + M%YX**2 + M%YY**2)
   end function Mat2DNorm

   PetscReal function MatS2DNorm(M)
      type(MatS2D), intent(IN)                     :: M

      MatS2DNorm = sqrt(M%XX**2 + 2.0_kr * M%XY**2 + M%YY**2)
   end function MatS2DNorm

   PetscReal function Mat3DNorm(M)
      type(Mat3D), intent(IN)                      :: M

      Mat3DNorm = sqrt(M%XX**2 + M%XY**2 + M%XZ**2 + M%YX**2 + M%YY**2 + M%YZ**2 + M%ZX**2 + M%ZY**2 + M%ZZ**2)
   end function Mat3DNorm

   PetscReal function MatS3DNorm(M)
      type(MatS3D), intent(IN)                     :: M

      MatS3DNorm = sqrt(M%XX**2 + 2.0_kr * M%XY**2 + 2.0_kr * M%XZ**2 + M%YY**2 + 2.0_kr * M%YZ**2 + M%ZZ**2)
   end function MatS3DNorm

   function Symmetrize2D(M1)
      type(Mat2D), intent(IN)                      :: M1
      type(MatS2D)                                 :: Symmetrize2D

      Symmetrize2D%XX = M1%XX
      Symmetrize2D%YY = M1%YY
      Symmetrize2D%XY = (M1%XY + M1%YX) * 0.5_kr
   end function Symmetrize2D

   function Symmetrize3D(M1)
      type(Mat3D), intent(IN)                      :: M1
      type(MatS3D)                                 :: Symmetrize3D

      Symmetrize3D%XX = M1%XX
      Symmetrize3D%YY = M1%YY
      Symmetrize3D%ZZ = M1%ZZ

      Symmetrize3D%YZ = (M1%YZ + M1%ZY) * 0.5_kr
      Symmetrize3D%XZ = (M1%XZ + M1%ZX) * 0.5_kr
      Symmetrize3D%XY = (M1%XY + M1%YX) * 0.5_kr
   end function Symmetrize3D

   subroutine MatS2DToMat2D(M1, M2)
      type(Mat2D), intent(OUT)                     :: M1
      type(MatS2D), intent(IN)                     :: M2

      M1%XX = M2%XX
      M1%XY = M2%XY
      M1%YX = M2%XY
      M1%YY = M2%YY
   end subroutine MatS2DToMat2D

   subroutine MatS3DToMat3D(M1, M2)
      type(Mat3D), intent(OUT)                     :: M1
      type(MatS3D), intent(IN)                     :: M2

      M1%XX = M2%XX
      M1%XY = M2%XY
      M1%XZ = M2%XZ
      M1%YX = M2%XY
      M1%YY = M2%YY
      M1%YZ = M2%YZ
      M1%ZX = M2%XZ
      M1%ZY = M2%YZ
      M1%ZZ = M2%ZZ
   end subroutine MatS3DToMat3D

   subroutine Mat2DToMatS2D(M1, M2)
      type(MatS2D), intent(OUT)                    :: M1
      type(Mat2D), intent(IN)                      :: M2

      M1%XX = M2%XX
      M1%XY = M2%XY
      M1%YY = M2%YY
   end subroutine Mat2DToMatS2D

   subroutine Mat3DToMatS3D(M1, M2)
      type(MatS3D), intent(OUT)                    :: M1
      type(Mat3D), intent(IN)                      :: M2

      M1%XX = M2%XX
      M1%YY = M2%YY
      M1%ZZ = M2%ZZ
      M1%YZ = M2%YZ
      M1%XZ = M2%XZ
      M1%XY = M2%XY
   end subroutine Mat3DToMatS3D

   function RARtMat2D(A, R)
      ! A <- R.A.R^T
      type(Mat2D), intent(IN)                      :: A
      type(Mat2D), intent(IN)                      :: R
      type(Mat2D)                                  :: RARtMat2D

      RARtMat2D = R * A * transpose(R)
   end function RaRtMat2D

   function RARtMatS2D(A, R)
      ! A <- R.A.R^T
      type(MatS2D), intent(IN)                     :: A
      type(Mat2D), intent(IN)                      :: R
      type(MatS2D)                                 :: RARtMatS2D

      RARtMatS2D = R * A * transpose(R)
   end function RaRtMatS2D

   function RARtMat3D(A, R)
      ! A <- R.A.R^T
      type(Mat3D), intent(IN)                      :: A
      type(Mat3D), intent(IN)                      :: R
      type(Mat3D)                                  :: RARtMat3D

      RARtMat3D = R * A * transpose(R)
   end function RaRtMat3D

   function RARtMatS3D(A, R)
      ! A <- R.A.R^T
      type(MatS3D), intent(IN)                     :: A
      type(Mat3D), intent(IN)                      :: R
      type(MatS3D)                                 :: RARtMatS3D

      RARtMatS3D = R * A * transpose(R)
   end function RaRtMatS3D

   function RtARMat2D(A, R)
      ! A <- Rt.A.R
      type(Mat2D), intent(IN)                      :: A
      type(Mat2D), intent(IN)                      :: R
      type(Mat2D)                                  :: RtARMat2D

      RtARMat2D = transpose(R) * A * R
   end function RtaRMat2D

   function RtaRMatS2D(A, R)
      ! A <- R^T.A.R
      type(MatS2D), intent(IN)                     :: A
      type(Mat2D), intent(IN)                      :: R
      type(MatS2D)                                 :: RtaRMatS2D

      RtaRMatS2D = transpose(R) * A * R
   end function RtaRMatS2D

   function RtaRMat3D(A, R)
      ! A <- R^T.A.R
      type(Mat3D), intent(IN)                      :: A
      type(Mat3D), intent(IN)                      :: R
      type(Mat3D)                                  :: RtaRMat3D

      RtaRMat3D = transpose(R) * A * R
   end function RtaRMat3D

   function RtaRMatS3D(A, R)
      ! A <- R^T.A.R
      type(MatS3D), intent(IN)                     :: A
      type(Mat3D), intent(IN)                      :: R
      type(MatS3D)                                 :: RtaRMatS3D

      RtaRMatS3D = transpose(R) * A * R
   end function RtaRMatS3D

   function DeviatoricPart2D(M1)
      type(Mat2D), intent(IN)                      :: M1
      type(Mat2D)                                  :: DeviatoricPart2D

      PetscReal                                   :: M1_Trace

      M1_Trace = Trace(M1)
      DeviatoricPart2D%XX = M1%XX - M1_Trace * 0.5_kr
      DeviatoricPart2D%XY = M1%XY
      DeviatoricPart2D%YX = M1%YX
      DeviatoricPart2D%YY = M1%YY - M1_Trace * 0.5_kr
   end function DeviatoricPart2D

   function DeviatoricPart2DS(M1)
      type(MatS2D), intent(IN)                     :: M1
      type(MatS2D)                                 :: DeviatoricPart2DS

      PetscReal                                    :: M1_Trace

      M1_Trace = Trace(M1)
      DeviatoricPart2DS%XX = M1%XX - M1_Trace * 0.5_kr
      DeviatoricPart2DS%YY = M1%YY - M1_Trace * 0.5_kr
      DeviatoricPart2DS%XY = M1%XY
   end function DeviatoricPart2DS

   function DeviatoricPart3D(M1)
      type(Mat3D), intent(IN)                      :: M1
      type(Mat3D)                                  :: DeviatoricPart3D

      PetscReal                                    :: M1_Trace

      M1_Trace = Trace(M1)
      DeviatoricPart3D%XX = M1%XX - M1_Trace / 3.0_kr
      DeviatoricPart3D%XY = M1%XY
      DeviatoricPart3D%XZ = M1%XZ

      DeviatoricPart3D%YX = M1%YX
      DeviatoricPart3D%YY = M1%YY - M1_Trace / 3.0_kr
      DeviatoricPart3D%YZ = M1%YZ

      DeviatoricPart3D%ZX = M1%ZX
      DeviatoricPart3D%ZY = M1%ZY
      DeviatoricPart3D%ZZ = M1%ZZ - M1_Trace / 3.0_kr
   end function DeviatoricPart3D

   function DeviatoricPart3DS(M1)
      type(MatS3D), intent(IN)                     :: M1
      type(MatS3D)                                 :: DeviatoricPart3DS

      PetscReal                                    :: M1_Trace

      M1_Trace = Trace(M1)
      DeviatoricPart3DS%XX = M1%XX - M1_Trace / 3.0_kr
      DeviatoricPart3DS%YY = M1%YY - M1_Trace / 3.0_kr
      DeviatoricPart3DS%ZZ = M1%ZZ - M1_Trace / 3.0_kr

      DeviatoricPart3DS%YZ = M1%YZ
      DeviatoricPart3DS%XZ = M1%XZ
      DeviatoricPart3DS%XY = M1%XY
   end function DeviatoricPart3DS

   function HydrostaticPart2D(M1)
      type(Mat2D), intent(IN)                      :: M1
      type(Mat2D)                                  :: HydrostaticPart2D

      PetscReal                                    :: M1_Trace

      M1_Trace = Trace(M1)
      HydrostaticPart2D%XX = M1_Trace * 0.5_kr
      HydrostaticPart2D%XY = 0.0_kr
      HydrostaticPart2D%YX = 0.0_kr
      HydrostaticPart2D%YY = M1_Trace * 0.5_kr
   end function HydrostaticPart2D

   function HydrostaticPart2DS(M1)
      type(MatS2D), intent(IN)                     :: M1
      type(MatS2D)                                 :: HydrostaticPart2DS

      PetscReal                                    :: M1_Trace

      M1_Trace = Trace(M1)
      HydrostaticPart2DS%XX = M1_Trace * 0.5_kr
      HydrostaticPart2DS%YY = M1_Trace * 0.5_kr
      HydrostaticPart2DS%XY = 0.0_kr
   end function HydrostaticPart2DS

   function HydrostaticPart3D(M1)
      type(Mat3D), intent(IN)                      :: M1
      type(Mat3D)                                  :: HydrostaticPart3D

      PetscReal                                    :: M1_Trace

      M1_Trace = Trace(M1)
      HydrostaticPart3D%XX = M1_Trace / 3.0_kr
      HydrostaticPart3D%XY = 0.0_kr
      HydrostaticPart3D%XZ = 0.0_kr

      HydrostaticPart3D%YX = 0.0_kr
      HydrostaticPart3D%YY = M1_Trace / 3.0_kr
      HydrostaticPart3D%YZ = 0.0_kr

      HydrostaticPart3D%ZX = 0.0_kr
      HydrostaticPart3D%ZY = 0.0_kr
      HydrostaticPart3D%ZZ = M1_Trace / 3.0_kr
   end function HydrostaticPart3D

   function HydrostaticPart3DS(M1)
      type(MatS3D), intent(IN)                     :: M1
      type(MatS3D)                                 :: HydrostaticPart3DS

      PetscReal                                    :: M1_Trace

      M1_Trace = Trace(M1)
      HydrostaticPart3DS%XX = M1_Trace / 3.0_kr
      HydrostaticPart3DS%YY = M1_Trace / 3.0_kr
      HydrostaticPart3DS%ZZ = M1_Trace / 3.0_kr

      HydrostaticPart3DS%YZ = 0.0_kr
      HydrostaticPart3DS%XZ = 0.0_kr
      HydrostaticPart3DS%XY = 0.0_kr
   end function HydrostaticPart3DS

!====================================================================
!             END OF OPERATOR OVERLOADING
!====================================================================

   function Vol_Tetra_3D(V1, V2, V3, V4)
      type(Vect3D), intent(IN)                     :: V1, V2, V3, V4
      PetscReal                                    :: Vol_Tetra_3D

      type(Vect3D)                                 :: C1, C2, C3

      C1 = V1 - V4
      C2 = V2 - V4
      C3 = V3 - V4
      Vol_Tetra_3D = abs(C1%X * (C2%Y * C3%Z - C2%Z * C3%Y) - C1%Y * (C2%X * C3%Z - C2%Z * C3%X) &
                         + C1%Z * (C2%X * C3%Y - C2%Y * C3%X)) / 6.0_kr
   end function Vol_Tetra_3D

   function Area_Tri_2D(S1, S2, S3)
      type(Vect2D), intent(IN)                     :: S1, S2, S3
      PetscReal                                    :: Area_Tri_2D

      type(Vect2D)                                 :: C1, C2

      C1 = S2 - S1
      C2 = S3 - S1
      Area_Tri_2D = abs(C1%X * C2%Y - C1%Y * C2%X) * 0.5_kr
   end function Area_Tri_2D

   function Ht_Min_Tri_2D(S1, S2, S3)
      type(Vect2D), intent(IN)                    :: S1, S2, S3
      PetscReal                                   :: Ht_Min_Tri_2D

      type(Vect2D)                                :: C1, C2, C3
      PetscReal                                   :: H1, H2, H3, AreaX2

      C1 = S2 - S3
      C2 = S3 - S1
      C3 = S1 - S2

      AreaX2 = abs(C1%X * C2%Y - C1%Y * C2%X)

      H1 = AreaX2 / sqrt((C1.DotP.C1))
      H2 = AreaX2 / sqrt((C2.DotP.C2))
      H3 = AreaX2 / sqrt((C3.DotP.C3))

      Ht_Min_Tri_2D = min(H1, H2, H3)
   end function Ht_Min_Tri_2D

   function DetMat2D(M)
      type(Mat2D), intent(IN)                     :: M
      PetscReal                                   :: DetMat2D

      DetMat2D = M%XX * M%YY - M%XY * M%YX
   end function DetMat2D

   function DetMatS2D(M)
      type(MatS2D), intent(IN)                    :: M
      PetscReal                                   :: DetMatS2D

      DetMatS2D = M%XX * M%YY - M%XY * M%XY
   end function DetMatS2D

   function DetMat3D(M)
      type(Mat3D), intent(IN)                     :: M
      PetscReal                                   :: DetMat3D

      DetMat3D = M%XX * (M%YY * M%ZZ - M%ZY * M%YZ) - M%YX * (M%XY * M%ZZ - M%ZY * M%XZ) + M%ZX * (M%XY * M%YZ - M%YY * M%XZ)
   end function DetMat3D

   function DetMatS3D(M)
      type(MatS3D), intent(IN)                    :: M
      PetscReal                                   :: DetMatS3D

      DetMatS3D = M%XX * (M%YY * M%ZZ - M%YZ * M%YZ) - M%XY * (M%XY * M%ZZ - M%YZ * M%XZ) + M%XZ * (M%XY * M%YZ - M%YY * M%XZ)
   end function DetMatS3D

   function InvertMat2D(M)
      type(Mat2D), intent(IN)                     :: M
      type(Mat2D)                                 :: InvertMat2D

      type(Mat2D)                                 :: CofMt
      PetscReal                                   :: DetM

      DetM = M%XX * M%YY - M%XY * M%YX
      CofMt%XX = M%YY
      CofMt%XY = -M%XY
      CofMt%YX = -M%YX
      CofMt%YY = M%XX

      InvertMat2D = CofMt / DetM
   end function InvertMat2D

   function InvertMatS2D(M)
      type(MatS2D), intent(IN)                    :: M
      type(MatS2D)                                :: InvertMatS2D

      type(MatS2D)                                :: CofMt
      PetscReal                                   :: DetM

      DetM = M%XX * M%YY - M%XY**2
      CofMt%XX = M%YY
      CofMt%XY = -M%XY
      CofMt%YY = M%XX

      InvertMatS2D = CofMt / DetM
   end function InvertMatS2D

   function InvertMat3D(M)
      type(Mat3D), intent(IN)                     :: M
      type(Mat3D)                                 :: InvertMat3D

      type(Mat3D)                                 :: CofMt
      PetscReal                                   :: DetM

      DetM = M%XX * (M%YY * M%ZZ - M%ZY * M%YZ) - M%YX * (M%XY * M%ZZ - M%ZY * M%XZ) + M%ZX * (M%XY * M%YZ - M%YY * M%XZ)

      CofMt%XX = M%YY * M%ZZ - M%ZY * M%YZ
      CofMt%YX = -(M%YX * M%ZZ - M%ZX * M%YZ)
      CofMt%ZX = M%YX * M%ZY - M%ZX * M%YY
      CofMt%XY = -(M%XY * M%ZZ - M%ZY * M%XZ)
      CofMt%YY = M%XX * M%ZZ - M%ZX * M%XZ
      CofMt%ZY = -(M%XX * M%ZY - M%ZX * M%XY)
      CofMt%XZ = M%XY * M%YZ - M%YY * M%XZ
      CofMt%YZ = -(M%XX * M%YZ - M%YX * M%XZ)
      CofMt%ZZ = M%XX * M%YY - M%YX * M%XY

      InvertMat3D = CofMt / DetM
   end function InvertMat3D

   function InvertMatS3D(M)
      type(MatS3D), intent(IN)                    :: M
      type(MatS3D)                                :: InvertMatS3D

      type(MatS3D)                                :: CofMt
      PetscReal                                   :: DetM

      DetM = M%XX * (M%YY * M%ZZ - M%YZ * M%YZ) - M%XY * (M%XY * M%ZZ - M%YZ * M%XZ) + M%XZ * (M%XY * M%YZ - M%YY * M%XZ)

      CofMt%XX = M%YY * M%ZZ - M%YZ * M%YZ
      CofMt%XY = -(M%XY * M%ZZ - M%XZ * M%YZ)
      CofMt%XZ = M%XY * M%YZ - M%XZ * M%YY
      CofMt%YY = M%XX * M%ZZ - M%XZ * M%XZ
      CofMt%YZ = -(M%XX * M%YZ - M%XY * M%XZ)
      CofMt%ZZ = M%XX * M%YY - M%XY * M%XY

      InvertMatS3D = CofMt / DetM
   end function InvertMatS3D

#undef __FUNCT__
#define __FUNCT__ "simplexNormal2D"
   !!!
   !!!
   !!!  simplexNormal2D: Compute the normal to a simplex in 2D
   !!!
   !!!  (c) 2013 Blaise Bourdin bourdin@lsu.edu
   !!!
   subroutine simplexNormal2D(Coord, n, ierr)
      type(Vect2D), dimension(:), pointer              :: Coord
      type(Vect2D), intent(OUT)                        :: n
      PetscErrorCode, intent(INOUT)                    :: ierr

      n = [Coord(1)%Y - Coord(2)%Y, Coord(2)%X - Coord(1)%X]
      n = n / norm(n)
   end subroutine simplexNormal2D

#undef __FUNCT__
#define __FUNCT__ "simplexNormal3D"
   !!!
   !!!
   !!!  simplexNormal3D: Compute the normal to a simplex in 3D
   !!!
   !!!  (c) 2013 Blaise Bourdin bourdin@lsu.edu
   !!!
   subroutine simplexNormal3D(Coord, n, ierr)
      type(Vect3D), dimension(:), pointer              :: Coord
      type(Vect3D), intent(OUT)                        :: n
      PetscErrorCode, intent(INOUT)                    :: ierr

      n = (Coord(2) - Coord(1)) .crossP. (Coord(1) - Coord(3))
      n = n / norm(n)
      ierr = 0
   end subroutine simplexNormal3D

   function InvertTens4OS2D(T)
      type(Tens4OS2D), intent(IN)                  :: T
      type(Tens4OS2D)                              :: InvertTens4OS2D

      integer                                      :: ierr

      PetscReal, dimension(3, 3)                   :: TmpArray
      PetscInt, dimension(3)                       :: ipiv
      PetscReal, dimension(3)                      :: work

      !!! We convert T in a matrix using Mandel notations,invert the matrix then write back in a tensor
      TmpArray = T
      call DGETRF(3, 3, TmpArray, 3, ipiv, ierr)
      call DGETRI(3, TmpArray, 3, ipiv, work, 3, ierr)
      InvertTens4OS2D = TmpArray
   end function InvertTens4OS2D

   function InvertTens4OS3D(T)
      type(Tens4OS3D), intent(IN)                  :: T
      type(Tens4OS3D)                              :: InvertTens4OS3D

      integer                                      :: ierr

      PetscReal, dimension(6, 6)                   :: TmpArray
      PetscInt, dimension(6)                       :: ipiv
      PetscReal, dimension(6)                      :: work

      !!! We convert T in a matrix using Mandel notations,invert the matrix then write back in a tensor
      TmpArray = T
      call DGETRF(6, 6, TmpArray, 6, ipiv, ierr)
      call DGETRI(6, TmpArray, 6, ipiv, work, 6, ierr)
      InvertTens4OS3D = TmpArray
   end function InvertTens4OS3D

   subroutine Tens4OS2D2Array4(A, T)
      PetscReal, dimension(2, 2, 2, 2), intent(OUT)   :: A
      type(Tens4OS2D), intent(IN)                     :: T

      A(1, 1, 1, 1) = T%XXXX
      A(1, 1, 1, 2) = T%XXXY
      A(1, 1, 2, 1) = T%XXXY
      A(1, 1, 2, 2) = T%XXYY
      A(1, 2, 1, 1) = T%XXXY
      A(1, 2, 1, 2) = T%XYXY
      A(1, 2, 2, 1) = T%XYXY
      A(1, 2, 2, 2) = T%YYXY
      A(2, 1, 1, 1) = T%XXXY
      A(2, 1, 1, 2) = T%XYXY
      A(2, 1, 2, 1) = T%XYXY
      A(2, 1, 2, 2) = T%YYXY
      A(2, 2, 1, 1) = T%XXYY
      A(2, 2, 1, 2) = T%YYXY
      A(2, 2, 2, 1) = T%YYXY
      A(2, 2, 2, 2) = T%YYYY
   end subroutine Tens4OS2D2Array4

   subroutine Tens4OS3D2Array4(A, T)
      PetscReal, dimension(3, 3, 3, 3), intent(OUT)   :: A
      type(Tens4OS3D), intent(IN)                     :: T

      A(1, 1, 1, 1) = T%XXXX
      A(1, 1, 2, 2) = T%XXYY; A(2, 2, 1, 1) = T%XXYY
      A(1, 1, 3, 3) = T%XXZZ; A(3, 3, 1, 1) = T%XXZZ
      A(1, 1, 2, 3) = T%XXYZ; A(1, 1, 3, 2) = T%XXYZ; A(2, 3, 1, 1) = T%XXYZ; A(3, 2, 1, 1) = T%XXYZ
      A(1, 1, 1, 3) = T%XXXZ; A(1, 1, 3, 1) = T%XXXZ; A(1, 3, 1, 1) = T%XXXZ; A(3, 1, 1, 1) = T%XXXZ
      A(1, 1, 1, 2) = T%XXXY; A(1, 1, 2, 1) = T%XXXY; A(1, 2, 1, 1) = T%XXXY; A(2, 1, 1, 1) = T%XXXY

      A(2, 2, 2, 2) = T%YYYY
      A(2, 2, 3, 3) = T%YYZZ; A(3, 3, 2, 2) = T%YYZZ
      A(2, 2, 2, 3) = T%YYYZ; A(2, 2, 3, 2) = T%YYYZ; A(2, 3, 2, 2) = T%YYYZ; A(3, 2, 2, 2) = T%YYYZ
      A(2, 2, 1, 3) = T%YYXZ; A(2, 2, 3, 1) = T%YYXZ; A(1, 3, 2, 2) = T%YYXZ; A(3, 1, 2, 2) = T%YYXZ
      A(2, 2, 1, 2) = T%YYXY; A(2, 2, 2, 1) = T%YYXY; A(1, 2, 2, 2) = T%YYXY; A(2, 1, 2, 2) = T%YYXY

      A(3, 3, 3, 3) = T%ZZZZ
      A(3, 3, 2, 3) = T%ZZYZ; A(3, 3, 3, 2) = T%ZZYZ; A(2, 3, 3, 3) = T%ZZYZ; A(3, 2, 3, 3) = T%ZZYZ
      A(3, 3, 1, 3) = T%ZZXZ; A(3, 3, 3, 1) = T%ZZXZ; A(1, 3, 3, 3) = T%ZZXZ; A(3, 1, 3, 3) = T%ZZXZ
      A(3, 3, 1, 2) = T%ZZXY; A(3, 3, 2, 1) = T%ZZXY; A(1, 2, 3, 3) = T%ZZXY; A(2, 1, 3, 3) = T%ZZXY

      A(2, 3, 2, 3) = T%YZYZ; A(2, 3, 3, 2) = T%YZYZ; A(3, 2, 2, 3) = T%YZYZ; A(3, 2, 3, 2) = T%YZYZ
      A(2, 3, 1, 3) = T%YZXZ; A(2, 3, 3, 1) = T%YZXZ; A(3, 2, 1, 3) = T%YZXZ; A(3, 2, 3, 1) = T%YZXZ
      A(1, 3, 2, 3) = T%YZXZ; A(1, 3, 3, 2) = T%YZXZ; A(3, 1, 2, 3) = T%YZXZ; A(3, 1, 3, 2) = T%YZXZ
      A(2, 3, 1, 2) = T%YZXY; A(2, 3, 2, 1) = T%YZXY; A(3, 2, 1, 2) = T%YZXY; A(3, 2, 2, 1) = T%YZXY
      A(1, 2, 2, 3) = T%YZXY; A(1, 2, 3, 2) = T%YZXY; A(2, 1, 2, 3) = T%YZXY; A(2, 1, 3, 2) = T%YZXY

      A(1, 2, 1, 2) = T%XYXY; A(1, 2, 2, 1) = T%XYXY; A(2, 1, 1, 2) = T%XYXY; A(2, 1, 2, 1) = T%XYXY
      A(1, 3, 1, 2) = T%XZXY; A(1, 3, 2, 1) = T%XZXY; A(3, 1, 1, 2) = T%XZXY; A(3, 1, 2, 1) = T%XZXY
      A(1, 2, 1, 3) = T%XZXY; A(1, 2, 3, 1) = T%XZXY; A(2, 1, 1, 3) = T%XZXY; A(2, 1, 3, 1) = T%XZXY
      A(1, 3, 1, 3) = T%XZXZ; A(1, 3, 3, 1) = T%XZXZ; A(3, 1, 1, 3) = T%XZXZ; A(3, 1, 3, 1) = T%XZXZ
   end subroutine Tens4OS3D2Array4

   subroutine Array42Tens4OS2D(T, A)
      type(Tens4OS2D), intent(OUT)                    :: T
      PetscReal, dimension(2, 2, 2, 2), intent(IN)    :: A

      T%XXXX = A(1, 1, 1, 1)
      T%XXXY = A(1, 1, 1, 2)
      T%XXYY = A(1, 1, 2, 2)
      T%XYXY = A(1, 2, 1, 2)
      T%YYXY = A(2, 2, 1, 2)
      T%YYYY = A(2, 2, 2, 2)
   end subroutine Array42Tens4OS2D

   subroutine Array42Tens4OS3D(T, A)
      type(Tens4OS3D), intent(OUT)                    :: T
      PetscReal, dimension(3, 3, 3, 3), intent(IN)    :: A

      T%XXXX = A(1, 1, 1, 1); T%XXYY = A(1, 1, 2, 2); T%XXZZ = A(1, 1, 3, 3); T%XXYZ = A(1, 1, 2, 3); T%XXXZ = A(1, 1, 1, 3); T%XXXY = A(1, 1, 1, 2)
      T%YYYY = A(2, 2, 2, 2); T%YYZZ = A(2, 2, 3, 3); T%YYYZ = A(2, 2, 2, 3); T%YYXZ = A(2, 2, 1, 3); T%YYXY = A(2, 2, 1, 2)
      T%ZZZZ = A(3, 3, 3, 3); T%ZZYZ = A(3, 3, 2, 3); T%ZZXZ = A(3, 3, 1, 3); T%ZZXY = A(3, 3, 1, 2)
      T%YZYZ = A(2, 3, 2, 3); T%YZXZ = A(2, 3, 1, 3); T%YZXY = A(2, 3, 1, 2)
      T%XZXZ = A(1, 3, 1, 3); T%XZXY = A(1, 3, 1, 2)
      T%XYXY = A(1, 2, 1, 2)
   end subroutine Array42Tens4OS3D

   function Tens4OS2DTransform(T, M)
      !!! Apply the transformation given by the matrix R to a 4th order tensor
      !!! i.e. C_{ijkl} = R_{ip}.R_{jq}.R_{kr}.R{ls} A_{pqrs}
      type(Tens4OS2D), intent(IN)                  :: T
      type(Mat2D), intent(IN)                      :: M

      type(Tens4OS2D)                              :: Tens4OS2DTransform

      PetscReal, dimension(2, 2, 2, 2)             :: TT, C
      PetscReal, dimension(2, 2)                   :: MM
      integer                                      :: i, j, k, l
      integer                                      :: p, q, r, s

      TT = T
      MM = M
      C = 0.0_kr
      do i = 1, 2
         do j = 1, 2
            do k = 1, 2
               do l = 1, 2
                  do p = 1, 2
                     do q = 1, 2
                        do r = 1, 2
                           do s = 1, 2
                              C(i, j, k, l) = C(i, j, k, l) + MM(i, p) * MM(j, q) * MM(k, r) * MM(l, s) * TT(p, q, r, s)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      Tens4OS2DTransform = C

   end function Tens4OS2DTransform

   function Tens4OS3DTransform(T, M)
      !!! Apply the transformation given by the matrix R to a 4th order tensor
      !!! i.e. C_{ijkl} = R_{ip}.R_{jq}.R_{kr}.R{ls} A_{pqrs}
      type(Tens4OS3D), intent(IN)                 :: T
      type(Mat3D), intent(IN)                     :: M

      type(Tens4OS3D)                             :: Tens4OS3DTransform

      PetscReal, dimension(3, 3, 3, 3)            :: TT, C
      PetscReal, dimension(3, 3)                  :: MM
      integer                                     :: i, j, k, l
      integer                                     :: p, q, r, s

      TT = T
      MM = M
      C = 0.0_kr
      do i = 1, 3
         do j = 1, 3
            do k = 1, 3
               do l = 1, 3
                  do p = 1, 3
                     do q = 1, 3
                        do r = 1, 3
                           do s = 1, 3
                              C(i, j, k, l) = C(i, j, k, l) + MM(i, p) * MM(j, q) * MM(k, r) * MM(l, s) * TT(p, q, r, s)
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do

      Tens4OS3DTransform = C
   end function Tens4OS3DTransform

   function Tens4OS2DSquareRoot(T)
      type(Tens4OS2D), intent(IN)                  :: T
      type(Tens4OS2D)                              :: Tens4OS2DSquareRoot

      integer, parameter                           :: n = 3
      integer                                      :: i, j
      PetscReal, dimension(n, n)                   :: A, Pt
      PetscReal, dimension(n)                      :: lmbda
      PetscReal                                    :: d
      PetscInt                                     :: lwork = 2 * n**2 + 6 * n + 1
      PetscReal, dimension(2*n**2 + 6*n + 1)       :: work
      PetscInt                                     :: liwork = 5 * n + 3
      PetscInt, dimension(5*n + 3)                 :: iwork
      PetscInt                                     :: info

      A = T
      call DSYEVD('V', 'L', n, A, n, lmbda, work, lwork, iwork, liwork, info)
      Pt = transpose(A)
      do i = 1, n

         if (lmbda(i) < 0.0_kr) then
            write (*, *) 'ERROR in Tens4OSDSquareRoot, negative eigenvalue ', lmbda(i)
         else
            d = sqrt(lmbda(i))
            do j = 1, n
               Pt(i, j) = A(j, i) * d
            end do
         end if
      end do
      Tens4OS2DSquareRoot = matmul(A, Pt)
   end function Tens4OS2DSquareRoot

   function Tens4OS3DSquareRoot(T)
      type(Tens4OS3D), intent(IN)                  :: T
      type(Tens4OS3D)                              :: Tens4OS3DSquareRoot

      integer, parameter                           :: n = 6
      integer                                      :: i, j
      PetscReal, dimension(n, n)                   :: A, Pt
      PetscReal, dimension(n)                      :: lmbda
      PetscReal                                    :: d
      PetscInt                                     :: lwork = 2 * n**2 + 6 * n + 1
      PetscReal, dimension(2*n**2 + 6*n + 1)       :: work
      PetscInt                                     :: liwork = 5 * n + 3
      PetscInt, dimension(5*n + 3)                 :: iwork
      PetscInt                                     :: info

      A = T
      call DSYEVD('V', 'L', n, A, n, lmbda, work, lwork, iwork, liwork, info)
      Pt = transpose(A)
      do i = 1, n
         if (lmbda(i) < 0.0_kr) then
            write (*, *) 'ERROR in Tens4OSDSquareRoot, negative eigenvalue ', lmbda(i)
         else
            d = sqrt(lmbda(i))
            do j = 1, n
               Pt(i, j) = A(j, i) * d
            end do
         end if
      end do
      Tens4OS3DSquareRoot = matmul(A, Pt)
   end function Tens4OS3DSquareRoot

   subroutine MatS2DSpectralDecomposition(M, ppleValues, ppleDirections)
      type(MatS2D), intent(IN)                      :: M
      PetscReal, dimension(2), intent(OUT)          :: ppleValues
      type(MatS2D), dimension(2), intent(OUT)       :: ppleDirections

      integer, parameter                            :: n = 2
      PetscReal, dimension(n, n)                    :: A
      PetscInt                                      :: i
      PetscInt                                      :: lwork = 2 * n**2 + 6 * n + 1
      PetscReal, dimension(2*n**2 + 6*n + 1)        :: work
      PetscInt                                      :: liwork = 5 * n + 3
      PetscInt, dimension(5*n + 3)                  :: iwork
      PetscInt                                      :: info

      A = M
      call DSYEVD('V', 'L', n, A, n, ppleValues, work, lwork, iwork, liwork, info)
      do i = 1, n
         ppleDirections(i)%XX = A(1, i)**2
         ppleDirections(i)%YY = A(2, i)**2
         ppleDirections(i)%XY = A(1, i) * A(2, i)
      end do
   end subroutine MatS2DSpectralDecomposition

   subroutine MatS3DSpectralDecomposition(M, ppleValues, ppleDirections)
      type(MatS3D), intent(IN)                      :: M
      PetscReal, dimension(3), intent(OUT)          :: ppleValues
      type(MatS3D), dimension(3), intent(OUT)       :: ppleDirections

      integer, parameter                            :: n = 3
      PetscReal, dimension(n, n)                    :: A
      PetscInt                                      :: i
      PetscInt                                      :: lwork = 2 * n**2 + 6 * n + 1
      PetscReal, dimension(2*n**2 + 6*n + 1)        :: work
      PetscInt                                      :: liwork = 5 * n + 3
      PetscInt, dimension(5*n + 3)                  :: iwork
      PetscInt                                      :: info

      A = M
      call DSYEVD('V', 'L', n, A, n, ppleValues, work, lwork, iwork, liwork, info)
      do i = 1, n
         ppleDirections(i)%XX = A(1, i)**2
         ppleDirections(i)%YY = A(2, i)**2
         ppleDirections(i)%ZZ = A(3, i)**2
         ppleDirections(i)%YZ = A(2, i) * A(3, i)
         ppleDirections(i)%XZ = A(1, i) * A(3, i)
         ppleDirections(i)%XY = A(1, i) * A(2, i)
      end do
   end subroutine MatS3DSpectralDecomposition

   subroutine MatS3DEigenVectorValues(M, MatProj, MatDiag)
      type(MatS3D), intent(IN)                     :: M
      PetscReal, dimension(3)                      :: ppleValues
      type(Mat3D), intent(OUT)                     :: MatProj
      type(MatS3D), intent(OUT)                    :: MatDiag

      integer, parameter                           :: n = 3
      PetscReal, dimension(n, n)                   :: A
      PetscInt                                     :: lwork = 2 * n**2 + 6 * n + 1
      PetscReal, dimension(2*n**2 + 6*n + 1)       :: work
      PetscInt                                     :: liwork = 5 * n + 3
      PetscInt, dimension(5*n + 3)                 :: iwork
      PetscInt                                     :: info

      A = M
      call DSYEVD('V', 'L', n, A, n, ppleValues, work, lwork, iwork, liwork, info)

      MatDiag = 0.0_kr
      MatDiag%XX = ppleValues(1)
      MatDiag%YY = ppleValues(2)
      MatDiag%ZZ = ppleValues(3)

      MatProj%XX = A(1, 1)
      MatProj%YX = A(2, 1)
      MatProj%ZX = A(3, 1)
      MatProj%XY = A(1, 2)
      MatProj%YY = A(2, 2)
      MatProj%ZY = A(3, 2)
      MatProj%XZ = A(1, 3)
      MatProj%YZ = A(2, 3)
      MatProj%ZZ = A(3, 3)
   end subroutine MatS3DEigenVectorValues

   subroutine MatS2DEigenVectorValues(M, MatProj, MatDiag)
      type(MatS2D), intent(IN)                     :: M
      PetscReal, dimension(2)                      :: ppleValues
      type(Mat2D), intent(OUT)                     :: MatProj
      type(MatS2D), intent(OUT)                    :: MatDiag

      integer, parameter                           :: n = 2
      PetscReal, dimension(n, n)                   :: A
      PetscInt                                     :: lwork = 2 * n**2 + 6 * n + 1
      PetscReal, dimension(2*n**2 + 6*n + 1)       :: work
      PetscInt                                     :: liwork = 5 * n + 3
      PetscInt, dimension(5*n + 3)                 :: iwork
      PetscInt                                     :: info

      A = M
      call DSYEVD('V', 'L', n, A, n, ppleValues, work, lwork, iwork, liwork, info)
      if (info /= 0) then
         write (*, *) 'DSYEVD failed with info=', info
         write (*, *) 'A: ', A
         write (*, *) 'ppleValues: ', ppleValues
         SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_LIB, "DSYEVD failed: "//__FUNCT__)
      end if

      MatDiag = 0.0_kr
      MatDiag%XX = ppleValues(1)
      MatDiag%YY = ppleValues(2)

      MatProj%XX = A(1, 1)
      MatProj%XY = A(1, 2)
      MatProj%YX = A(2, 1)
      MatProj%YY = A(2, 2)
   end subroutine MatS2DEigenVectorValues

#undef __FUNCT__
#define __FUNCT__ "Mat2DMoment"
!!!
!!!
!!!  Mat2DMoment: k-th moment of a 2x2 matrix, i.e tr(A^k)
!!!
!!!  (c) 2019 Blaise Bourdin bourdin@lsu.edu
!!!

   function Mat2DMoment(k, A)
      integer, intent(IN)                          :: k
      type(Mat2D), intent(IN)                      :: A
      PetscReal                                    :: Mat2DMoment

      integer                                      :: i
      type(Mat2D)                                  :: Ak

      select case (k)
      case (1)
         Mat2DMoment = trace(A)
      case (2)
         Mat2DMoment = (A%XX**2 + 2 * A%XY * A%YX + A%YY**2) / 2.0_kr
      case default
         Ak = A
         do i = 1, k - 1
            Ak = Ak * A
         end do
         Mat2DMoment = trace(Ak) / k
      end select
   end function Mat2DMoment

#undef __FUNCT__
#define __FUNCT__ "MatS2DMoment"
!!!
!!!
!!!  MatS2DMoment: k-th moment of a 2x2 symmetric matrix, i.e tr(A^k)
!!!
!!!  (c) 2019 Blaise Bourdin bourdin@lsu.edu
!!!

   function MatS2DMoment(k, A)
      integer, intent(IN)                         :: k
      type(MatS2D), intent(IN)                    :: A
      PetscReal                                   :: MatS2DMoment

      integer                                     :: i
      type(MatS2D)                                :: Ak

      select case (k)
      case (1)
         MatS2DMoment = trace(A)
      case (2)
         MatS2DMoment = (A%XX**2 + 2 * A%XY**2 + A%YY**2) / 2.0_kr
      case default
         Ak = A
         do i = 1, k - 1
            Ak = Ak * A
         end do
         MatS2DMoment = trace(Ak) / k
      end select
   end function MatS2DMoment

#undef __FUNCT__
#define __FUNCT__ "Mat3DMoment"
!!!
!!!
!!!  Mat3DMoment: k-th moment of a 3x3 matrix, i.e tr(A^k)
!!!
!!!  (c) 2019 Blaise Bourdin bourdin@lsu.edu
!!!

   function Mat3DMoment(k, A)
      integer, intent(IN)                         :: k
      type(Mat3D), intent(IN)                     :: A
      PetscReal                                   :: Mat3DMoment

      integer                                     :: i
      type(Mat3D)                                 :: Ak

      select case (k)
      case (1)
         Mat3DMoment = trace(A)
      case (2)
         Mat3DMoment = (A%XX**2 + A%YY**2 + A%ZZ**2 + 2 * A%XY * A%YX + 2 * A%XZ * A%ZX + 2 * A%YZ * A%ZY) / 2.0_kr
      case default
         Ak = A
         do i = 1, k - 1
            Ak = Ak * A
         end do
         Mat3DMoment = trace(Ak) / k
      end select
   end function Mat3DMoment

#undef __FUNCT__
#define __FUNCT__ "MatS3DMoment"
!!!
!!!
!!!  MatS3DMoment: k-th moment of a 3x3 symmetric matrix, i.e tr(A^k)
!!!
!!!  (c) 2019 Blaise Bourdin bourdin@lsu.edu
!!!

   function MatS3DMoment(k, A)
      integer, intent(IN)                         :: k
      type(MatS3D), intent(IN)                    :: A
      PetscReal                                   :: MatS3DMoment

      integer                                     :: i
      type(MatS3D)                                :: Ak

      select case (k)
      case (1)
         MatS3DMoment = trace(A)
      case (2)
         MatS3DMoment = (A%XX**2 + A%YY**2 + A%ZZ**2 + 2 * A%YZ**2 + 2 * A%XZ**2 + 2 * A%XY**2) / 2.0_kr
      case default
         Ak = A
         do i = 1, k - 1
            Ak = Ak * A
         end do
         MatS3DMoment = trace(Ak) / k
      end select
   end function MatS3DMoment
end module m_MEF90_LinAlg
