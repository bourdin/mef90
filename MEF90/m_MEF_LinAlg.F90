Module m_MEF_LinAlg
#include "finclude/petscdef.h"
   Use m_MEF_Parameters
   Use petsc
   
   IMPLICIT NONE
 
   Type Vect2D
      Sequence
      PetscReal          :: X
      PetscReal          :: Y
   End Type Vect2D
 
   Type Vect3D
      Sequence
      PetscReal          :: X
      PetscReal          :: Y
      PetscReal          :: Z
   End Type Vect3D
 
   Type Mat2D
      Sequence
      PetscReal          :: XX
      PetscReal          :: XY
      PetscReal          :: YX
      PetscReal          :: YY
   End Type Mat2D
 
   Type MatS2D
      Sequence
      PetscReal          :: XX
      PetscReal          :: YY
      PetscReal          :: XY
   End Type MatS2D
 
   Type Mat3D
      Sequence
      PetscReal          :: XX
      PetscReal          :: XY
      PetscReal          :: XZ
      PetscReal          :: YX
      PetscReal          :: YY
      PetscReal          :: YZ
      PetscReal          :: ZX
      PetscReal          :: ZY
      PetscReal          :: ZZ
   End Type Mat3D
 
   Type MatS3D
      Sequence
      PetscReal          :: XX
      PetscReal          :: YY
      PetscReal          :: ZZ
      PetscReal          :: YZ
      PetscReal          :: XZ
      PetscReal          :: XY
   End Type MatS3D
 
 !! After much hesitation,
 !! - the terms are numbered in alphabetical order (i.e. XXYX and not XYXX) 
 !! - the terms are stored in alphabetical order
   Type Tens4OS2D
      Sequence
      PetscReal          :: XXXX
      PetscReal          :: XXXY
      PetscReal          :: XXYY
      
      PetscReal          :: XYXY
      PetscReal          :: XYYY
      
      PetscReal          :: YYYY
   End Type Tens4OS2D
 
   Type Tens4OS3D
      Sequence
      PetscReal          :: XXXX
      PetscReal          :: XXXY
      PetscReal          :: XXXZ
      PetscReal          :: XXYY
      PetscReal          :: XXYZ
      PetscReal          :: XXZZ
      
      PetscReal          :: XYXY
      PetscReal          :: XYXZ
      PetscReal          :: XYYY
      PetscReal          :: XYYZ
      PetscReal          :: XYZZ
 
      PetscReal          :: XZXZ
      PetscReal          :: XZYY
      PetscReal          :: XZYZ
      PetscReal          :: XZZZ
      
      PetscReal          :: YYYY
      PetscReal          :: YYYZ
      PetscReal          :: YYZZ
      
      PetscReal          :: YZYZ
      PetscReal          :: YZZZ
      
      PetscReal          :: ZZZZ
   End Type Tens4OS3D

  Interface Operator (+)
     Module Procedure SumVect2D, SumVect3D, SumMat2D, SumMat3D, SumMatS2D, SumMatS3D, SumTens4OS2D, SumTens4OS3D
  End Interface

  Interface Operator (-)
     Module Procedure DifVect2D, DifVect3D, DifMat2D, DifMat3D,DifMatS2D, DifMatS3D, DifTens4OS2D, DifTens4OS3D
  End Interface

  Interface Operator (*)
     Module Procedure DbleXVect2D, Vect2DXDble, DbleXVect3D, Vect3DXDble, &
            DbleXMat2D, Mat2DXDble, DbleXMat3D, Mat3DXDble,               &
            DbleXMatS2D, MatS2DXDble, DbleXMatS3D, MatS3DXDble,           &
            MatXVect2D, MatXVect3D, MatXVect2DS, MatXVect3DS,             &
            DbleXTens4OS2D, Tens4OS2DXDble, Tens4OS2DXMatS2D,             &
            DbleXTens4OS3D, Tens4OS3DXDble, Tens4OS3DXMatS3D
  End Interface

  Interface Operator (/)
     Module Procedure Vect2DQuot, Vect3DQuot, Mat2DQuot, Mat3DQuot, MatS2DQuot, MatS3DQuot, Tens4OS2DQuot, Tens4OS3DQuot
  End Interface

  Interface Operator (.DotP.)
     Module Procedure DotP2D, DotP3D, ContP2D, ContP3D, ContP2DS, ContP3DS
  End Interface

  Interface Operator (.VectP.)
     Module Procedure VectP3D
  End Interface

  Interface Transpose
     Module Procedure Transpose2D, Transpose3D
  End Interface

  Interface Operator (.TensP.)
     Module Procedure TensPVect2D, TensPVect3D, TensPMat2D, TensPMat3D, TensPMatS2D, TensPMatS3D
  End Interface

  Interface Operator (.SymP.)
     Module Procedure SymPVect2D, SymPVect3D, SymPMat2D, SymPMat3D, SymPMatS2D, SymPMatS3D
  End Interface

  Interface Trace
     Module Procedure Trace2D, Trace3D, Trace2DS, Trace3DS
  End Interface

  Interface ValP
     Module Procedure ValP2D, ValP2DS
  End Interface

  Interface Assignment (=)
     Module Procedure Vect2D_Get_Real, Vect3D_Get_Real,                 &
            Vect2D_Get_VectR, Vect3D_Get_VectR,                         &
            Vect2DEQ, Vect3DEQ, Mat2D_Get_Real, Mat3D_Get_Real,         &
            Mat2DEQ, Mat3DEQ, MatS2D_Get_Real, MatS3D_Get_Real,         &
            MatS2DEQ, MatS3DEQ, Tens4OS2D_Get_Real, Tens4OS2DEQ,        &
            Tens4OS3D_Get_Real, Tens4OS3DEQ
  End Interface
  
  Interface Symmetrize
     Module Procedure Symmetrize2D, Symmetrize3D
  End Interface

!!$  Type(Vect2D), Parameter       :: e1_2D = (/ 1.0_Kr, 0.0_Kr /)
!!$  Type(Vect2D), Parameter       :: e2_2D = (/ 0.0_Kr, 1.0_Kr /)
!!$
!!$  Type(Vect3D), Parameter       :: e1_3D = (/ 1.0_Kr, 0.0_Kr, 0.0_Kr /)
!!$  Type(Vect3D), Parameter       :: e2_3D = (/ 0.0_Kr, 1.0_Kr, 0.0_Kr /)
!!$  Type(Vect3D), Parameter       :: e3_3D = (/ 0.0_Kr, 0.0_Kr, 1.0_Kr /)


Contains

  ! Fonctions pour la surcharge de l'operateur +
   Function SumVect2D (V1, V2)
      Type (Vect2D), intent(IN)                   :: V1
      Type (Vect2D), intent(IN)                   :: V2
      Type (Vect2D)                               :: SumVect2D
      
      SumVect2D%X = V1%X + V2%X
      SumVect2D%Y = V1%Y + V2%Y
   End Function SumVect2D

   Function SumVect3D (V1, V2)
      Type (Vect3D), intent(IN)                   :: V1, V2
      Type (Vect3D)                               :: SumVect3D
      
      SumVect3D%X = V1%X + V2%X 
      SumVect3D%Y = V1%Y + V2%Y
      SumVect3D%Z = V1%Z + V2%Z
   End Function SumVect3D

  Function SumMat2D (M1, M2)
    Type (Mat2D), Intent(IN)                    :: M1, M2
    Type (Mat2D)                                :: SumMat2D

    SumMat2D%XX = M1%XX + M2%XX
    SumMat2D%XY = M1%XY + M2%XY
    SumMat2D%YX = M1%YX + M2%YX
    SumMat2D%YY = M1%YY + M2%YY
  End Function SumMat2D

  Function SumMatS2D (M1, M2)
    Type (MatS2D), Intent(IN)                   :: M1, M2
    Type (MatS2D)                               :: SumMatS2D

    SumMatS2D%XX = M1%XX + M2%XX
    SumMatS2D%XY = M1%XY + M2%XY
    SumMatS2D%YY = M1%YY + M2%YY
  End Function SumMatS2D

  Function SumMat3D (M1, M2)
    Type (Mat3D), Intent(IN)                    :: M1, M2
    Type (Mat3D)                                :: SumMat3D

    SumMat3D%XX = M1%XX + M2%XX
    SumMat3D%XY = M1%XY + M2%XY
    SumMat3D%XZ = M1%XZ + M2%XZ
    SumMat3D%YX = M1%YX + M2%YX
    SumMat3D%YY = M1%YY + M2%YY
    SumMat3D%YZ = M1%YZ + M2%YZ
    SumMat3D%ZX = M1%ZX + M2%ZX
    SumMat3D%ZY = M1%ZY + M2%ZY
    SumMat3D%ZZ = M1%ZZ + M2%ZZ
  End Function SumMat3D

  Function SumMatS3D (M1, M2)
    Type (MatS3D), Intent(IN)                   :: M1, M2
    Type (MatS3D)                               :: SumMatS3D

    SumMatS3D%XX = M1%XX + M2%XX
    SumMatS3D%YY = M1%YY + M2%YY
    SumMatS3D%ZZ = M1%ZZ + M2%ZZ
    SumMatS3D%YZ = M1%YZ + M2%YZ
    SumMatS3D%XZ = M1%XZ + M2%XZ
    SumMatS3D%XY = M1%XY + M2%XY
  End Function SumMatS3D
 
  Function SumTens4OS2D (T1, T2)
    Type (Tens4OS2D), Intent(IN)                :: T1, T2
    Type (Tens4OS2D)                            :: SumTens4OS2D

    SumTens4OS2D%XXXX = T1%XXXX + T2%XXXX
    SumTens4OS2D%XXXY = T1%XXXY + T2%XXXY
    SumTens4OS2D%XXYY = T1%XXYY + T2%XXYY

    SumTens4OS2D%XYXY = T1%XYXY + T2%XYXY
    SumTens4OS2D%XYYY = T1%XYYY + T2%XYYY

    SumTens4OS2D%YYYY = T1%YYYY + T2%YYYY
  End Function SumTens4OS2D
  
  Function SumTens4OS3D (T1, T2)
    Type (Tens4OS3D), Intent(IN)                :: T1, T2
    Type (Tens4OS3D)                            :: SumTens4OS3D

    SumTens4OS3D%XXXX = T1%XXXX + T2%XXXX  
    SumTens4OS3D%XXXY = T1%XXXY + T2%XXXY
    SumTens4OS3D%XXXZ = T1%XXXZ + T2%XXXZ  
    SumTens4OS3D%XXYY = T1%XXYY + T2%XXYY  
    SumTens4OS3D%XXYZ = T1%XXYZ + T2%XXYZ  
    SumTens4OS3D%XXZZ = T1%XXZZ + T2%XXZZ  
          
    SumTens4OS3D%XYXY = T1%XYXY + T2%XYXY  
    SumTens4OS3D%XYXZ = T1%XYXZ + T2%XYXZ  
    SumTens4OS3D%XYYY = T1%XYYY + T2%XYYY  
    SumTens4OS3D%XYYZ = T1%XYYZ + T2%XYYZ  
    SumTens4OS3D%XYZZ = T1%XYZZ + T2%XYZZ  
          
    SumTens4OS3D%XZXZ = T1%XZXZ + T2%XZXZ  
    SumTens4OS3D%XZYY = T1%XZYY + T2%XZYY  
    SumTens4OS3D%XZYZ = T1%XZYZ + T2%XZYZ  
    SumTens4OS3D%XZZZ = T1%XZZZ + T2%XZZZ  
          
    SumTens4OS3D%YYYY = T1%YYYY + T2%YYYY  
    SumTens4OS3D%YYYZ = T1%YYYZ + T2%YYYZ  
    SumTens4OS3D%YYZZ = T1%YYZZ + T2%YYZZ  
          
    SumTens4OS3D%YZYZ = T1%YZYZ + T2%YZYZ  
    SumTens4OS3D%YZZZ = T1%YZZZ + T2%YZZZ  
          
    SumTens4OS3D%ZZZZ = T1%ZZZZ + T2%ZZZZ  
  End Function SumTens4OS3D

  ! Fonctions pour la surcharge de l'operateur -
  Function DifVect2D (V1, V2)
    Type (Vect2D), intent(IN)                   :: V1
    Type (Vect2D), intent(IN)                   :: V2
    Type (Vect2D)                               :: DifVect2D

    DifVect2D%X = V1%X-V2%X
    DifVect2D%Y = V1%Y-V2%Y
  End Function DifVect2D

  Function DifVect3D (V1, V2)
    Type (Vect3D), intent(IN)                   :: V1, V2
    Type (Vect3D)                               :: DifVect3D

    DifVect3D%X = V1%X - V2%X 
    DifVect3D%Y = V1%Y - V2%Y
    DifVect3D%Z = V1%Z - V2%Z
  End Function DifVect3D

  Function DifMat2D (M1, M2)
    Type (Mat2D), Intent(IN)                    :: M1, M2
    Type (Mat2D)                                :: DifMat2D

    DifMat2D%XX = M1%XX - M2%XX
    DifMat2D%XY = M1%XY - M2%XY
    DifMat2D%YX = M1%YX - M2%YX
    DifMat2D%YY = M1%YY - M2%YY
  End Function DifMat2D

  Function DifMatS2D (M1, M2)
    Type (MatS2D), Intent(IN)                   :: M1, M2
    Type (MatS2D)                               :: DifMatS2D

    DifMatS2D%XX = M1%XX - M2%XX
    DifMatS2D%YY = M1%YY - M2%YY
    DifMatS2D%XY = M1%XY - M2%XY
  End Function DifMatS2D

  Function DifMat3D (M1, M2)
    Type (Mat3D), Intent(IN)                    :: M1, M2
    Type (Mat3D)                                :: DifMat3D

    DifMat3D%XX = M1%XX - M2%XX
    DifMat3D%XY = M1%XY - M2%XY
    DifMat3D%XZ = M1%XZ - M2%XZ
    DifMat3D%YX = M1%YX - M2%YX
    DifMat3D%YY = M1%YY - M2%YY
    DifMat3D%YZ = M1%YZ - M2%YZ
    DifMat3D%ZX = M1%ZX - M2%ZX
    DifMat3D%ZY = M1%ZY - M2%ZY
    DifMat3D%ZZ = M1%ZZ - M2%ZZ
  End Function DifMat3D

  Function DifMatS3D (M1, M2)
    Type (MatS3D), Intent(IN)                   :: M1, M2
    Type (MatS3D)                               :: DifMatS3D

    DifMatS3D%XX = M1%XX - M2%XX
    DifMatS3D%YY = M1%YY - M2%YY
    DifMatS3D%ZZ = M1%ZZ - M2%ZZ
    DifMatS3D%YZ = M1%YZ - M2%YZ
    DifMatS3D%XZ = M1%XZ - M2%XZ
    DifMatS3D%XY = M1%XY - M2%XY
  End Function DifMatS3D

  Function DifTens4OS2D (T1, T2)
    Type (Tens4OS2D), Intent(IN)                :: T1, T2
    Type (Tens4OS2D)                            :: DifTens4OS2D

    DifTens4OS2D%XXXX = T1%XXXX - T2%XXXX
    DifTens4OS2D%XXXY = T1%XXXY - T2%XXXY
    DifTens4OS2D%XXYY = T1%XXYY - T2%XXYY

    DifTens4OS2D%XYXY = T1%XYXY - T2%XYXY
    DifTens4OS2D%XYYY = T1%XYYY - T2%XYYY

    DifTens4OS2D%YYYY = T1%YYYY - T2%YYYY
  End Function DifTens4OS2D

  Function DifTens4OS3D (T1, T2)
    Type (Tens4OS3D), Intent(IN)                :: T1, T2
    Type (Tens4OS3D)                            :: DifTens4OS3D

    DifTens4OS3D%XXXX = T1%XXXX - T2%XXXX  
    DifTens4OS3D%XXXY = T1%XXXY - T2%XXXY
    DifTens4OS3D%XXXZ = T1%XXXZ - T2%XXXZ  
    DifTens4OS3D%XXYY = T1%XXYY - T2%XXYY  
    DifTens4OS3D%XXYZ = T1%XXYZ - T2%XXYZ  
    DifTens4OS3D%XXZZ = T1%XXZZ - T2%XXZZ  
          
    DifTens4OS3D%XYXY = T1%XYXY - T2%XYXY  
    DifTens4OS3D%XYXZ = T1%XYXZ - T2%XYXZ  
    DifTens4OS3D%XYYY = T1%XYYY - T2%XYYY  
    DifTens4OS3D%XYYZ = T1%XYYZ - T2%XYYZ  
    DifTens4OS3D%XYZZ = T1%XYZZ - T2%XYZZ  
          
    DifTens4OS3D%XZXZ = T1%XZXZ - T2%XZXZ  
    DifTens4OS3D%XZYY = T1%XZYY - T2%XZYY  
    DifTens4OS3D%XZYZ = T1%XZYZ - T2%XZYZ  
    DifTens4OS3D%XZZZ = T1%XZZZ - T2%XZZZ  
          
    DifTens4OS3D%YYYY = T1%YYYY - T2%YYYY  
    DifTens4OS3D%YYYZ = T1%YYYZ - T2%YYYZ  
    DifTens4OS3D%YYZZ = T1%YYZZ - T2%YYZZ  
          
    DifTens4OS3D%YZYZ = T1%YZYZ - T2%YZYZ  
    DifTens4OS3D%YZZZ = T1%YZZZ - T2%YZZZ  
          
    DifTens4OS3D%ZZZZ = T1%ZZZZ - T2%ZZZZ  
  End Function DifTens4OS3D

  ! Surcharge de l'operateur *
  Function DbleXVect2D(D1, V1)
    PetscReal,       intent(IN)                 :: D1
    Type (Vect2D), intent(IN)                   :: V1
    Type (Vect2D)                               :: DbleXVect2D

    DbleXVect2D%X = D1*V1%X
    DbleXVect2D%Y = D1*V1%Y
  End Function DbleXVect2D

  Function Vect2DXDble(V1,D1)
    PetscReal,       intent(IN)                 :: D1
    Type (Vect2D), intent(IN)                   :: V1
    Type (Vect2D)                               :: Vect2DXDble

    Vect2DXDble%X = D1*V1%X
    Vect2DXDble%Y = D1*V1%Y
  End Function Vect2DXDble

  Function DbleXVect3D(D1,V1)
    PetscReal,       intent(IN)                 :: D1
    Type (Vect3D), intent(IN)                   :: V1
    Type (Vect3D)                               :: DbleXVect3D

    DbleXVect3D%X = D1*V1%X
    DbleXVect3D%Y = D1*V1%Y
    DbleXVect3D%Z = D1*V1%Z
  End Function DbleXVect3D

  Function Vect3DXDble(V1,D1)
    PetscReal,       intent(IN)                 :: D1
    Type (Vect3D), intent(IN)                   :: V1
    Type (Vect3D)                               :: Vect3DXDble

    Vect3DXDble%X = D1*V1%X
    Vect3DXDble%Y = D1*V1%Y
    Vect3DXDble%Z = D1*V1%Z
  End Function Vect3DXDble

  Function DbleXMat2D(D1, M1)
    PetscReal,       intent(IN)                 :: D1
    Type (Mat2D), intent(IN)                    :: M1
    Type (Mat2D)                                :: DbleXMat2D

    DbleXMat2D%XX = D1*M1%XX
    DbleXMat2D%XY = D1*M1%XY
    DbleXMat2D%YX = D1*M1%YX
    DbleXMat2D%YY = D1*M1%YY
  End Function DbleXMat2D

  Function Mat2DXDble(M1,D1)
    PetscReal,       intent(IN)                 :: D1
    Type (Mat2D), intent(IN)                    :: M1
    Type (Mat2D)                                :: Mat2DXDble

    Mat2DXDble%XX = D1*M1%XX
    Mat2DXDble%XY = D1*M1%XY
    Mat2DXDble%YX = D1*M1%YX
    Mat2DXDble%YY = D1*M1%YY
  End Function Mat2DXDble

  Function DbleXMatS2D(D1, M1)
    PetscReal,       intent(IN)                 :: D1
    Type (MatS2D), intent(IN)                   :: M1
    Type (MatS2D)                               :: DbleXMatS2D

    DbleXMatS2D%XX = D1*M1%XX
    DbleXMatS2D%YY = D1*M1%YY
    DbleXMatS2D%XY = D1*M1%XY
  End Function DbleXMatS2D

  Function MatS2DXDble(M1,D1)
    PetscReal,       intent(IN)                 :: D1
    Type (MatS2D), intent(IN)                   :: M1
    Type (MatS2D)                               :: MatS2DXDble

    MatS2DXDble%XX = D1*M1%XX
    MatS2DXDble%YY = D1*M1%YY
    MatS2DXDble%XY = D1*M1%XY
  End Function MatS2DXDble

  Function DbleXMat3D(D1,M1)
    PetscReal,       intent(IN)                 :: D1
    Type (Mat3D), intent(IN)                    :: M1
    Type (Mat3D)                                :: DbleXMat3D

    DbleXMat3D%XX = D1*M1%XX
    DbleXMat3D%XY = D1*M1%XY
    DbleXMat3D%XZ = D1*M1%XZ
    DbleXMat3D%YX = D1*M1%YX
    DbleXMat3D%YY = D1*M1%YY
    DbleXMat3D%YZ = D1*M1%YZ
    DbleXMat3D%ZX = D1*M1%ZX
    DbleXMat3D%ZY = D1*M1%ZY
    DbleXMat3D%ZZ = D1*M1%ZZ
  End Function DbleXMat3D

  Function Mat3DXDble(M1,D1)
    PetscReal,       intent(IN)                 :: D1
    Type (Mat3D), intent(IN)                    :: M1
    Type (Mat3D)                                :: Mat3DXDble

    Mat3DXDble%XX = D1*M1%XX
    Mat3DXDble%XY = D1*M1%XY
    Mat3DXDble%XZ = D1*M1%XZ
    Mat3DXDble%YX = D1*M1%YX
    Mat3DXDble%YY = D1*M1%YY
    Mat3DXDble%YZ = D1*M1%YZ
    Mat3DXDble%ZX = D1*M1%ZX
    Mat3DXDble%ZY = D1*M1%ZY
    Mat3DXDble%ZZ = D1*M1%ZZ
  End Function Mat3DXDble

  Function DbleXMatS3D(D1,M1)
    PetscReal,       intent(IN)                 :: D1
    Type (MatS3D), intent(IN)                   :: M1
    Type (MatS3D)                               :: DbleXMatS3D

    DbleXMatS3D%XX = D1*M1%XX
    DbleXMatS3D%YY = D1*M1%YY
    DbleXMatS3D%ZZ = D1*M1%ZZ
    DbleXMatS3D%YZ = D1*M1%YZ
    DbleXMatS3D%XZ = D1*M1%XZ
    DbleXMatS3D%XY = D1*M1%XY
  End Function DbleXMatS3D

  Function MatS3DXDble(M1,D1)
    Real (Kind = Kr), intent(IN)                :: D1
    Type (MatS3D), intent(IN)                   :: M1
    Type (MatS3D)                               :: MatS3DXDble

    MatS3DXDble%XX = D1*M1%XX
    MatS3DXDble%YY = D1*M1%YY
    MatS3DXDble%ZZ = D1*M1%ZZ
    MatS3DXDble%YZ = D1*M1%YZ
    MatS3DXDble%XZ = D1*M1%XZ
    MatS3DXDble%XY = D1*M1%XY
  End Function MatS3DXDble

  Function DbleXTens4OS2D (D1, T1)
    PetscReal,       Intent(IN)                 :: D1
    Type (Tens4OS2D), Intent(IN)                :: T1
    Type (Tens4OS2D)                            :: DbleXTens4OS2D

    DbleXTens4OS2D%XXXX = D1 * T1%XXXX
    DbleXTens4OS2D%XXXY = D1 * T1%XXXY
    DbleXTens4OS2D%XXYY = D1 * T1%XXYY

    DbleXTens4OS2D%XYXY = D1 * T1%XYXY
    DbleXTens4OS2D%XYYY = D1 * T1%XYYY

    DbleXTens4OS2D%YYYY = D1 * T1%YYYY
  End Function DbleXTens4OS2D

  Function DbleXTens4OS3D (D1, T1)
    PetscReal,       Intent(IN)                 :: D1
    Type (Tens4OS3D), Intent(IN)                :: T1
    Type (Tens4OS3D)                            :: DbleXTens4OS3D

    DbleXTens4OS3D%XXXX = D1 * T1%XXXX  
    DbleXTens4OS3D%XXXY = D1 * T1%XXXY
    DbleXTens4OS3D%XXXZ = D1 * T1%XXXZ  
    DbleXTens4OS3D%XXYY = D1 * T1%XXYY  
    DbleXTens4OS3D%XXYZ = D1 * T1%XXYZ  
    DbleXTens4OS3D%XXZZ = D1 * T1%XXZZ  
          
    DbleXTens4OS3D%XYXY = D1 * T1%XYXY  
    DbleXTens4OS3D%XYXZ = D1 * T1%XYXZ  
    DbleXTens4OS3D%XYYY = D1 * T1%XYYY  
    DbleXTens4OS3D%XYYZ = D1 * T1%XYYZ  
    DbleXTens4OS3D%XYZZ = D1 * T1%XYZZ  
          
    DbleXTens4OS3D%XZXZ = D1 * T1%XZXZ  
    DbleXTens4OS3D%XZYY = D1 * T1%XZYY  
    DbleXTens4OS3D%XZYZ = D1 * T1%XZYZ  
    DbleXTens4OS3D%XZZZ = D1 * T1%XZZZ  
          
    DbleXTens4OS3D%YYYY = D1 * T1%YYYY  
    DbleXTens4OS3D%YYYZ = D1 * T1%YYYZ  
    DbleXTens4OS3D%YYZZ = D1 * T1%YYZZ  
          
    DbleXTens4OS3D%YZYZ = D1 * T1%YZYZ  
    DbleXTens4OS3D%YZZZ = D1 * T1%YZZZ  
          
    DbleXTens4OS3D%ZZZZ = D1 * T1%ZZZZ  
  End Function DbleXTens4OS3D


  Function Tens4OS2DXDble (T1, D1)
    PetscReal,       Intent(IN)                 :: D1
    Type (Tens4OS2D), Intent(IN)                :: T1
    Type (Tens4OS2D)                            :: Tens4OS2DXDble

    Tens4OS2DXDble%XXXX = D1 * T1%XXXX
    Tens4OS2DXDble%XXXY = D1 * T1%XXXY
    Tens4OS2DXDble%XXYY = D1 * T1%XXYY

    Tens4OS2DXDble%XYXY = D1 * T1%XYXY
    Tens4OS2DXDble%XYYY = D1 * T1%XYYY

    Tens4OS2DXDble%YYYY = D1 * T1%YYYY
  End Function Tens4OS2DXDble

  Function Tens4OS3DXDble (T1, D1)
    Type (Tens4OS3D), Intent(IN)                :: T1
    PetscReal,       Intent(IN)                 :: D1
    Type (Tens4OS3D)                            :: Tens4OS3DXDble

    Tens4OS3DXDble%XXXX = D1 * T1%XXXX  
    Tens4OS3DXDble%XXXY = D1 * T1%XXXY
    Tens4OS3DXDble%XXXZ = D1 * T1%XXXZ  
    Tens4OS3DXDble%XXYY = D1 * T1%XXYY  
    Tens4OS3DXDble%XXYZ = D1 * T1%XXYZ  
    Tens4OS3DXDble%XXZZ = D1 * T1%XXZZ  
          
    Tens4OS3DXDble%XYXY = D1 * T1%XYXY  
    Tens4OS3DXDble%XYXZ = D1 * T1%XYXZ  
    Tens4OS3DXDble%XYYY = D1 * T1%XYYY  
    Tens4OS3DXDble%XYYZ = D1 * T1%XYYZ  
    Tens4OS3DXDble%XYZZ = D1 * T1%XYZZ  
          
    Tens4OS3DXDble%XZXZ = D1 * T1%XZXZ  
    Tens4OS3DXDble%XZYY = D1 * T1%XZYY  
    Tens4OS3DXDble%XZYZ = D1 * T1%XZYZ  
    Tens4OS3DXDble%XZZZ = D1 * T1%XZZZ  
          
    Tens4OS3DXDble%YYYY = D1 * T1%YYYY  
    Tens4OS3DXDble%YYYZ = D1 * T1%YYYZ  
    Tens4OS3DXDble%YYZZ = D1 * T1%YYZZ  
          
    Tens4OS3DXDble%YZYZ = D1 * T1%YZYZ  
    Tens4OS3DXDble%YZZZ = D1 * T1%YZZZ  
          
    Tens4OS3DXDble%ZZZZ = D1 * T1%ZZZZ  
  End Function Tens4OS3DXDble

  Function  MatXVect2D(M1, V1)
    Type (Mat2D), intent(IN)                    :: M1
    Type (Vect2D), intent(IN)                   :: V1
    Type (Vect2D)                               :: MatXVect2D

    MatXVect2D%X = M1%XX * V1%X + M1%XY * V1%Y
    MatXVect2D%Y = M1%YX * V1%X + M1%YY * V1%Y
  End Function MatXVect2D

  Function  MatXVect2DS(M1, V1)
    Type (MatS2D), intent(IN)                   :: M1
    Type (Vect2D), intent(IN)                   :: V1
    Type (Vect2D)                               :: MatXVect2DS

    MatXVect2DS%X = M1%XX * V1%X + M1%XY * V1%Y
    MatXVect2DS%Y = M1%XY * V1%X + M1%YY * V1%Y
  End Function MatXVect2DS

  Function  MatXVect3D(M1, V1)
    Type (Mat3D), intent(IN)                    :: M1
    Type (Vect3D), intent(IN)                   :: V1
    Type (Vect3D)                               :: MatXVect3D

    MatXVect3D%X = M1%XX * V1%X + M1%XY * V1%Y + M1%XZ * V1%Z
    MatXVect3D%Y = M1%YX * V1%X + M1%YY * V1%Y + M1%YZ * V1%Z
    MatXVect3D%Z = M1%ZX * V1%X + M1%ZY * V1%Y + M1%ZZ * V1%Z
  End Function MatXVect3D

  Function  MatXVect3DS(M1, V1)
    Type (MatS3D), intent(IN)                   :: M1
    Type (Vect3D), intent(IN)                   :: V1
    Type (Vect3D)                               :: MatXVect3DS

    MatXVect3DS%X = M1%XX * V1%X + M1%XY * V1%Y + M1%XZ * V1%Z
    MatXVect3DS%Y = M1%XY * V1%X + M1%YY * V1%Y + M1%YZ * V1%Z
    MatXVect3DS%Z = M1%XZ * V1%X + M1%YZ * V1%Y + M1%ZZ * V1%Z
  End Function MatXVect3DS

  Function Tens4OS2DXMatS2D(T1, M1)
    Type(Tens4OS2D), Intent(IN)                 :: T1
    Type(MatS2D), Intent(IN)                    :: M1
    Type(MatS2D)                                :: Tens4OS2DXMatS2D

    Tens4OS2DXMatS2D%XX = T1%XXXX * M1%XX + T1%XXYY * M1%YY + T1%XXXY * M1%XY * 2.0_Kr
    Tens4OS2DXMatS2D%YY = T1%XXYY * M1%XX + T1%YYYY * M1%YY + T1%XYYY * M1%XY * 2.0_Kr
    Tens4OS2DXMatS2D%XY = T1%XXXY * M1%XX + T1%XYYY * M1%YY + T1%XYXY * M1%XY * 2.0_Kr
  End Function Tens4OS2DXMatS2D

  Function Tens4OS3DXMatS3D(T1, M1)
    Type(Tens4OS3D), Intent(IN)                 :: T1
    Type(MatS3D), Intent(IN)                    :: M1
    Type(MatS3D)                                :: Tens4OS3DXMatS3D

    Tens4OS3DXMatS3D%XX = T1%XXXX * M1%XX + T1%XXXY * M1%XY * 2.0_Kr + T1%XXXZ * M1%XZ * 2.0_Kr                                    &
                                          + T1%XXYY * M1%YY          + T1%XXYZ * M1%YZ * 2.0_Kr                                    &
                                                                     + T1%XXZZ * M1%ZZ 
    Tens4OS3DXMatS3D%XY = T1%XXXY * M1%XX + T1%XYXY * M1%XY * 2.0_Kr + T1%XYXZ * M1%XZ * 2.0_Kr                                    &
                                          + T1%XYYY * M1%YY          + T1%XYYZ * M1%YZ * 2.0_Kr                                    &
                                                                     + T1%XYZZ * M1%ZZ 
    Tens4OS3DXMatS3D%XZ = T1%XXXZ * M1%XX + T1%XYXZ * M1%XY * 2.0_Kr + T1%XZXZ * M1%XZ * 2.0_Kr                                    &
                                          + T1%XZYY * M1%YY          + T1%XZYZ * M1%YZ * 2.0_Kr                                    &
                                                                     + T1%XZZZ * M1%ZZ 

    Tens4OS3DXMatS3D%YY = T1%XXYY * M1%XX + T1%XYYY * M1%XY * 2.0_Kr + T1%XZYY * M1%XZ * 2.0_Kr                                    &
                                          + T1%YYYY * M1%YY          + T1%YYYZ * M1%YZ * 2.0_Kr                                    &
                                                                     + T1%YYZZ * M1%ZZ 
    Tens4OS3DXMatS3D%YZ = T1%XXYZ * M1%XX + T1%XYYZ * M1%XY * 2.0_Kr + T1%XZYZ * M1%XZ * 2.0_Kr                                    &
                                          + T1%YYYZ * M1%YY          + T1%YZYZ * M1%YZ * 2.0_Kr                                    &
                                                                     + T1%YZZZ * M1%ZZ 

    Tens4OS3DXMatS3D%ZZ = T1%XXZZ * M1%XX + T1%XYZZ * M1%XY * 2.0_Kr + T1%XZZZ * M1%XZ * 2.0_Kr                                    &
                                          + T1%YYZZ * M1%YY          + T1%YZZZ * M1%YZ * 2.0_Kr                                    &
                                                                     + T1%ZZZZ * M1%ZZ 
  End Function Tens4OS3DXMatS3D

  ! Surcharge de l'operateur /
  Function Vect2DQuot(V1,D1)
    PetscReal,       intent(IN)                 :: D1
    Type (Vect2D), intent(IN)                   :: V1
    Type (Vect2D)                               :: Vect2DQuot

    Vect2DQuot%X = V1%X / D1
    Vect2DQuot%Y = V1%Y / D1
  End Function Vect2DQuot


  Function Vect3DQuot(V1,D1)
    PetscReal,       intent(IN)                 :: D1
    Type (Vect3D), intent(IN)                   :: V1
    Type (Vect3D)                               :: Vect3DQuot

    Vect3DQuot%X = V1%X / D1
    Vect3DQuot%Y = V1%Y / D1
    Vect3DQuot%Z = V1%Z / D1
  End Function Vect3DQuot

  Function Mat2DQuot(M1,D1)
    PetscReal,       intent(IN)                 :: D1
    Type (Mat2D), intent(IN)                    :: M1
    Type (Mat2D)                                :: Mat2DQuot

    Mat2DQuot%XX = M1%XX / D1
    Mat2DQuot%XY = M1%XY / D1
    Mat2DQuot%YX = M1%YX / D1
    Mat2DQuot%YY = M1%YY / D1
  End Function Mat2DQuot

  Function MatS2DQuot(M1,D1)
    PetscReal,       intent(IN)                 :: D1
    Type (MatS2D), intent(IN)                   :: M1
    Type (MatS2D)                               :: MatS2DQuot

    MatS2DQuot%XX = M1%XX / D1
    MatS2DQuot%YY = M1%YY / D1
    MatS2DQuot%XY = M1%XY / D1
  End Function MatS2DQuot

  Function Mat3DQuot(M1,D1)
    PetscReal,       intent(IN)                 :: D1
    Type (Mat3D), intent(IN)                    :: M1
    Type (Mat3D)                                :: Mat3DQuot

    Mat3DQuot%XX = M1%XX / D1
    Mat3DQuot%XY = M1%XY / D1
    Mat3DQuot%XZ = M1%XZ / D1
    Mat3DQuot%YX = M1%YX / D1
    Mat3DQuot%YY = M1%YY / D1
    Mat3DQuot%YZ = M1%YZ / D1
    Mat3DQuot%ZX = M1%ZX / D1
    Mat3DQuot%ZY = M1%ZY / D1
    Mat3DQuot%ZZ = M1%ZZ / D1
  End Function Mat3DQuot

  Function MatS3DQuot(M1,D1)
    PetscReal,       intent(IN)                 :: D1
    Type (MatS3D), intent(IN)                   :: M1
    Type (MatS3D)                               :: MatS3DQuot

    MatS3DQuot%XX = M1%XX / D1
    MatS3DQuot%YY = M1%YY / D1
    MatS3DQuot%ZZ = M1%ZZ / D1
    MatS3DQuot%YZ = M1%YZ / D1
    MatS3DQuot%XZ = M1%XZ / D1
    MatS3DQuot%XY = M1%XY / D1
  End Function MatS3DQuot

  Function Tens4OS2DQuot(T1, D1)
    Type(Tens4OS2D), Intent(IN)                 :: T1
    PetscReal,       Intent(IN)                 :: D1
    Type(Tens4OS2D)                             :: Tens4OS2DQuot

    Tens4OS2DQuot%XXXX = T1%XXXX / D1
    Tens4OS2DQuot%XXXY = T1%XXXY / D1
    Tens4OS2DQuot%XXYY = T1%XXYY / D1

    Tens4OS2DQuot%XYXY = T1%XYXY / D1
    Tens4OS2DQuot%XYYY = T1%XYYY / D1

    Tens4OS2DQuot%YYYY = T1%YYYY / D1
  End Function Tens4OS2DQuot
    

  Function Tens4OS3DQuot (T1, D1)
    Type (Tens4OS3D), Intent(IN)                :: T1
    PetscReal,       Intent(IN)                 :: D1
    Type (Tens4OS3D)                            :: Tens4OS3DQuot

    Tens4OS3DQuot%XXXX = T1%XXXX / D1  
    Tens4OS3DQuot%XXXY = T1%XXXY / D1
    Tens4OS3DQuot%XXXZ = T1%XXXZ / D1  
    Tens4OS3DQuot%XXYY = T1%XXYY / D1  
    Tens4OS3DQuot%XXYZ = T1%XXYZ / D1  
    Tens4OS3DQuot%XXZZ = T1%XXZZ / D1  
          
    Tens4OS3DQuot%XYXY = T1%XYXY / D1  
    Tens4OS3DQuot%XYXZ = T1%XYXZ / D1  
    Tens4OS3DQuot%XYYY = T1%XYYY / D1  
    Tens4OS3DQuot%XYYZ = T1%XYYZ / D1  
    Tens4OS3DQuot%XYZZ = T1%XYZZ / D1  
          
    Tens4OS3DQuot%XZXZ = T1%XZXZ / D1  
    Tens4OS3DQuot%XZYY = T1%XZYY / D1  
    Tens4OS3DQuot%XZYZ = T1%XZYZ / D1  
    Tens4OS3DQuot%XZZZ = T1%XZZZ / D1  
          
    Tens4OS3DQuot%YYYY = T1%YYYY / D1  
    Tens4OS3DQuot%YYYZ = T1%YYYZ / D1  
    Tens4OS3DQuot%YYZZ = T1%YYZZ / D1  
          
    Tens4OS3DQuot%YZYZ = T1%YZYZ / D1  
    Tens4OS3DQuot%YZZZ = T1%YZZZ / D1  
          
    Tens4OS3DQuot%ZZZZ = T1%ZZZZ / D1  
  End Function Tens4OS3DQuot
  
  ! dot product 2D et 3D
  Function DotP2D(V1, V2)
    Type (Vect2D), Intent(IN)                   :: V1, V2
    PetscReal                                   :: DotP2D
    DotP2D = V1%X * V2%X + V1%Y * V2%Y
  End Function DotP2D

  Function DotP3D(V1, V2)
    Type (Vect3D), Intent(IN)                   :: V1, V2
    PetscReal                                   :: DotP3D
    DotP3D = V1%X * V2%X + V1%Y * V2%Y + V1%Z * V2%Z
  End Function DotP3D

  Function ContP2D(M1, M2)
    ! tr(A^t x B)
    Type(Mat2D), Intent(IN)                       :: M1, M2
    PetscReal                                     :: ContP2D
    ContP2D = M1%XX * M2%XX + M1%XY * M2%XY + &
         &          M1%YX * M2%YX + M1%YY * M2%YY
  End Function ContP2D

  Function ContP2DS(M1, M2)
    ! tr(A^t x B)
    Type(MatS2D), Intent(IN)                      :: M1, M2
    PetscReal                                     :: ContP2DS
    ContP2DS = M1%XX * M2%XX + M1%YY * M2%YY + &
         &         2.0_Kr * M1%XY * M2%XY
  End Function ContP2DS

  Function ContP3D(M1, M2)
    ! tr(A^t x B)
    Type(Mat3D), Intent(IN)                       :: M1, M2
    PetscReal                                     :: ContP3D
    ContP3D = M1%XX * M2%XX + M1%XY * M2%XY + M1%XZ * M2%XZ + &
         &          M1%YX * M2%YX + M1%YY * M2%YY + M1%YZ * M2%YZ + &
         &          M1%ZX * M2%ZX + M1%ZY * M2%ZY + M1%ZZ * M2%ZZ
  End Function ContP3D

  Function ContP3DS(M1, M2)
    ! tr(A^t x B)
    Type(MatS3D), Intent(IN)                      :: M1, M2
    PetscReal                                     :: ContP3DS
    ContP3DS = M1%XX * M2%XX + M1%YY * M2%YY + M1%ZZ * M2%ZZ &
         &      + 2.0_Kr * M1%YZ * M2%YZ &
         &      + 2.0_Kr * M1%XZ * M2%XZ &
         &      + 2.0_Kr * M1%XY * M2%XY
  End Function ContP3DS

  ! cross product 3D
  Function VectP3D(V1, V2)
    Type (Vect3D), Intent(IN)                   :: V1, V2
    Type (Vect3D)                               :: VectP3D
    VectP3D%X =  V1%Y * V2%Z - V1%Z * V2%Y
    VectP3D%Y = -V1%X * V2%Z + V1%Z * V2%X
    VectP3D%Z =  V1%X * V2%Y - V1%Y * V2%X
  End Function VectP3D

  ! Transpose
  Function Transpose2D(M1)
    Type (Mat2D), Intent(IN)                    :: M1
    Type (Mat2D)                                :: Transpose2D
    Transpose2D%XX = M1%XX
    Transpose2D%XY = M1%YX
    Transpose2D%YX = M1%XY
    Transpose2D%YY = M1%YY
  End Function Transpose2D

  Function Transpose3D(M1)
    Type (Mat3D), Intent(IN)                    :: M1
    Type (Mat3D)                                :: Transpose3D
    Transpose3D%XX = M1%XX
    Transpose3D%XY = M1%YX
    Transpose3D%XZ = M1%ZX
    Transpose3D%YX = M1%XY
    Transpose3D%YY = M1%YY
    Transpose3D%YZ = M1%ZY
    Transpose3D%ZX = M1%XZ
    Transpose3D%ZY = M1%YZ
    Transpose3D%ZZ = M1%ZZ
  End Function Transpose3D

  ! Tensorial product
  Function TensPVect2D (V1, V2)
    Type (Vect2D), intent(IN)                   :: V1
    Type (Vect2D), intent(IN)                   :: V2
    Type (Mat2D)                                :: TensPVect2D

    TensPVect2D%XX = V1%X * V2%X
    TensPVect2D%XY = V1%X * V2%Y
    TensPVect2D%YX = V1%Y * V2%X
    TensPVect2D%YY = V1%Y * V2%Y
  End Function TensPVect2D

  Function TensPVect3D (V1, V2)
    Type (Vect3D), intent(IN)                   :: V1, V2
    Type (Mat3D)                                :: TensPVect3D

    TensPVect3D%XX = V1%X * V2%X 
    TensPVect3D%XY = V1%X * V2%Y
    TensPVect3D%XZ = V1%X * V2%Z
    TensPVect3D%YX = V1%Y * V2%X 
    TensPVect3D%YY = V1%Y * V2%Y
    TensPVect3D%YZ = V1%Y * V2%Z
    TensPVect3D%ZX = V1%Z * V2%X 
    TensPVect3D%ZY = V1%Z * V2%Y
    TensPVect3D%ZZ = V1%Z * V2%Z
  End Function TensPVect3D

  Function TensPMat2D (M1, M2)
    Type (Mat2D), Intent(IN)                    :: M1, M2
    Type (Mat2D)                                :: TensPMat2D

    TensPMat2D%XX = M1%XX * M2%XX
    TensPMat2D%XY = M1%XY * M2%XY
    TensPMat2D%YX = M1%YX * M2%YX
    TensPMat2D%YY = M1%YY * M2%YY
  End Function TensPMat2D

  Function TensPMatS2D (M1, M2)
    Type (MatS2D), Intent(IN)                   :: M1, M2
    Type (MatS2D)                               :: TensPMatS2D

    TensPMatS2D%XX = M1%XX * M2%XX
    TensPMatS2D%YY = M1%YY * M2%YY
    TensPMatS2D%XY = M1%XY * M2%XY
  End Function TensPMatS2D

  Function TensPMat3D (M1, M2)
    Type (Mat3D), Intent(IN)                    :: M1, M2
    Type (Mat3D)                                :: TensPMat3D

    TensPMat3D%XX = M1%XX * M2%XX
    TensPMat3D%XY = M1%XY * M2%XY
    TensPMat3D%XZ = M1%XZ * M2%XZ
    TensPMat3D%YX = M1%YX * M2%YX
    TensPMat3D%YY = M1%YY * M2%YY
    TensPMat3D%YZ = M1%YZ * M2%YZ
    TensPMat3D%ZX = M1%ZX * M2%ZX
    TensPMat3D%ZY = M1%ZY * M2%ZY
    TensPMat3D%ZZ = M1%ZZ * M2%ZZ
  End Function TensPMat3D

  Function TensPMatS3D (M1, M2)
    Type (MatS3D), Intent(IN)                   :: M1, M2
    Type (MatS3D)                               :: TensPMatS3D

    TensPMatS3D%XX = M1%XX * M2%XX
    TensPMatS3D%YY = M1%YY * M2%YY
    TensPMatS3D%ZZ = M1%ZZ * M2%ZZ
    TensPMatS3D%YZ = M1%YZ * M2%YZ
    TensPMatS3D%XZ = M1%XZ * M2%XZ
    TensPMatS3D%XY = M1%XY * M2%XY
  End Function TensPMatS3D

  ! Symmetrized product
  Function SymPVect2D (V1, V2)
    Type (Vect2D), intent(IN)                   :: V1, V2
    Type (MatS2D)                               :: SymPVect2D

    SymPVect2D = Symmetrize(V1 .TensP. V2)
  End Function SymPVect2D

  Function SymPVect3D (V1, V2)
    Type (Vect3D), intent(IN)                   :: V1, V2
    Type (MatS3D)                               :: SymPVect3D

    SymPVect3D = Symmetrize(V1 .TensP. V2)
  End Function SymPVect3D

  Function SymPMat2D(M1, M2)
    Type (Mat2D), Intent(IN)                    :: M1, M2
    Type (Mat2D)                                :: SymPMat2D

    SymPMat2D = ((M1 .TensP. M2) + (M2 .TensP. M1)) * 0.5_Kr
  End Function SymPMat2D

  Function SymPMatS2D(M1, M2)
    Type (MatS2D), Intent(IN)                   :: M1, M2
    Type (MatS2D)                               :: SymPMatS2D

    SymPMatS2D = M1 .TensP. M2
  End Function SymPMatS2D

  Function SymPMat3D(M1, M2)
    Type (Mat3D), Intent(IN)                    :: M1, M2
    Type (Mat3D)                                :: SymPMat3D

    SymPMat3D = ((M1 .TensP. M2) + (M2 .TensP. M1)) * 0.5_Kr
  End Function SymPMat3D

  Function SymPMatS3D(M1, M2)
    Type (MatS3D), Intent(IN)                   :: M1, M2
    Type (MatS3D)                               :: SymPMatS3D

    SymPMatS3D = M1 .TensP. M2
  End Function SymPMatS3D

  Function Trace2D(M1)
    Type (Mat2D), Intent(IN)                      :: M1
    PetscReal                                     :: Trace2D

    Trace2D = M1%XX + M1%YY
  End Function Trace2D

  Function Trace2DS(M1)
    Type (MatS2D), Intent(IN)                     :: M1
    PetscReal                                     :: Trace2DS

    Trace2DS = M1%XX + M1%YY
  End Function Trace2DS

  Function Trace3D(M1)
    Type (Mat3D), Intent(IN)                      :: M1
    PetscReal                                     :: Trace3D

    Trace3D = M1%XX + M1%YY + M1%ZZ
  End Function Trace3D

  Function Trace3DS(M1)
    Type (MatS3D), Intent(IN)                     :: M1
    PetscReal                                     :: Trace3DS

    Trace3DS = M1%XX + M1%YY + M1%ZZ
  End Function Trace3DS

  Subroutine Vect2D_Get_Real(V1, R1)
    Type (Vect2D), intent(OUT)                    :: V1
    PetscReal,       Intent(IN)                   :: R1

    V1%X = R1
    V1%Y = R1
  End Subroutine Vect2D_Get_Real

  Subroutine Vect3D_Get_Real(V1, R1)
    Type (Vect3D), intent(OUT)                    :: V1
    PetscReal,       Intent(IN)                   :: R1

    V1%X = R1
    V1%Y = R1
    V1%Z = R1
  End Subroutine Vect3D_Get_Real

  Subroutine Vect2D_Get_VectR(V1, R1)
    Type (Vect2D), intent(OUT)                    :: V1
    PetscReal,       Dimension(2), Intent(IN)     :: R1

    V1%X = R1(1)
    V1%Y = R1(2)
  End Subroutine Vect2D_Get_VectR

  Subroutine Vect3D_Get_VectR(V1, R1)
    Type (Vect3D), intent(OUT)                    :: V1
    PetscReal,       Dimension(3), Intent(IN)     :: R1

    V1%X = R1(1)
    V1%Y = R1(2)
    V1%Z = R1(3)
  End Subroutine Vect3D_Get_VectR

  Subroutine Vect2DEQ(V1, V2)
    Type (Vect2D), intent(OUT)                    :: V1
    Type (Vect2D), Intent(IN)                     :: V2

    V1%X = V2%X
    V1%Y = V2%Y
  End Subroutine Vect2DEQ

  Subroutine Vect3DEQ(V1, V2)
    Type (Vect3D), intent(OUT)                    :: V1
    Type (Vect3D), Intent(IN)                     :: V2

    V1%X = V2%X
    V1%Y = V2%Y
    V1%Z = V2%Z
  End Subroutine Vect3DEQ

  Subroutine Mat2D_Get_Real(M1, R1)
    Type (Mat2D), intent(OUT)                     :: M1
    PetscReal,       Intent(IN)                   :: R1

    M1%XX = R1 ; M1%XY = R1
    M1%YX = R1 ; M1%YY = R1
  End Subroutine Mat2D_Get_Real

  Subroutine Mat3D_Get_Real(M1, R1)
    Type (Mat3D), intent(OUT)                     :: M1
    PetscReal,       Intent(IN)                   :: R1

    M1%XX = R1 ; M1%XY = R1 ; M1%XZ = R1
    M1%YX = R1 ; M1%YY = R1 ; M1%YZ = R1
    M1%ZX = R1 ; M1%ZY = R1 ; M1%ZZ = R1
  End Subroutine Mat3D_Get_Real

  Subroutine Mat2DEQ(M1, M2)
    Type (Mat2D), Intent(OUT)                     :: M1
    Type (Mat2D), Intent(IN)                      :: M2

    M1%XX = M2%XX  ;  M1%XY = M2%XY
    M1%YX = M2%YX  ;  M1%YY = M2%YY
  End Subroutine Mat2DEQ

  Subroutine Mat3DEQ(M1, M2)
    Type (Mat3D), Intent(OUT)                     :: M1
    Type (Mat3D), Intent(IN)                      :: M2

    M1%XX = M2%XX  ;  M1%XY = M2%XY  ;  M1%XZ = M2%XZ
    M1%YX = M2%YX  ;  M1%YY = M2%YY  ;  M1%YZ = M2%YZ
    M1%ZX = M2%ZX  ;  M1%ZY = M2%ZY  ;  M1%ZZ = M2%ZZ
  End Subroutine Mat3DEQ

  Subroutine MatS2D_Get_Real(M1, R1)
    Type (MatS2D), intent(OUT)                    :: M1
    PetscReal,       Intent(IN)                   :: R1

    M1%XX = R1 ; M1%YY = R1 ; M1%XY = R1 
  End Subroutine MatS2D_Get_Real

  Subroutine MatS3D_Get_Real(M1, R1)
    Type (MatS3D), intent(OUT)                     :: M1
    PetscReal,       Intent(IN)                    :: R1

    M1%XX = R1 ; M1%YY = R1 ; M1%ZZ = R1
    M1%YZ = R1 ; M1%XZ = R1 ; M1%XY = R1
  End Subroutine MatS3D_Get_Real

  Subroutine MatS2DEQ(M1, M2)
    Type (MatS2D), Intent(OUT)                     :: M1
    Type (MatS2D), Intent(IN)                      :: M2

    M1%XX = M2%XX  ;  M1%YY = M2%YY ; M1%XY = M2%XY
  End Subroutine MatS2DEQ

  Subroutine MatS3DEQ(M1, M2)
    Type (MatS3D), Intent(OUT)                     :: M1
    Type (MatS3D), Intent(IN)                      :: M2

    M1%XX = M2%XX  ;  M1%YY = M2%YY  ;  M1%ZZ = M2%ZZ
    M1%YZ = M2%YZ  ;  M1%XZ = M2%XZ  ;  M1%XY = M2%XY
  End Subroutine MatS3DEQ

  Subroutine Tens4OS2DEQ(T1, T2)
    Type(Tens4OS2D), Intent(OUT)                        :: T1
    Type(Tens4OS2D), Intent(IN)                         :: T2

    T1%XXXX = T2%XXXX
    T1%XXXY = T2%XXXY
    T1%XXYY = T2%XXYY

    T1%XYXY = T2%XYXY
    T1%XYYY = T2%XYYY

    T1%YYYY = T2%YYYY
  End Subroutine Tens4OS2DEQ

  Subroutine Tens4OS3DEQ(T1, T2)
    Type(Tens4OS3D), Intent(OUT)                        :: T1
    Type(Tens4OS3D), Intent(IN)                         :: T2

    T1%XXXX = T2%XXXX  
    T1%XXXY = T2%XXXY
    T1%XXXZ = T2%XXXZ  
    T1%XXYY = T2%XXYY  
    T1%XXYZ = T2%XXYZ  
    T1%XXZZ = T2%XXZZ  
          
    T1%XYXY = T2%XYXY  
    T1%XYXZ = T2%XYXZ
    T1%XYYY = T2%XYYY
    T1%XYYZ = T2%XYYZ
    T1%XYZZ = T2%XYZZ
          
    T1%XZXZ = T2%XZXZ  
    T1%XZYY = T2%XZYY
    T1%XZYZ = T2%XZYZ
    T1%XZZZ = T2%XZZZ
          
    T1%YYYY = T2%YYYY
    T1%YYYZ = T2%YYYZ
    T1%YYZZ = T2%YYZZ
          
    T1%YZYZ = T2%YZYZ  
    T1%YZZZ = T2%YZZZ
          
    T1%ZZZZ = T2%ZZZZ 
  End Subroutine Tens4OS3DEQ

  Subroutine Tens4OS2D_Get_Real(T1, D1)
    Type(Tens4OS2D), Intent(OUT)                        :: T1
    PetscReal,       Intent(IN)                         :: D1

    T1%XXXX = D1
    T1%XXXY = D1
    T1%XXYY = D1

    T1%XYXY = D1
    T1%XYYY = D1

    T1%YYYY = D1
  End Subroutine Tens4OS2D_Get_Real

  Subroutine Tens4OS3D_Get_Real(T1, D1)
    Type(Tens4OS3D), Intent(OUT)                        :: T1
    PetscReal,       Intent(IN)                         :: D1

    T1%XXXX = D1  
    T1%XXXY = D1
    T1%XXXZ = D1  
    T1%XXYY = D1  
    T1%XXYZ = D1  
    T1%XXZZ = D1  
          
    T1%XYXY = D1  
    T1%XYXZ = D1  
    T1%XYYY = D1  
    T1%XYYZ = D1  
    T1%XYZZ = D1  
          
    T1%XZXZ = D1  
    T1%XZYY = D1  
    T1%XZYZ = D1  
    T1%XZZZ = D1  
          
    T1%YYYY = D1  
    T1%YYYZ = D1  
    T1%YYZZ = D1  
          
    T1%YZYZ = D1  
    T1%YZZZ = D1  
          
    T1%ZZZZ = D1  
  End Subroutine Tens4OS3D_Get_Real

  Function Symmetrize2D(M1)
    Type(Mat2D), Intent(IN)                             :: M1
    Type(MatS2D)                                        :: Symmetrize2D

    Symmetrize2D%XX = M1%XX
    Symmetrize2D%YY = M1%YY
    Symmetrize2D%XY = (M1%XY + M1%YX) * 0.5_Kr
  End Function Symmetrize2D

  Function Symmetrize3D(M1)
    Type(Mat3D), Intent(IN)                             :: M1
    Type(MatS3D)                                        :: Symmetrize3D

    Symmetrize3D%XX = M1%XX
    Symmetrize3D%YY = M1%YY
    Symmetrize3D%ZZ = M1%YY
    
    Symmetrize3D%YZ = (M1%YZ + M1%ZY) * 0.5_Kr
    Symmetrize3D%XZ = (M1%XZ + M1%ZX) * 0.5_Kr
    Symmetrize3D%XY = (M1%XY + M1%YX) * 0.5_Kr
  End Function Symmetrize3D
!====================================================================
!             END OF OPERATOR OVERLOADING
!====================================================================

  Function Vol_Tetra_3D(V1, V2, V3, V4)
    Type (Vect3D), Intent(IN)                   :: V1, V2, V3, V4
    PetscReal                                   :: Vol_Tetra_3D

    Type(Vect3D)                                :: C1, C2, C3

    C1 = V1 - V4
    C2 = V2 - V4
    C3 = V3 - V4
    Vol_Tetra_3D = ABS(C1%X * (C2%Y * C3%Z - C2%Z * C3%Y) -                   &
         &             C1%Y * (C2%X * C3%Z - C2%Z * C3%X) +                   &
         &             C1%Z * (C2%X * C3%Y - C2%Y * C3%X) ) / 6.0_Kr
  End Function Vol_Tetra_3D

  Function Area_Tri_2D(S1, S2, S3)
    Type (Vect2D), Intent(IN)                   :: S1, S2, S3
    PetscReal                                   :: Area_Tri_2D

    Type(Vect2D)                                :: C1, C2

    C1 = S2-S1
    C2 = S3-S1
    Area_Tri_2D = ABS(C1%X * C2%Y - C1%Y * C2%X) * 0.5_Kr
  End Function Area_Tri_2D

  Function Ht_Min_Tri_2D(S1, S2, S3)
    Type (Vect2D), Intent(IN)                      :: S1, S2, S3
    PetscReal                                      :: Ht_Min_Tri_2D

    Type (Vect2D)                                  :: C1, C2, C3
    PetscReal                                      :: H1, H2, H3, AreaX2

    C1 = S2-S3
    C2 = S3-S1
    C3 = S1-S2

    AreaX2 = abs(C1%X * C2%Y - C1%Y * C2%X)

    H1 = AreaX2 / SQRT( (C1 .DotP. C1) )
    H2 = AreaX2 / SQRT( (C2 .DotP. C2) )
    H3 = AreaX2 / SQRT( (C3 .DotP. C3) )

    Ht_Min_Tri_2D = Min (H1, H2, H3)
  End Function Ht_Min_Tri_2D

  Function ValP2D(M)
    Type (Mat2D), Intent(IN)                      :: M
    Type (Vect2D)                                 :: ValP2D, VPTemp
    PetscReal                                     :: Tmp1

    Tmp1 = (M%XX - M%YY)**2+ 4.0_Kr * M%XY * M%YX
    If (Tmp1 < 0.0_Kr) Then
       Print*, 'Error in Valp2D : non diagonalizable matrix'
       Print*, 'Tmp1 = ',Tmp1
       Print*, 'M :'
       Print*, M%XX, M%XY
       Print*, M%YX, M%YY
       Print*, '(M%XX - M%YY)**2',  (M%XX - M%YY)**2
       Print*, '4.0_Kr * M%XY * M%YX', 4.0_Kr * M%XY * M%YX
       Print*, 'Result is (0.0, 0.0)'
       VPTemp%X = 0.0_Kr
       VPTemp%Y = 0.0_Kr
    Else
       VPTemp%X = (M%XX + M%YY  + SQRT ( Tmp1)) * 0.5_Kr
       VPTemp%Y = (M%XX + M%YY  - SQRT ( Tmp1)) * 0.5_Kr
    EndIf
    ValP2D = VPTemp
  End Function ValP2D

  Function ValP2DS(M)
    Type (MatS2D), Intent(IN)                     :: M
    Type (Vect2D)                                 :: ValP2DS, VPTemp
    PetscReal                                     :: Tmp1

    Tmp1 = (M%XX - M%YY)**2+ 4.0_Kr * M%XY**2

    VPTemp%X = (M%XX + M%YY  + SQRT ( Tmp1)) * 0.5_Kr
    VPTemp%Y = (M%XX + M%YY  - SQRT ( Tmp1)) * 0.5_Kr

    ValP2DS= VPTemp
  End Function ValP2DS

End Module m_MEF_LinAlg
