Module m_constantes

   IMPLICIT NONE
   PRIVATE
#include "finclude/petsc.h"

   
   ! Get the numerical value for differents kinds.
   ! Do _NOT_ change this but the parameters Ki, Ks and Kr
   !
   
   Integer, Parameter    :: spr = Selected_Real_Kind(1,1)
   ! spr is the smallest real kind i.e. single precision
   
   Integer, Parameter    :: dpr = Selected_Real_Kind(2*Precision(1.0_spr))
   ! dpr is the F77 called double precision real if exists
   
   Integer, Parameter    :: qpr = Selected_Real_Kind(2*Precision(1.0_dpr))
   ! qpr is the F77 called quad precision (?) real if exists, double otherwise
   
   Integer, Parameter    :: spi = Selected_Int_Kind(1)
   ! spi is the smallest integer kind
   
   Integer, Parameter    :: dpi = Selected_Int_Kind(Range(1_spi)+1)
   ! dpi is usually the standard F77 integer
   
   Integer, Parameter    :: qpi = Selected_Int_Kind(Range(1_dpi)+1)
   
   
   
   ! Now set the default precision for the MEF90 Modules
   ! This shouldn't be changed
   ! This way, the kinds are not hard coded and adapt automatically
   ! to any machine ...
   ! thank's to Michael Metcalf in comp.lang.fortran
   
   
   ! The following ensures that mef90 and PETSC real types are compatible:
   PetscReal, Parameter                  :: PReal = 1.0
   Integer, Parameter, Public            :: Kr = Selected_Real_Kind(Precision(PReal))
                                         
   PetscInt, Parameter                   :: PInt = 1
   Integer, Parameter, Public            :: Ki = Selected_Int_Kind(PInt)
                                         
   Integer, Parameter, Public            :: KiS = spi
   
   
   Real(Kind = Kr), Parameter, Public    :: Dble0  = 0._Kr
   Real(Kind = Kr), Parameter, Public    :: Dble1  = 1.0_Kr
   Real(Kind = Kr), Parameter, Public    :: DbleM1 = -1.0_Kr
   Real(Kind = Kr), Parameter, Public    :: Dble2  = 2.0_Kr
   Real(Kind = Kr), Parameter, Public    :: DbleM2 = -2.0_Kr
   !Real(Kind = Kr), Parameter, Public    :: CL_Pen = 10.0D+30
   
   Real(Kind = Kr), Parameter, Public   :: Pi     = 3.14159265358979323846264_Kr
   
   Real(Kind = Kr), Parameter, Public   :: DbleE  = 2.71828182845904523536029_Kr
   
   Real(Kind = Kr), Parameter, Public   :: InvOf2 = 1.0_Kr / 2.0_Kr
   Real(Kind = Kr), Parameter, Public   :: InvOf3 = 1.0_Kr / 3.0_Kr
   Real(Kind = Kr), Parameter, Public   :: InvOf4 = 1.0_Kr / 4.0_Kr
   Real(Kind = Kr), Parameter, Public   :: InvOf5 = 1.0_Kr / 5.0_Kr
   Real(Kind = Kr), Parameter, Public   :: InvOf6 = 1.0_Kr / 6.0_Kr
   Real(Kind = Kr), Parameter, Public   :: InvOf7 = 1.0_Kr / 7.0_Kr
   Real(Kind = Kr), Parameter, Public   :: InvOf8 = 1.0_Kr / 8.0_Kr
   Real(Kind = Kr), Parameter, Public   :: InvOf9 = 1.0_Kr / 9.0_Kr
   
   Integer, Parameter, Public           :: F_In   = 60
   ! Default unit used to read in files
   
   Integer, Parameter, Public           :: F_Out  = 61
   ! Default unit to write files
   
   Integer, Parameter, Public           :: Numbering_PerNodes = 13
   Integer, Parameter, Public           :: Numbering_PerCoord = 11
   
   Integer, Public                      :: MEF90_GaussOrder
   Real(Kind = Kr), Parameter, Public   :: MEF90_VLV = 1.0E+30
   Integer, Public                      :: MEF90_MyRank
   Integer, Public                      :: MEF90_NumProcs
   

End Module m_Constantes
