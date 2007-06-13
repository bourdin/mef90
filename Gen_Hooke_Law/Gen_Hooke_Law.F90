Program Gen_Hooke_Law
   Use m_MEF90
   Use m_Rupt_Struct
   
   Implicit NONE
   
   
   
   Real(Kind = Kr)      :: lambda, mu
   Real(Kind = Kr)      :: mu1, mu2, Theta
   Real(Kind = Kr)      :: E, nu
   
   Integer              :: iDim
   Type (Tens4OS2D)     :: A2D
   Type (Tens4OS3D)     :: A3D
   Integer, Parameter   :: HL_Iso2DStrain   = 1
   Integer, Parameter   :: HL_Iso2DStress   = 2
   Integer, Parameter   :: HL_Iso2DlambdaMu = 3
   Integer, Parameter   :: HL_Ortho2D       = 4
   Integer, Parameter   :: HL_Iso3DEnu      = 5
   Integer, Parameter   :: HL_Iso3DLambdaMu = 6
   Integer              :: HL_Type
   
   
   Write(*, 110) HL_Iso2DStrain,    'Isotropic 2D plain strain (E, nu)'   
   Write(*, 110) HL_Iso2DStress,    'Isotropic 2D plain stress (E, nu)'   
   Write(*, 110) HL_Iso2Dlambdamu,  'Isotropic 2D (lambda, mu)'   
   Write(*, 110) HL_Ortho2D    ,    'Orthotropic 2D (lambda, mu1, mu2, theta)'
   Write(*, 110) HL_Iso3DEnu      , 'Isotropic 3D (E, nu)'
   Write(*, 110) HL_Iso3DLambdaMu , 'Isotropic 3D (lambda, mu)'
   Write(*, 100, advance = 'no')    'Type of law: ' ; Read(*, *) HL_Type

   Select Case (HL_Type)
   Case(HL_Iso2DStrain)
      Write(*, 100, advance = 'no') '   E:      '; Read(*, *) E
      Write(*, 100, advance = 'no') '   nu:     '; Read(*, *) nu
      lambda = 1.
      mu     = 1.
      A2D = iso2D(Lambda, mu)
      Write(*, *)
      Write(*,130)
      Write(*,120) A2D
   Case(HL_Iso2DStress)
      Write(*, 100, advance = 'no') '   E:      '; Read(*, *) E
      Write(*, 100, advance = 'no') '   nu:     '; Read(*, *) nu
      lambda = E * nu / (1.0_Kr - nu**2)
      mu     = .5_Kr / (1.0_Kr+ nu)
      A2D = iso2D(Lambda, mu)
      Write(*, *)
      Write(*,130)
      Write(*,120) A2D
   Case(HL_Iso2Dlambdamu)
      Write(*, 100, advance = 'no') '   lambda: '; Read(*, *) lambda
      Write(*, 100, advance = 'no') '   mu:     '; Read(*, *) mu
      A2D = iso2D(Lambda, mu)
      Write(*, *)
      Write(*,130)
      Write(*,120) A2D
   Case(HL_Ortho2D)
      Write(*, 100, advance = 'no') '   lambda: '; Read(*, *) lambda
      Write(*, 100, advance = 'no') '   mu1:    '; Read(*, *) mu1
      Write(*, 100, advance = 'no') '   mu2:    '; Read(*, *) mu2
      Write(*, 100, advance = 'no') '   theta:  '; Read(*, *) theta
      theta = theta / 180_Kr * Pi
      A2D = Ortho2D(Lambda, mu1, mu2, theta)
      Write(*, *)
      Write(*,130)
      Write(*,120) A2D
   Case Default
      Write(*, 100) 'Wrong choice or not implemented yet!'
   End Select
   
100 Format(A)
110 Format('   [',I1,'] ',A)
120 Format(6(ES12.5,' '))
130 Format(' A_1111       A_1112       A_1122       A_1212       A_1222       A_2222')
Contains

   Function Iso2D(lambda, mu)
      Real(Kind = Kr), Intent(IN)         :: lambda, mu
      Type(Tens4OS2D)                     :: Iso2D
      Iso2D = 0.0_Kr
      
      Iso2D%XXXX = lambda + 2.0_Kr * mu
      Iso2D%XXYY = lambda
      Iso2D%XYXY = mu
      Iso2D%YYYY = lambda + 2.0_Kr * mu
   End Function Iso2D         
   
   Function Ortho2D(lambda, mu1, mu2, theta)
      Real(Kind = Kr), Intent(IN)         :: Lambda, mu1, mu2, theta
      Type(Tens4OS2D)                     :: Ortho2D
      
      Ortho2D = 0.0_Kr
      Ortho2D%XXXX = lambda + mu1 * .5_Kr + mu2 * .5_Kr+ mu1 * (cos(theta))**2 + mu2 * (sin(theta))**2
      Ortho2D%XXXY = (mu1-mu2) * cos(theta) * sin(theta)
      Ortho2D%XXYY = lambda + (mu1-mu2) * (.5_Kr - (cos(theta))**2)

      Ortho2D%XYXY = mu1 * (sin(theta))**2 + mu2 * (cos(theta))**2
      Ortho2D%XYYY = -Ortho2D%XXXY

      Ortho2D%YYYY = Ortho2D%XXXX
   End Function Ortho2D
   
End Program Gen_Hooke_Law   
