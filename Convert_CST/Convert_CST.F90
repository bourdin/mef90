Program Convert_CST
   Use m_MEF90
   Use m_Rupt_Struct
   Implicit NONE
   
   Integer                                        :: iDim
   Type(Tens4OS2D)                                :: HL2D
   Type(Tens4OS3D)                                :: HL3D
   Real(Kind = Kr)                                :: E, nu, Gc, alpha
   Integer                                        :: iBlk, NumBlk, iJunk
   Character(len = 256)                           :: CST_IN, CST_OUT2D, CST_OUT3D
   
   Write(*,100, advance = 'no') 'Prefix: '
   Read(*,100) CST_IN
   
   CST_OUT2D = Trim(CST_IN) // '.CST2D'
   CST_OUT3D = Trim(CST_IN) // '.CST3D'
   CST_IN    = Trim(CST_IN)//'.CST'
   
   Print*, 'CST_IN', CST_IN
   Print*, 'CST_OUT2D', CST_OUT2D
   Print*, 'CST_OUT3D', CST_OUT3D

   Open(File = CST_IN, Unit = 90, Status = 'Old')
   Rewind(90)
   Open(File = CST_OUT2D, Unit = 91, Status = 'Unknown')
   Rewind(91)
   Open(File = CST_OUT3D, Unit = 92, Status = 'Unknown')
   Rewind(92)
   
   Read(90,*) NumBlk

   Write(91,110) NumBlk
   Write(92,130) NumBlk
   
   Do iBlk = 1, NumBlk
      Read(90, *) ijunk, Gc, E, nu, alpha
      Call GenHL_Iso2D_EnuPlaneStress(E, nu, HL2D)
      Write(91,120) iBlk, Gc, HL2D, Alpha
      Call GenHL_Iso3D_Enu(E, nu, HL3D)
      Write(92,140) iBlk, Gc, HL3D, Alpha
   End Do
   Close(90)
   Close(91)
   Close(92)
  
  
100 Format(A)
110 Format(I6,' Toughness    A1111        A1112        A1122        A1212        A1222        A2222        Alpha')
120 Format(I6, 8(ES12.5,' '))
130 Format(I6,' Toughness    A_1111       A_1112       A_1113       A_1122       A_1123       A_1133       A_1212       A_1213       A_1222       A_1223       A_12133      A_1313       A_1322       A_1323       A_1333       A_2222       A_2223       A_2233       A_2323       A_2333       A_3333       Alpha')
140 Format(I6, 23(ES12.5,' '))

End Program Convert_CST