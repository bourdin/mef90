Program ModEner

   Implicit NONE
  
   Character(len=128)                                 :: Prefix, Ener_Str, Mod_Str
   Integer, Parameter                                 :: Mod_Unit = 98
   Integer, Parameter                                 :: Ener_Unit = 99
   Real                                               :: junk
   Real, Dimension(:), Pointer                        :: Load
   Real, Dimension(:), Pointer                        :: Bulk_Ener
   Real, Dimension(:), Pointer                        :: Surf_Ener
   Real, Dimension(:), Pointer                        :: Tot_Ener
   Integer                                            :: Flg, iErr, iTS, TimeStep



   Write(*,100, advance = 'no') 'Prefix: '
   Read(*,100) Prefix
   
   Ener_Str = Trim(Prefix)//'.ener'
   Mod_Str  = Trim(Prefix)//'.enermod'
   Open (File = Ener_Str, Unit = Ener_Unit, status = 'Old',               &
         & Action = 'Read')
   Rewind(Ener_Unit)
   Do
      Read(Ener_Unit, *, end=50, err=50) iTS, junk, junk, junk, junk
      CYCLE
50    EXIT
   End Do

   Print*, 'Final Time step is ', iTS
   Allocate (Load(iTS))
   Allocate (Bulk_Ener(iTS))
   Allocate (Surf_Ener(iTS))
   Allocate ( Tot_Ener(iTS))

   Rewind(Ener_Unit)
   Do
      Read(Ener_Unit, *, end=60, err=60) TimeStep, Load(TimeStep), Bulk_Ener(TimeStep), &
         & Surf_Ener(TimeStep),  Tot_Ener(TimeStep)
      CYCLE
60    EXIT
   End Do
   Close(Ener_Unit)
   Print*, 'Done reading '//Ener_Str
  
   Open (File = Mod_Str, Unit = Mod_Unit, status = 'Unknown')
   Rewind(Mod_Unit)
   Do TimeStep = 1, iTS
      Write(Mod_Unit, 70) TimeStep, Load(TimeStep), Bulk_Ener(TimeStep), &
         & Surf_Ener(TimeStep),  Tot_Ener(TimeStep)
   End Do
   Close(Mod_Unit)
   Print*, 'Done writing '//Mod_Str
   
   DeAllocate(Load, Bulk_Ener, Surf_Ener, Tot_Ener)
      
70 Format(I4, 4(ES13.5,'  '))
100 Format(A)      
End Program ModEner  