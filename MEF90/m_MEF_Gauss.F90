Module m_MEF_Gauss
  ! Blaise Bourdin, 1996-1998
  ! Merci de faire parvenir toutes remarques et bugs 
  ! eventuels aux adresses suivantes :
  ! bourdin@lpmtm.univ-paris13.fr
  ! bourdin@mat.dtu.dk
  !
  ! En cas de modifications, les noms de  modules et de 
  ! fichiers _DOIVENT_ etre renommes
  !          ^^^^^^^^^
  !
  Use m_MEF_Types
  Use m_Utils
  
  IMPLICIT NONE
  
  ! Allocation des composants NB_Gauss, BF et Der_BF de Elem_db
  !
  ! Syntaxe InitGauss_XXXX(Elem_bd, Nodes_bd, Order) avec
  !         XXXX = Nom de l'element
  !         Order = Order d'integration du pb EF
  
  ! DANS TOUS CE QUI SUIT, ON SUPPOSE QUE LA TRANSFO ELEMT COURANT ->
  ! ELEMT DE REF. EST AFFINE (Elt diff = Area de l'elt courant)
  
  ! Actuellement :
  !     1D_P1           Order 3
  !     2DA_Tri_P1      Order 2
  
  ! Blaise Bourdin 04-97
  ! Version 3   09-97
  
  Interface Init_Gauss_P1
     Module Procedure InitGauss_1D_P1, InitGauss_2D_Elast_P1,                 &
          & InitGauss_2D_Scal_P1
  End Interface

  Interface Init_Gauss_EXO
     Module Procedure Init_Gauss_EXO_2D_Scal, Init_Gauss_EXO_2D,              &
          & Init_Gauss_EXO_2D_Elast, Init_Gauss_EXO_3D_Scal,                  &
          & Init_Gauss_EXO_3D, Init_Gauss_EXO_3D_Elast
  End Interface

  Interface Destroy_Gauss_EXO
     Module Procedure Destroy_Gauss_EXO_3D, Destroy_Gauss_EXO_3D_Scal,        &
          & Destroy_Gauss_EXO_3D_Elast, Destroy_Gauss_EXO_2D,                 &
          & Destroy_Gauss_EXO_2D_Scal, Destroy_Gauss_EXO_2D_Elast
  End Interface
Contains
  
  
  Subroutine InitGauss_1D_P1(Elem_bd, Nodes_bd, Order)
    Type (Element1D), Dimension(:), Pointer     :: Elem_bd
    Type (Node1D), Dimension(:), Pointer        :: Nodes_bd
    Integer                                     :: Order
    
    Integer                                     :: NE, iE
    Integer                                     :: Nb_Gauss, NB_BF
    Real(Kind = Kr)                             :: Area, InvOfArea
    
    ! On verifie d'abord que les bd sont allouees :
    If ((.NOT. Associated(Elem_bd)) .OR. (.NOT. Associated(Nodes_bd))) Then
       Print*, 'One of the required database is not allocated.'
       Stop
    EndIf
    
    NE = Size(Elem_bd)
    
    
    OrderG : Select Case(Order)
    Case (1)
       Nb_Gauss = 2
       Nb_BF = 2
       DoiE1: Do iE = 1, NE
          Elem_bd(iE)%NB_Gauss = Nb_Gauss
          Allocate(Elem_bd(iE)%BF(Nb_BF,Nb_Gauss))
          Allocate(Elem_bd(iE)%Der_BF(NB_BF,Nb_Gauss))
          Elem_bd(iE)%BF(1,1) = 1.0_Kr
          Elem_bd(iE)%BF(1,2) = 0.0_Kr
          Elem_bd(iE)%BF(2,1) = 0.0_Kr
          Elem_bd(iE)%BF(2,2) = 1.0_Kr

          Area = ABS(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord)
          InvOfArea = 1.0_Kr / Area

          Elem_bd(iE)%Der_BF(1,1) = -1.0_Kr * InvOfArea
          Elem_bd(iE)%Der_BF(1,2) = -1.0_Kr * InvOfArea
          Elem_bd(iE)%Der_BF(2,1) = InvOfArea
          Elem_bd(iE)%Der_BF(2,2) = InvOfArea

          Allocate(Elem_bd(iE)%Gauss_C(Nb_Gauss))
          Elem_bd(iE)%Gauss_C(1) = Area * 0.5_Kr
          Elem_bd(iE)%Gauss_C(2) = Area * 0.5_Kr
       EndDo DoiE1
       
    Case (2)
       Nb_Gauss = 3
       NB_BF = 2
       DoiE2: Do iE = 1, NE
          Elem_bd(iE)%NB_Gauss = Nb_Gauss
          Allocate(Elem_bd(iE)%BF(Nb_BF,Nb_Gauss))
          Allocate(Elem_bd(iE)%Der_BF(NB_BF,Nb_Gauss))
          Elem_bd(iE)%BF(1,1) = 1.0_Kr
          Elem_bd(iE)%BF(1,2) = 0.5_Kr
          Elem_bd(iE)%BF(1,3) = 0.0_Kr

          Elem_bd(iE)%BF(2,1) = 0.0_Kr
          Elem_bd(iE)%BF(2,2) = 0.5_Kr
          Elem_bd(iE)%BF(2,3) = 1.0_Kr

          Area = ABS(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord)
          InvOfArea = 1.0_Kr / Area

          Elem_bd(iE)%Der_BF(1,1) = -1.0_Kr * InvOfArea
          Elem_bd(iE)%Der_BF(1,2) = -1.0_Kr * InvOfArea
          Elem_bd(iE)%Der_BF(1,3) = -1.0_Kr * InvOfArea

          Elem_bd(iE)%Der_BF(2,1) = InvOfArea
          Elem_bd(iE)%Der_BF(2,2) = InvOfArea
          Elem_bd(iE)%Der_BF(2,3) = InvOfArea

          Allocate(Elem_bd(iE)%Gauss_C(Nb_Gauss))

          Elem_bd(iE)%Gauss_C(1) = Area * InvOf6
          Elem_bd(iE)%Gauss_C(2) = 2.0_Kr * Area * InvOf3
          Elem_bd(iE)%Gauss_C(3) = Area * InvOf6

       EndDo DoiE2

    Case (3)
       Nb_Gauss = 4
       NB_BF = 2
       DoiE3: Do iE = 1, NE
          Elem_bd(iE)%NB_Gauss = Nb_Gauss
          Allocate(Elem_bd(iE)%BF(Nb_BF,Nb_Gauss))
          Allocate(Elem_bd(iE)%Der_BF(NB_BF,Nb_Gauss))
          Elem_bd(iE)%BF(1,1) = 1.0_Kr
          Elem_bd(iE)%BF(1,2) = 2.0_Kr * InvOf3
          Elem_bd(iE)%BF(1,3) = 1.0_Kr * InvOf3
          Elem_bd(iE)%BF(1,4) = 0.0_Kr

          Elem_bd(iE)%BF(2,1) = 0.0_Kr
          Elem_bd(iE)%BF(2,2) = 1.0_Kr * InvOf3
          Elem_bd(iE)%BF(2,3) = 2.0_Kr * InvOf3
          Elem_bd(iE)%BF(2,4) = 1.0_Kr

          Area = ABS(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord)

          InvOfArea = 1.0_Kr / Area

          Elem_bd(iE)%Der_BF(1,1) = -1.0_Kr * InvOfArea
          Elem_bd(iE)%Der_BF(1,2) = -1.0_Kr * InvOfArea
          Elem_bd(iE)%Der_BF(1,3) = -1.0_Kr * InvOfArea
          Elem_bd(iE)%Der_BF(1,4) = -1.0_Kr * InvOfArea

          Elem_bd(iE)%Der_BF(2,1) = InvOfArea
          Elem_bd(iE)%Der_BF(2,2) = InvOfArea
          Elem_bd(iE)%Der_BF(2,3) = InvOfArea
          Elem_bd(iE)%Der_BF(2,4) = InvOfArea

          Allocate(Elem_bd(iE)%Gauss_C(Nb_Gauss))

          Elem_bd(iE)%Gauss_C(1) = Area * InvOf8 
          Elem_bd(iE)%Gauss_C(2) = 3.0_Kr * Area * InvOf8
          Elem_bd(iE)%Gauss_C(3) = 3.0_Kr * Area * InvOf8
          Elem_bd(iE)%Gauss_C(4) = Area * InvOf8           
          
       EndDo DoiE3
       
    Case Default
       Print*, 'This order is not implemented yet...', Order
       Stop
       
    End Select OrderG
  End Subroutine InitGauss_1D_P1
  
  Subroutine Init_Gauss_EXO_2D_Scal(Elems_db, Nodes_db, Geom, Order, Elem,    &
       & Elem_List)
    Type (Element2D_Scal), Dimension(:), Pointer        :: Elems_db
    Type (Node2D), Dimension(:), Pointer                :: Nodes_db
    Type (EXO_Geom_Info), intent(IN)                    :: Geom
    Integer                                             :: Order
    Integer, Intent(IN), Optional                       :: Elem
    Integer, Dimension(:), Pointer, Optional            :: Elem_List
    
    Integer                                             :: iE, iBlk, eNum
    Integer                                             :: NB_BF, NB_Gauss
    Real(Kind = Kr)                                     :: Jac, Jac_Inv
    Integer, Dimension(:), Pointer                      :: Elements

    ! Make sure the databases are allocated
    If ((.NOT. Associated(Elems_db)) .OR. (.NOT. Associated(Nodes_db))) Then
       Print*, 'One of the required databases is not allocated...'
       Stop
    EndIf
!!$    If (Geom%Numbering /= Numbering_PerNodes) Then
!!$       Write(*,*) 'Error, incompatible numbering scheme for local DoF'
!!$    End If
!!$
    !!! Test args here, and allocate elem_list accordingly, if necessary
    If ( Present(Elem) .AND. Present(Elem_List) ) Then
       Print*, 'Error: Init_Gauss_EXO_2D_Scal, cannot specify Elem and',      &
            &  'Elem_List at the same time'
       STOP
    ElseIf ( Present(Elem) ) Then
       Allocate (Elements(1))
       Elements = Elem
    ElseIf ( .NOT. Present(Elem_List) ) Then
       Allocate(Elements(Geom%Num_Elems))
       Elements = (/ (iE, iE=1,Geom%Num_Elems) /)
    Else
       Allocate(Elements(Size(ELem_List)))
       Elements = Elem_List
    End If
    
    Do_eNum: Do eNum = 1, Size(Elements)
       iE = Elements(eNum)
       iBlk = Elems_db(iE)%ID_EL
       Select case (Trim(Geom%Elem_Blk(iBlk)%Type))
       Case('TRI3')
          Nb_BF = 3
          ! Base functions: X, Y, 1-X-Y
          Jac = Area_Tri_2D(Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord, &
               & Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord,&
               & Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord) * 2.0_Kr
          Jac_Inv = 1.0_Kr / Jac
          Select Case (Order)
          Case(1)
             Nb_Gauss = 1
             Elems_db(iE)%Nb_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(NB_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF, NB_Gauss))
             Allocate(Elems_db(iE)%Grad_BF(Nb_BF, NB_Gauss))
             Elems_db(iE)%Gauss_C(1) = Jac / 2.0_Kr
             Elems_db(iE)%BF(1,1) = InvOf3
             Elems_db(iE)%BF(2,1) = InvOf3
             Elems_db(iE)%BF(3,1) = InvOf3
             
          Case(2)
             Nb_Gauss = 3
             Elems_db(iE)%Nb_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(NB_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF, NB_Gauss))
             Allocate(Elems_db(iE)%Grad_BF(Nb_BF, NB_Gauss))
             Elems_db(iE)%Gauss_C = Jac * InvOf6

             Elems_db(iE)%BF(1,1) = InvOf6
             Elems_db(iE)%BF(1,2) = 2.0_Kr * InvOf3
             Elems_db(iE)%BF(1,3) = InvOf6
                          
             Elems_db(iE)%BF(2,1) = InvOf6
             Elems_db(iE)%BF(2,2) = InvOf6
             Elems_db(iE)%BF(2,3) = 2.0_Kr * InvOf3
                          
             Elems_db(iE)%BF(3,:) = 1.0_Kr - Elems_db(iE)%BF(1,:) -           &
                  & Elems_db(iE)%BF(2,:) 

          Case(3)
             Nb_Gauss = 4
             Elems_db(iE)%Nb_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(NB_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF, NB_Gauss))
             Allocate(Elems_db(iE)%Grad_BF(Nb_BF, NB_Gauss))
             Elems_db(iE)%Gauss_C = Jac * 25.0_Kr / 48.0_Kr / 2.0_Kr
             Elems_db(iE)%Gauss_C(1) = -Jac * 9.0_Kr / 16.0_Kr / 2.0_Kr

             Elems_db(iE)%BF(1,1) = InvOf3
             Elems_db(iE)%BF(1,2) = 3.0_Kr * InvOf5
             Elems_db(iE)%BF(1,3) = InvOf5
             Elems_db(iE)%BF(1,4) = InvOf5

             Elems_db(iE)%BF(2,1) = InvOf3
             Elems_db(iE)%BF(2,2) = InvOf5
             Elems_db(iE)%BF(2,3) = 3.0_Kr * InvOf5
             Elems_db(iE)%BF(2,4) = InvOf5

             Elems_db(iE)%BF(3,:) = 1.0_Kr - Elems_db(iE)%BF(1,:) -           &
                  & Elems_db(iE)%BF(2,:) 
          Case Default
             Write(*,*) 'Order not implemented for this element type ',       &
                  & Trim(Geom%Elem_Blk(iBlk)%Type), Order
             STOP
          End Select
       Case Default
          Write(*,*) 'Element type not implemented yet...',                   &
               & Geom%Elem_Blk(iBlk)%Type
          STOP
       End Select

       Elems_db(iE)%Grad_BF(:,:)%X = 0.0_Kr
       Elems_db(iE)%Grad_BF(:,:)%Y = 0.0_Kr

       Elems_db(iE)%Grad_BF(1,:)%X =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y) * Jac_Inv
       Elems_db(iE)%Grad_BF(1,:)%Y =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X) * Jac_Inv

       Elems_db(iE)%Grad_BF(2,:)%X =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y) * Jac_Inv
       Elems_db(iE)%Grad_BF(2,:)%Y =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X) * Jac_Inv

       Elems_db(iE)%Grad_BF(3,:)%X =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y) * Jac_Inv
       Elems_db(iE)%Grad_BF(3,:)%Y =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X) * Jac_Inv

    End Do Do_eNum
    DeAllocate(Elements)
  End Subroutine Init_Gauss_EXO_2D_Scal
    
  Subroutine Init_Gauss_EXO_2D(Elems_db, Nodes_db, Geom, Order, Elem,         &
       & Elem_List)
    Type (Element2D), Dimension(:), Pointer             :: Elems_db
    Type (Node2D), Dimension(:), Pointer                :: Nodes_db
    Type (EXO_Geom_Info), intent(IN)                    :: Geom
    Integer                                             :: Order
    Integer, Intent(IN), Optional                       :: Elem
    Integer, Dimension(:), Pointer, Optional            :: Elem_List
    
    Integer                                             :: iE, iBlk, eNum
    Integer                                             :: NB_BF, NB_Gauss
    Real(Kind = Kr)                                     :: Jac, Jac_Inv
    Integer, DImension(:), Pointer                      :: Elements

    ! Make sure the databases are allocated
    If ((.NOT. Associated(Elems_db)) .OR. (.NOT. Associated(Nodes_db))) Then
       Print*, 'One of the required databases is not allocated...'
       Stop
    EndIf

    !!! Test args here, and allocate elem_list accordingly, if necessary
    If ( Present(Elem) .AND. Present(Elem_List) ) Then
       Print*, 'Error: Init_Gauss_EXO_2D, cannot specify Elem and',           &
            &  'Elem_List at the same time'
       STOP
    ElseIf ( Present(Elem) ) Then
       Allocate (Elements(1))
       Elements = Elem
    ElseIf ( .NOT. Present(Elem_List) ) Then
       Allocate(Elements(Geom%Num_Elems))
       Elements = (/ (iE, iE=1,Geom%Num_Elems) /)
    Else
       Allocate(Elements(Size(ELem_List)))
       Elements = Elem_List
    End If
    
    Do_eNum: Do eNum = 1, Size(Elements)
       iE = Elements(eNum)
       iBlk = Elems_db(iE)%ID_EL
       Select case (Trim(Geom%Elem_Blk(iBlk)%Type))
       Case('TRI3')
          Nb_BF = 6
          ! Base functions: (X,0), (Y,0), (1-X-Y,0)
          !                 (0,X), (0,Y), (0,1-X-Y)  
          Jac = Area_Tri_2D(Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord,           &
               & Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord,                      &
               & Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord) * 2.0_Kr
          Jac_Inv = 1.0_Kr / Jac
          Select Case (Order)
          Case(1)
             Nb_Gauss = 1
             Elems_db(iE)%Nb_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(NB_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF, NB_Gauss))
             Allocate(Elems_db(iE)%Der_BF(Nb_BF, NB_Gauss)) 
             Elems_db(iE)%Gauss_C(1) = Jac / 2.0_Kr
             Elems_db(iE)%BF(1,1) = (/ InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(2,1) = (/ InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(3,1) = (/ InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(4,1) = (/ 0.0_Kr, InvOf3 /)
             Elems_db(iE)%BF(5,1) = (/ 0.0_Kr, InvOf3 /)
             Elems_db(iE)%BF(6,1) = (/ 0.0_Kr, InvOf3 /)
             
          Case(2)
             Nb_Gauss = 3
             Elems_db(iE)%Nb_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(NB_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF, NB_Gauss))
             Allocate(Elems_db(iE)%Der_BF(Nb_BF, NB_Gauss))
             Elems_db(iE)%Gauss_C = Jac * InvOf3 / 2.0_Kr

             Elems_db(iE)%BF(1,1) = (/ InvOf6, 0.0_Kr /)
             Elems_db(iE)%BF(1,2) = (/ 2.0_Kr * InvOf3, 0.0_Kr/)
             Elems_db(iE)%BF(1,3) = (/ InvOf6, 0.0_Kr/)
                          
             Elems_db(iE)%BF(2,1) = (/ InvOf6, 0.0_Kr /)
             Elems_db(iE)%BF(2,2) = (/ InvOf6, 0.0_Kr /)
             Elems_db(iE)%BF(2,3) = (/ 2.0_Kr * InvOf3, 0.0_Kr /)
                          
             Elems_db(iE)%BF(3,1) = (/ 2.0_Kr * InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(3,2) = (/ InvOf6, 0.0_Kr /)
             Elems_db(iE)%BF(3,3) = (/ InvOf6, 0.0_Kr /)

             Elems_db(iE)%BF(4,1) = (/ 0.0_Kr, InvOf6 /)
             Elems_db(iE)%BF(4,2) = (/ 0.0_Kr, 2.0_Kr * InvOf3 /)
             Elems_db(iE)%BF(4,3) = (/ 0.0_Kr, InvOf6 /)
                          
             Elems_db(iE)%BF(5,1) = (/ 0.0_Kr, InvOf6 /)
             Elems_db(iE)%BF(5,2) = (/ 0.0_Kr, InvOf6 /)
             Elems_db(iE)%BF(5,3) = (/ 0.0_Kr, 2.0_Kr * InvOf3 /)

             Elems_db(iE)%BF(6,1) = (/ 0.0_Kr, 2.0_Kr * InvOf3 /)
             Elems_db(iE)%BF(6,2) = (/ 0.0_Kr, InvOf6 /)
             Elems_db(iE)%BF(6,3) = (/ 0.0_Kr, InvOf6 /)
          Case(3)
             Nb_Gauss = 4
             Elems_db(iE)%Nb_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(NB_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF, NB_Gauss))
             Allocate(Elems_db(iE)%Der_BF(Nb_BF, NB_Gauss))
             Elems_db(iE)%Gauss_C = Jac * 25.0_Kr / 48.0_Kr / 2.0_Kr
             Elems_db(iE)%Gauss_C(1) = -Jac * 9.0_Kr / 16.0_Kr / 2.0_Kr

             Elems_db(iE)%BF(1,1) = (/ InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(1,2) = (/ 3.0_Kr * InvOf5, 0.0_Kr /)
             Elems_db(iE)%BF(1,3) = (/ InvOf5, 0.0_Kr /)
             Elems_db(iE)%BF(1,4) = (/ InvOf5, 0.0_Kr /)
                                    
             Elems_db(iE)%BF(2,1) = (/ InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(2,2) = (/ InvOf5, 0.0_Kr /)
             Elems_db(iE)%BF(2,3) = (/ 3.0_Kr * InvOf5, 0.0_Kr /)
             Elems_db(iE)%BF(2,4) = (/ InvOf5, 0.0_Kr /)

             Elems_db(iE)%BF(3,1) = (/ InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(3,2) = (/ InvOf5, 0.0_Kr /)
             Elems_db(iE)%BF(3,3) = (/ InvOf5, 0.0_Kr /)
             Elems_db(iE)%BF(3,4) = (/ 3.0_Kr * InvOf5, 0.0_Kr /)


             Elems_db(iE)%BF(4,1) = (/ 0.0_Kr, InvOf3 /)
             Elems_db(iE)%BF(4,2) = (/ 0.0_Kr, 3.0_Kr * InvOf5 /)
             Elems_db(iE)%BF(4,3) = (/ 0.0_Kr, InvOf5 /)
             Elems_db(iE)%BF(4,4) = (/ 0.0_Kr, InvOf5 /)
                                    
             Elems_db(iE)%BF(5,1) = (/ 0.0_Kr, InvOf3 /)
             Elems_db(iE)%BF(5,2) = (/ 0.0_Kr, InvOf5 /)
             Elems_db(iE)%BF(5,3) = (/ 0.0_Kr, 3.0_Kr * InvOf5 /)
             Elems_db(iE)%BF(5,4) = (/ 0.0_Kr, InvOf5 /)

             Elems_db(iE)%BF(6,1) = (/ 0.0_Kr, InvOf3 /)
             Elems_db(iE)%BF(6,2) = (/ 0.0_Kr, InvOf5 /)
             Elems_db(iE)%BF(6,3) = (/ 0.0_Kr, InvOf5 /)
             Elems_db(iE)%BF(6,4) = (/ 0.0_Kr, 3.0_Kr * InvOf5 /)
          Case Default
             Write(*,*) 'Order not implemented for this element type ',       &
                  & Trim(Geom%Elem_Blk(iBlk)%Type), Order
             STOP
          End Select
       Case Default
          Write(*,*) 'Element type not implemented yet...',                   &
               & Geom%Elem_Blk(iBlk)%Type
          STOP
       End Select

       Elems_db(iE)%Der_BF(:,:)%XX = 0.0_Kr
       Elems_db(iE)%Der_BF(:,:)%YY = 0.0_Kr
       Elems_db(iE)%Der_BF(:,:)%XY = 0.0_Kr
       Elems_db(iE)%Der_BF(:,:)%YX = 0.0_Kr

       Elems_db(iE)%Der_BF(1,:)%XX =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y) * Jac_Inv
       Elems_db(iE)%Der_BF(1,:)%XY =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X) * Jac_Inv

       Elems_db(iE)%Der_BF(2,:)%XX =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y) * Jac_Inv
       Elems_db(iE)%Der_BF(2,:)%XY =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X) * Jac_Inv

       Elems_db(iE)%Der_BF(3,:)%XX =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y) * Jac_Inv
       Elems_db(iE)%Der_BF(3,:)%XY =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X) * Jac_Inv

       Elems_db(iE)%Der_BF(4,:)%YX =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y) * Jac_Inv
       Elems_db(iE)%Der_BF(4,:)%YY =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X) * Jac_Inv

       Elems_db(iE)%Der_BF(5,:)%YX =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y) * Jac_Inv
       Elems_db(iE)%Der_BF(5,:)%YY =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X) * Jac_Inv

       Elems_db(iE)%Der_BF(6,:)%YX =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y) * Jac_Inv
       Elems_db(iE)%Der_BF(6,:)%YY =                                          &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X) * Jac_Inv
    End Do Do_eNum
    DeAllocate(Elements)
  End Subroutine Init_Gauss_EXO_2D
    
  Subroutine Init_Gauss_EXO_2D_Elast(Elems_db, Nodes_db, Geom, Order, Elem,   &
       & Elem_List)
    Type (Element2D_Elast), Dimension(:), Pointer       :: Elems_db
    Type (Node2D), Dimension(:), Pointer                :: Nodes_db
    Type (EXO_Geom_Info), intent(IN)                    :: Geom
    Integer                                             :: Order
    Integer, Intent(IN), Optional                       :: Elem
    Integer, Dimension(:), Pointer, Optional            :: Elem_List
        
    Integer                                             :: iE, iBlk, eNum
    Integer                                             :: NB_BF, NB_Gauss
    Real(Kind = Kr)                                     :: Jac, Jac_Inv
    Integer, Dimension(:), Pointer                      :: Elements

    ! Make sure the databases are allocated
    If ((.NOT. Associated(Elems_db)) .OR. (.NOT. Associated(Nodes_db))) Then
       Print*, 'One of the required databases is not allocated...'
       Stop
    EndIf

    !!! Test args here, and allocate elem_list accordingly, if necessary
    If ( Present(Elem) .AND. Present(Elem_List) ) Then
       Print*, 'Error: Init_Gauss_EXO_2D_Elast, cannot specify Elem and',     &
            &  'Elem_List at the same time'
       STOP
    ElseIf ( Present(Elem) ) Then
       Allocate (Elements(1))
       Elements = Elem
    ElseIf ( .NOT. Present(Elem_List) ) Then
       Allocate(Elements(Geom%Num_Elems))
       Elements = (/ (iE, iE=1,Geom%Num_Elems) /)
    Else
       Allocate(Elements(Size(ELem_List)))
       Elements = Elem_List
    End If
    
    Do_eNum: Do eNum = 1, Size(Elements)
       iE = Elements(eNum)
       iBlk = Elems_db(iE)%ID_EL
       Select case (Trim(Geom%Elem_Blk(iBlk)%Type))
       Case('TRI3')
          Nb_BF = 6
          ! Base functions: (X,0), (Y,0), (1-X-Y,0)
          !                 (0,X), (0,Y), (0,1-X-Y)  
          Jac = Area_Tri_2D(Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord,           &
               & Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord,                      &
               & Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord) * 2.0_Kr
          Jac_Inv = 1.0_Kr / Jac
          Select Case (Order)
          Case(1)
             Nb_Gauss = 1
             Elems_db(iE)%Nb_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(NB_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF, NB_Gauss))
             Allocate(Elems_db(iE)%GradS_BF(Nb_BF, NB_Gauss))
             Elems_db(iE)%Gauss_C(1) = Jac / 2.0_Kr
             Elems_db(iE)%BF(1,1) = (/ InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(2,1) = (/ InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(3,1) = (/ InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(4,1) = (/ 0.0_Kr, InvOf3 /)
             Elems_db(iE)%BF(5,1) = (/ 0.0_Kr, InvOf3 /)
             Elems_db(iE)%BF(6,1) = (/ 0.0_Kr, InvOf3 /)
             
          Case(2)
             Nb_Gauss = 3
             Elems_db(iE)%Nb_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(NB_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF, NB_Gauss))
             Allocate(Elems_db(iE)%GradS_BF(Nb_BF, NB_Gauss))
             Elems_db(iE)%Gauss_C = Jac * InvOf6

             Elems_db(iE)%BF(1,1) = (/ InvOf6, 0.0_Kr /)
             Elems_db(iE)%BF(1,2) = (/ 2.0_Kr * InvOf3, 0.0_Kr/)
             Elems_db(iE)%BF(1,3) = (/ InvOf6, 0.0_Kr/)
                          
             Elems_db(iE)%BF(2,1) = (/ InvOf6, 0.0_Kr /)
             Elems_db(iE)%BF(2,2) = (/ InvOf6, 0.0_Kr /)
             Elems_db(iE)%BF(2,3) = (/ 2.0_Kr * InvOf3, 0.0_Kr /)
                          
             Elems_db(iE)%BF(3,1) = (/ 2.0_Kr * InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(3,2) = (/ InvOf6, 0.0_Kr /)
             Elems_db(iE)%BF(3,3) = (/ InvOf6, 0.0_Kr /)

             Elems_db(iE)%BF(4,1) = (/ 0.0_Kr, InvOf6 /)
             Elems_db(iE)%BF(4,2) = (/ 0.0_Kr, 2.0_Kr * InvOf3 /)
             Elems_db(iE)%BF(4,3) = (/ 0.0_Kr, InvOf6 /)
                          
             Elems_db(iE)%BF(5,1) = (/ 0.0_Kr, InvOf6 /)
             Elems_db(iE)%BF(5,2) = (/ 0.0_Kr, InvOf6 /)
             Elems_db(iE)%BF(5,3) = (/ 0.0_Kr, 2.0_Kr * InvOf3 /)

             Elems_db(iE)%BF(6,1) = (/ 0.0_Kr, 2.0_Kr * InvOf3 /)
             Elems_db(iE)%BF(6,2) = (/ 0.0_Kr, InvOf6 /)
             Elems_db(iE)%BF(6,3) = (/ 0.0_Kr, InvOf6 /)
          Case(3)
             Nb_Gauss = 4
             Elems_db(iE)%Nb_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(NB_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF, NB_Gauss))
             Allocate(Elems_db(iE)%GradS_BF(Nb_BF, NB_Gauss))
             Elems_db(iE)%Gauss_C = Jac * 25.0_Kr / 48.0_Kr / 2.0_Kr
             Elems_db(iE)%Gauss_C(1) = -Jac * 9.0_Kr / 16.0_Kr / 2.0_Kr

             Elems_db(iE)%BF(1,1) = (/ InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(1,2) = (/ 3.0_Kr * InvOf5, 0.0_Kr /)
             Elems_db(iE)%BF(1,3) = (/ InvOf5, 0.0_Kr /)
             Elems_db(iE)%BF(1,4) = (/ InvOf5, 0.0_Kr /)
                                    
             Elems_db(iE)%BF(2,1) = (/ InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(2,2) = (/ InvOf5, 0.0_Kr /)
             Elems_db(iE)%BF(2,3) = (/ 3.0_Kr * InvOf5, 0.0_Kr /)
             Elems_db(iE)%BF(2,4) = (/ InvOf5, 0.0_Kr /)

             Elems_db(iE)%BF(3,1) = (/ InvOf3, 0.0_Kr /)
             Elems_db(iE)%BF(3,2) = (/ InvOf5, 0.0_Kr /)
             Elems_db(iE)%BF(3,3) = (/ InvOf5, 0.0_Kr /)
             Elems_db(iE)%BF(3,4) = (/ 3.0_Kr * InvOf5, 0.0_Kr /)


             Elems_db(iE)%BF(4,1) = (/ 0.0_Kr, InvOf3 /)
             Elems_db(iE)%BF(4,2) = (/ 0.0_Kr, 3.0_Kr * InvOf5 /)
             Elems_db(iE)%BF(4,3) = (/ 0.0_Kr, InvOf5 /)
             Elems_db(iE)%BF(4,4) = (/ 0.0_Kr, InvOf5 /)
                                    
             Elems_db(iE)%BF(5,1) = (/ 0.0_Kr, InvOf3 /)
             Elems_db(iE)%BF(5,2) = (/ 0.0_Kr, InvOf5 /)
             Elems_db(iE)%BF(5,3) = (/ 0.0_Kr, 3.0_Kr * InvOf5 /)
             Elems_db(iE)%BF(5,4) = (/ 0.0_Kr, InvOf5 /)

             Elems_db(iE)%BF(6,1) = (/ 0.0_Kr, InvOf3 /)
             Elems_db(iE)%BF(6,2) = (/ 0.0_Kr, InvOf5 /)
             Elems_db(iE)%BF(6,3) = (/ 0.0_Kr, InvOf5 /)
             Elems_db(iE)%BF(6,4) = (/ 0.0_Kr, 3.0_Kr * InvOf5 /)
          Case Default
             Write(*,*) 'Order not implemented for this element type ',       &
                  & Trim(Geom%Elem_Blk(iBlk)%Type), Order
             STOP
          End Select
       Case Default
          Write(*,*) 'Element type not implemented yet...',                   &
               & Geom%Elem_Blk(iBlk)%Type
          STOP
       End Select

       Elems_db(iE)%GradS_BF(:,:)%XX = 0.0_Kr
       Elems_db(iE)%GradS_BF(:,:)%YY = 0.0_Kr
       Elems_db(iE)%GradS_BF(:,:)%XY = 0.0_Kr

       Elems_db(iE)%GradS_BF(1,:)%XX =                                       &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y) * Jac_Inv
       Elems_db(iE)%GradS_BF(1,:)%XY =                                       &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X) * Jac_Inv * .5_Kr

       Elems_db(iE)%GradS_BF(2,:)%XX =                                       &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y) * Jac_Inv
       Elems_db(iE)%GradS_BF(2,:)%XY =                                       &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X) * Jac_Inv * .5_Kr

       Elems_db(iE)%GradS_BF(3,:)%XX =                                       &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y) * Jac_Inv
       Elems_db(iE)%GradS_BF(3,:)%XY =                                       &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X) * Jac_Inv * .5_Kr

       Elems_db(iE)%GradS_BF(4,:)%XY =                                       &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y) * Jac_Inv * .5_Kr
       Elems_db(iE)%GradS_BF(4,:)%YY =                                       &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X) * Jac_Inv

       Elems_db(iE)%GradS_BF(5,:)%XY =                                       &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y) * Jac_Inv * .5_Kr
       Elems_db(iE)%GradS_BF(5,:)%YY =                                       &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X) * Jac_Inv

       Elems_db(iE)%GradS_BF(6,:)%XY =                                       &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y) * Jac_Inv * .5_Kr
       Elems_db(iE)%GradS_BF(6,:)%YY =                                       &
            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X                      &
            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X) * Jac_Inv
!!$       Elems_db(iE)%GradS_BF(1,:)%XX =                                        &
!!$            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y                      &
!!$            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y) * Jac_Inv
!!$       Elems_db(iE)%GradS_BF(1,:)%YY = 0.0_Kr
!!$       Elems_db(iE)%GradS_BF(1,:)%XY =                                        &
!!$            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X                      &
!!$            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X) * Jac_Inv * .5_Kr
!!$
!!$       Elems_db(iE)%GradS_BF(2,:)%XX =                                        &
!!$            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y                      &
!!$            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y) * Jac_Inv
!!$       Elems_db(iE)%GradS_BF(2,:)%YY = 0.0_Kr
!!$       Elems_db(iE)%GradS_BF(2,:)%XY =                                        &
!!$            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X                      &
!!$            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X) * Jac_Inv * .5_Kr
!!$
!!$       Elems_db(iE)%GradS_BF(3,:)%XX =                                        &
!!$            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y                      &
!!$            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y) * Jac_Inv
!!$       Elems_db(iE)%GradS_BF(3,:)%YY = 0.0_Kr
!!$       Elems_db(iE)%GradS_BF(3,:)%XY =                                        &
!!$            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X                      &
!!$            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X) * Jac_Inv * .5_Kr
!!$
!!$       Elems_db(iE)%GradS_BF(4,:)%XX = 0.0_Kr
!!$       Elems_db(iE)%GradS_BF(4,:)%YY =                                        &
!!$            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y                      &
!!$            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y) * Jac_Inv
!!$       Elems_db(iE)%GradS_BF(4,:)%XY =                                        &
!!$            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X                      &
!!$            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X) * Jac_Inv * .5_Kr
!!$
!!$       Elems_db(iE)%GradS_BF(5,:)%XX = 0.0_Kr
!!$       Elems_db(iE)%GradS_BF(5,:)%YY =                                        &
!!$            &  (Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%Y                      &
!!$            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y) * Jac_Inv
!!$       Elems_db(iE)%GradS_BF(5,:)%XY =                                        &
!!$            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X                      &
!!$            & - Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord%X) * Jac_Inv * .5_Kr
!!$
!!$       Elems_db(iE)%GradS_BF(6,:)%XX = 0.0_Kr
!!$       Elems_db(iE)%GradS_BF(6,:)%YY =                                        &
!!$            &  (Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%Y                      &
!!$            & - Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%Y) * Jac_Inv
!!$       Elems_db(iE)%GradS_BF(6,:)%XY =                                        &
!!$            &  (Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord%X                      &
!!$            & - Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord%X) * Jac_Inv * .5_Kr
    End Do Do_eNum
    DeAllocate(Elements)
  End Subroutine Init_Gauss_EXO_2D_Elast
  
  
  
  
  Subroutine InitGauss_2D_Elast_P1(Elem_bd, Nodes_bd, Order)
    Type (Element2D_Elast), Dimension(:), Pointer       :: Elem_bd
    Type (Node2D), Dimension(:), Pointer                :: Nodes_bd
    Integer                                             :: Order

    Integer(Kind = Ki)                                  :: NE, iE
    Integer(Kind = Ki)                                  :: iBF
    Integer(Kind = Ki)                                  :: iGauss, Nb_Gauss
    Integer(Kind = Ki)                                  :: NB_BF
    Type (MatS2D), Dimension(:),Pointer                 :: DerSBF
    Real(Kind = Kr)                                     :: Area, InvOfArea
    Real(Kind = Kr)                                     :: rA, rB, rG1, rG2
    Type (Vect2D)                                       :: Top, Bot


    ! On verifie d'abord que les bd sont allouees :
    If ((.NOT. Associated(Elem_bd)) .OR. (.NOT. Associated(Nodes_bd))) Then
       Print*, 'One of the required databases is not allocated...'
       Stop
    EndIf

    NE = Size(Elem_bd)
    Top%X = 1.0_Kr
    Top%Y = 0.0_Kr
    Bot%X = 0.0_Kr
    Bot%Y = 1.0_Kr

    OrderG : Select Case(Order)
    Case (2)
       Nb_Gauss = 3
       Nb_BF = 6
       Allocate(DerSBF(Nb_BF))
       Do_iE: Do iE = 1, NE
          Elem_bd(iE)%NB_Gauss = Nb_Gauss
          Allocate(Elem_bd(iE)%BF(Nb_BF,Nb_Gauss))
          Allocate(Elem_bd(iE)%GradS_BF(NB_BF,Nb_Gauss))

          Elem_bd(iE)%BF(1,1) = 0.0_Kr * Top
          Elem_bd(iE)%BF(1,2) = 0.5_Kr * Top
          Elem_bd(iE)%BF(1,3) = 0.5_Kr * Top

          Elem_bd(iE)%BF(2,1) = 0.0_Kr * Bot
          Elem_bd(iE)%BF(2,2) = 0.5_Kr * Bot
          Elem_bd(iE)%BF(2,3) = 0.5_Kr * Bot

          Elem_bd(iE)%BF(3,1) = 0.5_Kr * Top
          Elem_bd(iE)%BF(3,2) = 0.0_Kr * Top
          Elem_bd(iE)%BF(3,3) = 0.5_Kr * Top

          Elem_bd(iE)%BF(4,1) = 0.5_Kr * Bot
          Elem_bd(iE)%BF(4,2) = 0.0_Kr * Bot
          Elem_bd(iE)%BF(4,3) = 0.5_Kr * Bot

          Elem_bd(iE)%BF(5,1) = 0.5_Kr * Top
          Elem_bd(iE)%BF(5,2) = 0.5_Kr * Top
          Elem_bd(iE)%BF(5,3) = 0.0_Kr * Top

          Elem_bd(iE)%BF(6,1) = 0.5_Kr * Bot
          Elem_bd(iE)%BF(6,2) = 0.5_Kr * Bot
          Elem_bd(iE)%BF(6,3) = 0.0_Kr * Bot

          Area = Area_Tri_2D(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord, &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord,&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord)

          InvOfArea = 1.0_Kr / Area

          DerSBF(1)%XX = -(Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%Y - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y)            &
               & * InvOf2 * InvOfArea
          DerSBF(1)%YY = 0.0_Kr
          DerSBF(1)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%X)            &
               & * InvOf4 * InvOfArea

          DerSBF(2)%XX = 0.0_Kr
          DerSBF(2)%YY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%X)            &
               & * InvOf2 * InvOfArea
          DerSBF(2)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%Y - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y)            &
               & * InvOf4 * InvOfArea

          DerSBF(3)%XX = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%Y)            &
               & * InvOf2 * InvOfArea
          DerSBF(3)%YY = 0.0_Kr
          DerSBF(3)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%X - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X)            &
               & * InvOf4 * InvOfArea

          DerSBF(4)%XX = 0.0_Kr
          DerSBF(4)%YY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%X - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X)            &
               & * InvOf2 * InvOfArea
          DerSBF(4)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%Y)            &
               & * InvOf4 * InvOfArea

          DerSBF(5)%XX = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y)            &
               & * InvOf2 * InvOfArea
          DerSBF(5)%YY = 0.0_Kr
          DerSBF(5)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X)            &
               & * InvOf4 * InvOfArea

          DerSBF(6)%XX = 0.0_Kr
          DerSBF(6)%YY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X)            &
               & * InvOf2 * InvOfArea
          DerSBF(6)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y)            &
               & * InvOf4 * InvOfArea

          Do iGauss = 1, NB_Gauss
             Elem_bd(iE)%GradS_BF(:, iGauss) = DerSBF
          EndDo

          Allocate(Elem_bd(iE)%Gauss_C(Nb_Gauss))
          Elem_bd(iE)%Gauss_C = InvOf3 * Area


       EndDo Do_iE
       DeAllocate(DerSBF)

    Case (3)
       Nb_Gauss = 4
       Nb_BF = 6
       rA = -27.0_Kr / 48.0_Kr
       rB = 25.0_Kr / 48.0_Kr
       Allocate(DerSBF(Nb_BF))
       Do iE = 1, NE
          Elem_bd(iE)%NB_Gauss = Nb_Gauss
          Allocate(Elem_bd(iE)%BF(Nb_BF,Nb_Gauss))
          Allocate(Elem_bd(iE)%GradS_BF(NB_BF,Nb_Gauss))
          Elem_bd(iE)%BF(1,1) = InvOf3 * Top
          Elem_bd(iE)%BF(1,2) = 0.6_Kr * Top
          Elem_bd(iE)%BF(1,3) = 0.2_Kr * Top
          Elem_bd(iE)%BF(1,4) = 0.2_Kr * Top

          Elem_bd(iE)%BF(2,1) = InvOf3 * Bot
          Elem_bd(iE)%BF(2,2) = 0.6_Kr * Bot
          Elem_bd(iE)%BF(2,3) = 0.2_Kr * Bot
          Elem_bd(iE)%BF(2,4) = 0.2_Kr * Bot

          Elem_bd(iE)%BF(3,1) = InvOf3 * Top
          Elem_bd(iE)%BF(3,2) = 0.2_Kr * Top
          Elem_bd(iE)%BF(3,3) = 0.6_Kr * Top
          Elem_bd(iE)%BF(3,4) = 0.2_Kr * Top

          Elem_bd(iE)%BF(4,1) = InvOf3 * Bot
          Elem_bd(iE)%BF(4,2) = 0.2_Kr * Bot
          Elem_bd(iE)%BF(4,3) = 0.6_Kr * Bot
          Elem_bd(iE)%BF(4,4) = 0.2_Kr * Bot

          Elem_bd(iE)%BF(5,1) = InvOf3 * Top
          Elem_bd(iE)%BF(5,2) = 0.2_Kr * Top
          Elem_bd(iE)%BF(5,3) = 0.2_Kr * Top
          Elem_bd(iE)%BF(5,4) = 0.6_Kr * Top

          Elem_bd(iE)%BF(6,1) = InvOf3 * Bot
          Elem_bd(iE)%BF(6,2) = 0.2_Kr * Bot
          Elem_bd(iE)%BF(6,3) = 0.2_Kr * Bot
          Elem_bd(iE)%BF(6,4) = 0.6_Kr * Bot

          Area = Area_Tri_2D(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord, &
               &    Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord,&
               &    Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord)

          InvOfArea = 1.0_Kr / Area


          DerSBF(1)%XX = -(Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%Y - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y) * &
               & 0.5_Kr * InvOfArea
          DerSBF(1)%YY = 0.0_Kr
          DerSBF(1)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%X) *&
               & 0.5_Kr * InvOfArea * InvOf2

          DerSBF(2)%XX = 0.0_Kr
          DerSBF(2)%YY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%X) *&
               & 0.5_Kr * InvOfArea
          DerSBF(2)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%Y -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y) *&
               & 0.5_Kr * InvOfArea * InvOf2

          DerSBF(3)%XX = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%Y) *&
               & 0.5_Kr * InvOfArea
          DerSBF(3)%YY = 0.0_Kr
          DerSBF(3)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%X -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X) *&
               & 0.5_Kr * InvOfArea * InvOf2

          DerSBF(4)%XX = 0.0_Kr
          DerSBF(4)%YY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%X -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X) *&
               & 0.5_Kr * InvOfArea
          DerSBF(4)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%Y) *&
               & 0.5_Kr * InvOfArea * InvOf2

          DerSBF(5)%XX = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y) *&
               & 0.5_Kr * InvOfArea
          DerSBF(5)%YY = 0.0_Kr
          DerSBF(5)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X) * &
               & 0.5_Kr * InvOfArea * InvOf2

          DerSBF(6)%XX = 0.0_Kr
          DerSBF(6)%YY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X) *&
               & 0.5_Kr * InvOfArea
          DerSBF(6)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y) *&
               & 0.5_Kr * InvOfArea * InvOf2

          Do iGauss = 1, Nb_Gauss
             Elem_bd(iE)%GradS_BF(:,iGauss) = DerSBF
          EndDo
          Allocate(Elem_bd(iE)%Gauss_C(Nb_Gauss))
          Elem_bd(iE)%Gauss_C(1) = Area * rA
          Elem_bd(iE)%Gauss_C(2) = Area * rB
          Elem_bd(iE)%Gauss_C(3) = Area * rB
          Elem_bd(iE)%Gauss_C(4) = Area * rB
       EndDo
       DeAllocate(DerSBF)
       ! OK 
    Case (4)
       Nb_Gauss = 6
       Nb_BF = 6
       rA = 0.445948490915965_Kr
       rB = 0.091576213509771_Kr
       rG1 = 0.111690794839005_Kr * 2.0_Kr
       rG2 = 0.054975871827661_Kr * 2.0_Kr
       Allocate(DerSBF(Nb_BF))
       Do iE = 1, NE
          Elem_bd(iE)%NB_Gauss = Nb_Gauss
          Allocate(Elem_bd(iE)%BF(Nb_BF,Nb_Gauss))
          Allocate(Elem_bd(iE)%GradS_BF(NB_BF,Nb_Gauss))
          Elem_bd(iE)%BF(1,1) = (1.0_Kr - 2.0_Kr * rA)  * Top
          Elem_bd(iE)%BF(1,2) = rA                      * Top
          Elem_bd(iE)%BF(1,3) = rA                      * Top
          Elem_bd(iE)%BF(1,4) = (1.0_Kr - 2.0_Kr * rB)  * Top
          Elem_bd(iE)%BF(1,5) = rB                      * Top
          Elem_bd(iE)%BF(1,6) = rB                      * Top

          Elem_bd(iE)%BF(2,1) = (1.0_Kr - 2.0_Kr * rA)  * Bot
          Elem_bd(iE)%BF(2,2) = rA                      * Bot
          Elem_bd(iE)%BF(2,3) = rA                      * Bot
          Elem_bd(iE)%BF(2,4) = (1.0_Kr - 2.0_Kr * rB)  * Bot
          Elem_bd(iE)%BF(2,5) = rB                      * Bot
          Elem_bd(iE)%BF(2,6) = rB                      * Bot

          Elem_bd(iE)%BF(3,1) = rA                      * Top
          Elem_bd(iE)%BF(3,2) = (1.0_Kr - 2.0_Kr * rA)  * Top
          Elem_bd(iE)%BF(3,3) = rA                      * Top
          Elem_bd(iE)%BF(3,4) = rB                      * Top
          Elem_bd(iE)%BF(3,5) = (1.0_Kr - 2.0_Kr * rB)  * Top
          Elem_bd(iE)%BF(3,6) = rB                      * Top

          Elem_bd(iE)%BF(4,1) = rA                      * Bot
          Elem_bd(iE)%BF(4,2) = (1.0_Kr - 2.0_Kr * rA)  * Bot
          Elem_bd(iE)%BF(4,3) = rA                      * Bot
          Elem_bd(iE)%BF(4,4) = rB                      * Bot
          Elem_bd(iE)%BF(4,5) = (1.0_Kr - 2.0_Kr * rB)  * Bot
          Elem_bd(iE)%BF(4,6) = rB                      * Bot

          Elem_bd(iE)%BF(5,1) = rA                      * Top
          Elem_bd(iE)%BF(5,2) = rA                      * Top
          Elem_bd(iE)%BF(5,3) = (1.0_Kr - 2.0_Kr * rA)  * Top
          Elem_bd(iE)%BF(5,4) = rB                      * Top
          Elem_bd(iE)%BF(5,5) = rB                      * Top
          Elem_bd(iE)%BF(5,6) = (1.0_Kr - 2.0_Kr * rB)  * Top

          Elem_bd(iE)%BF(6,1) = rA                      * Bot
          Elem_bd(iE)%BF(6,2) = rA                      * Bot
          Elem_bd(iE)%BF(6,3) = (1.0_Kr - 2.0_Kr * rA)  * Bot
          Elem_bd(iE)%BF(6,4) = rB                      * Bot
          Elem_bd(iE)%BF(6,5) = rB                      * Bot
          Elem_bd(iE)%BF(6,6) = (1.0_Kr - 2.0_Kr * rB)  * Bot

          Area = Area_Tri_2D(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord, &
               &    Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord,&
               &    Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord)

          InvOfArea = 1.0_Kr / Area


          DerSBF(1)%XX = -(Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%Y - &
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y) * &
               & 0.5_Kr * InvOfArea
          DerSBF(1)%YY = 0.0_Kr
          DerSBF(1)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%X) *&
               & 0.5_Kr * InvOfArea * InvOf2

          DerSBF(2)%XX = 0.0_Kr
          DerSBF(2)%YY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%X) *&
               & 0.5_Kr * InvOfArea
          DerSBF(2)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%Y -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y) *&
               & 0.5_Kr * InvOfArea * InvOf2

          DerSBF(3)%XX = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%Y) *&
               & 0.5_Kr * InvOfArea
          DerSBF(3)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%X -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X) *&
               & 0.5_Kr * InvOfArea * InvOf2
          DerSBF(3)%YY = 0.0_Kr

          DerSBF(4)%XX = 0.0_Kr
          DerSBF(4)%YY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%X -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X) *&
               & 0.5_Kr * InvOfArea
          DerSBF(4)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(5))%Coord%Y) *&
               & 0.5_Kr * InvOfArea * InvOf2

          DerSBF(5)%XX = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y) *&
               & 0.5_Kr * InvOfArea
          DerSBF(5)%YY = 0.0_Kr
          DerSBF(5)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X) * &
               & 0.5_Kr * InvOfArea * InvOf2

          DerSBF(6)%XX = 0.0_Kr
          DerSBF(6)%YY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X) *&
               & 0.5_Kr * InvOfArea
          DerSBF(6)%XY = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y -&
               & Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y) *&
               & 0.5_Kr * InvOfArea * InvOf2

          Do iGauss = 1, Nb_Gauss
             Elem_bd(iE)%GradS_BF(:,iGauss) = DerSBF
          EndDo
          Allocate(Elem_bd(iE)%Gauss_C(Nb_Gauss))
          Elem_bd(iE)%Gauss_C(1) = rG1 * Area
          Elem_bd(iE)%Gauss_C(2) = rG1 * Area
          Elem_bd(iE)%Gauss_C(3) = rG1 * Area
          Elem_bd(iE)%Gauss_C(4) = rG2 * Area
          Elem_bd(iE)%Gauss_C(5) = rG2 * Area
          Elem_bd(iE)%Gauss_C(6) = rG2 * Area
       End Do
       DeAllocate(DerSBF)      

    Case Default
       Print*, 'This order is not implemented yet...', Order
       Stop

    End Select OrderG


  End Subroutine InitGauss_2D_Elast_P1

  Subroutine InitGauss_2D_Scal_P1(Elem_bd, Nodes_bd, Order)
    Type (Element2D_Scal), Dimension(:), Pointer       :: Elem_bd
    Type (Node2D), Dimension(:), Pointer                :: Nodes_bd
    Integer                                             :: Order

    Integer(Kind = Ki)                                  :: NE, iE
!    Integer(Kind = Ki)                                 :: iBF
    Integer(Kind = Ki)                                  :: iGauss, Nb_Gauss
    Integer(Kind = Ki)                                  :: NB_BF
    Type (Vect2D), Dimension(:),Pointer                 :: GradBF
    Real(Kind = Kr)                                     :: Area, InvOfArea
    Real(Kind = Kr)                                     :: rG1, rG2, rA, rB

    ! On verifie d'abord que les bd sont allouees :
    If ((.NOT. Associated(Elem_bd)) .OR. (.NOT. Associated(Nodes_bd))) Then
       Print*, 'One of the required databases is not allocated...'
       Stop
    EndIf

    NE = Size(Elem_bd)

    OrderG : Select Case(Order)
    Case (2)
       Nb_Gauss = 3
       Nb_BF = 3
       Allocate(GradBF(Nb_BF))
       Do iE = 1, NE
          Elem_bd(iE)%NB_Gauss = Nb_Gauss
          Allocate(Elem_bd(iE)%BF(Nb_BF,Nb_Gauss))
          Allocate(Elem_bd(iE)%Grad_BF(NB_BF,Nb_Gauss))
          Elem_bd(iE)%BF(1,1) = 0.0_Kr
          Elem_bd(iE)%BF(1,2) = 0.5_Kr
          Elem_bd(iE)%BF(1,3) = 0.5_Kr

          Elem_bd(iE)%BF(2,1) = 0.5_Kr
          Elem_bd(iE)%BF(2,2) = 0.0_Kr
          Elem_bd(iE)%BF(2,3) = 0.5_Kr

          Elem_bd(iE)%BF(3,1) = 0.5_Kr
          Elem_bd(iE)%BF(3,2) = 0.5_Kr
          Elem_bd(iE)%BF(3,3) = 0.0_Kr

          Area = Area_Tri_2D(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord, &
               &    Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord,&
               &    Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord)

          InvOfArea = 1.0_Kr / Area

          GradBF(1)%X = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord%Y) * &
               &                   0.5_Kr * InvOfArea

          GradBF(1)%Y = -(Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord%X - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X) * &
               &                   0.5_Kr * InvOfArea

          GradBF(2)%X = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y) * &
               &                   0.5_Kr * InvOfArea

          GradBF(2)%Y = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X) * &
               &                   0.5_Kr * InvOfArea

          GradBF(3)%X = -(Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord%Y - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y) * &
               &                   0.5_Kr * InvOfArea

          GradBF(3)%Y = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord%X) * &
               &                   0.5_Kr * InvOfArea


          Do iGauss = 1, Nb_Gauss
!             Do iBF = 1, Nb_BF
!                Elem_bd(iE)%Grad_BF(iBF, iGauss)%X = GradBF(iBF)%X
!                Elem_bd(iE)%Grad_BF(iBF, iGauss)%Y = GradBF(iBF)%Y
!             EndDo
             ELEM_BD(IE)%GRAD_BF(:,iGauss) = GradBF
          EndDo
          Allocate(Elem_bd(iE)%Gauss_C(Nb_Gauss))
!          Do iGauss = 1, Nb_Gauss
!             Elem_bd(iE)%Gauss_C(iGauss) = Area * InvOf3
!          EndDo
          Elem_bd(iE)%Gauss_C = Area * InvOf3

       EndDo


    Case (3)
       Nb_Gauss = 4
       Nb_BF = 3
       Allocate(GradBF(Nb_BF))
       Do iE = 1, NE
          Elem_bd(iE)%NB_Gauss = Nb_Gauss
          Allocate(Elem_bd(iE)%BF(Nb_BF,Nb_Gauss))
          Allocate(Elem_bd(iE)%Grad_BF(NB_BF,Nb_Gauss))
          Elem_bd(iE)%BF(1,1) = InvOf3
          Elem_bd(iE)%BF(1,2) = 0.6_Kr
          Elem_bd(iE)%BF(1,3) = 0.2_Kr
          Elem_bd(iE)%BF(1,4) = 0.2_Kr

          Elem_bd(iE)%BF(2,1) = InvOf3
          Elem_bd(iE)%BF(2,2) = 0.2_Kr
          Elem_bd(iE)%BF(2,3) = 0.6_Kr
          Elem_bd(iE)%BF(2,4) = 0.2_Kr

          Elem_bd(iE)%BF(3,1) = InvOf3
          Elem_bd(iE)%BF(3,2) = 0.2_Kr
          Elem_bd(iE)%BF(3,3) = 0.2_Kr
          Elem_bd(iE)%BF(3,4) = 0.6_Kr

          Area = Area_Tri_2D(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord, &
               &    Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord,&
               &    Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord)
          
          InvOfArea = 1.0_Kr / Area

          GradBF(1)%X = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord%Y) * &
               &                   0.5_Kr * InvOfArea

          GradBF(1)%Y = -(Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord%X - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X) * &
               &                   0.5_Kr * InvOfArea

          GradBF(2)%X = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y) * &
               &                   0.5_Kr * InvOfArea

          GradBF(2)%Y = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X) * &
               &                   0.5_Kr * InvOfArea

          GradBF(3)%X = -(Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord%Y - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y) * &
               &                   0.5_Kr * InvOfArea

          GradBF(3)%Y = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord%X) * &
               &                   0.5_Kr * InvOfArea


          Do iGauss = 1, Nb_Gauss
!             Do iBF = 1, Nb_BF
!                Elem_bd(iE)%Grad_BF(iBF, iGauss)%X = GradBF(iBF)%X
!                Elem_bd(iE)%Grad_BF(iBF, iGauss)%Y = GradBF(iBF)%Y
!             EndDo
             Elem_Bd(iE)%Grad_BF(:,iGauss) = GradBF
          EndDo
          Allocate(Elem_bd(iE)%Gauss_C(Nb_Gauss))
          rA = -27.0_Kr / 48.0_Kr
          rB = 25.0_Kr / 48.0_Kr
          Elem_bd(iE)%Gauss_C(1) =  Area * rA
          Elem_bd(iE)%Gauss_C(2) =  Area * rB
          Elem_bd(iE)%Gauss_C(3) =  Area * rB
          Elem_bd(iE)%Gauss_C(4) =  Area * rB


       EndDo
       DeAllocate(GradBF)

    Case (4)
       Nb_Gauss = 6
       Nb_BF = 3
       rA = 0.445948490915965_Kr
       rB = 0.091576213509771_Kr
       rG1 = 0.111690794839005_Kr * 2.0_Kr
       rG2 = 0.054975871827661_Kr * 2.0_Kr
       Allocate(GradBF(Nb_BF))
       Do iE = 1, NE
          Elem_bd(iE)%NB_Gauss = Nb_Gauss
          Allocate(Elem_bd(iE)%BF(Nb_BF,Nb_Gauss))
          Allocate(Elem_bd(iE)%Grad_BF(NB_BF,Nb_Gauss))
          Elem_bd(iE)%BF(1,1) = 1.0_Kr-2.0_Kr*rA
          Elem_bd(iE)%BF(1,2) = rA
          Elem_bd(iE)%BF(1,3) = rA
          Elem_bd(iE)%BF(1,4) = 1.0_Kr-2.0_Kr*rB
          Elem_bd(iE)%BF(1,5) = rB
          Elem_bd(iE)%BF(1,6) = rB

          Elem_bd(iE)%BF(2,1) = rA
          Elem_bd(iE)%BF(2,2) = 1.0_Kr-2.0_Kr*rA
          Elem_bd(iE)%BF(2,3) = rA
          Elem_bd(iE)%BF(2,4) = rB
          Elem_bd(iE)%BF(2,5) = 1.0_Kr-2.0_Kr*rB
          Elem_bd(iE)%BF(2,6) = rB

          Elem_bd(iE)%BF(3,1) = rA
          Elem_bd(iE)%BF(3,2) = rA
          Elem_bd(iE)%BF(3,3) = 1.0_Kr-2.0_Kr*rA
          Elem_bd(iE)%BF(3,4) = rB
          Elem_bd(iE)%BF(3,5) = rB
          Elem_bd(iE)%BF(3,6) = 1.0_Kr-2.0_Kr*rB

          Area = Area_Tri_2D(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord, &
               &    Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord,&
               &    Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord)
          
          InvOfArea = 1.0_Kr / Area

          GradBF(1)%X = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord%Y) * &
               &                   InvOfArea * 0.5_Kr

          GradBF(1)%Y = -(Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord%X - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X) * &
               &                   InvOfArea * 0.5_Kr

          GradBF(2)%X = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%Y) * &
               &                   InvOfArea * 0.5_Kr

          GradBF(2)%Y = -(Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord%X - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X) * &
               &                   InvOfArea * 0.5_Kr

          GradBF(3)%X = -(Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord%Y - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%Y) * &
               &                   InvOfArea * 0.5_Kr

          GradBF(3)%Y = -(Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord%X - &
               &                   Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord%X) * &
               &                   InvOfArea * 0.5_Kr


          Do iGauss = 1, Nb_Gauss
!             Do iBF = 1, Nb_BF
!                Elem_bd(iE)%Grad_BF(iBF, iGauss)%X = GradBF(iBF)%X
!                Elem_bd(iE)%Grad_BF(iBF, iGauss)%Y = GradBF(iBF)%Y
!             EndDo
             Elem_bd(iE)%Grad_BF(:,iGauss) = GradBF
          EndDo
          Allocate(Elem_bd(iE)%Gauss_C(Nb_Gauss))
          Elem_bd(iE)%Gauss_C(1) = rG1*Area
          Elem_bd(iE)%Gauss_C(2) = rG1*Area
          Elem_bd(iE)%Gauss_C(3) = rG1*Area
          Elem_bd(iE)%Gauss_C(4) = rG2*Area
          Elem_bd(iE)%Gauss_C(5) = rG2*Area
          Elem_bd(iE)%Gauss_C(6) = rG2*Area


       EndDo
       DeAllocate(GradBF)      

    Case Default
       Print*, 'This order is not implemented yet...', Order
       Stop

    End Select OrderG


  End Subroutine InitGauss_2D_Scal_P1

  Subroutine Init_Gauss_EXO_3D_Scal(Elems_db, Nodes_db, Geom, Order, Elem,    &
       & Elem_List)
    Type (Element3D_Scal), Dimension(:), Pointer        :: Elems_db
    Type (Node3D), Dimension(:), Pointer                :: Nodes_db
    Type (EXO_Geom_Info), intent(IN)                    :: Geom
    Integer                                             :: Order
    Integer, Intent(IN), Optional                       :: Elem
    Integer, Dimension(:), Pointer, Optional            :: Elem_List
    
    Integer                                             :: iE, iBlk, eNum
    Integer                                             :: NB_BF, NB_Gauss
    Real(Kind = Kr)                                     :: Jac, Jac_Inv

    Real(Kind = Kr)                                     :: r1, s1, r2, s2
    Real(Kind = Kr)                                     :: u1, v1

    Real(Kind = Kr), Dimension(:,:), Pointer            :: F, F_Inv
    ! F = linear transformation from the current element to the 3rd simplex
    Type(Vect3D)                                        :: Ve1, Ve2, Ve3, Ve4
    Logical                                             :: GJStat
    Type(Vect3D), Dimension(:), Pointer                 :: GaussP
    Integer, Dimension(:), Pointer                      :: Elements

    ! Make sure the databases are allocated
    If ((.NOT. Associated(Elems_db)) .OR. (.NOT. Associated(Nodes_db))) Then
       Print*, 'One of the required databases is not allocated...'
       Stop
    EndIf

    !!! Test args here, and allocate elem_list accordingly, if necessary
    If ( Present(Elem) .AND. Present(Elem_List) ) Then
       Print*, 'Error: Init_Gauss_EXO_3D_Scal, cannot specify Elem and',     &
            &  'Elem_List at the same time'
       STOP
    ElseIf ( Present(Elem) ) Then
       Allocate (Elements(1))
       Elements = Elem
    ElseIf ( .NOT. Present(Elem_List) ) Then
       Allocate(Elements(Geom%Num_Elems))
       Elements = (/ (iE, iE=1,Geom%Num_Elems) /)
    Else
       Allocate(Elements(Size(Elem_List)))
       Elements = Elem_List
    End If

    Do_eNum: Do eNum = 1, Size(Elements)
       iE = Elements(eNum)
       iBlk = Elems_db(iE)%ID_EL
       Elem_Type: Select case (Trim(Geom%Elem_Blk(iBlk)%Type))
       Case('TETRA4', 'TETRA')
          Allocate(F_Inv(3,3))
          Allocate(F(3, 3))
          Ve1 = Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord
          Ve2 = Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord
          Ve3 = Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord
          Ve4 = Nodes_db(Elems_db(iE)%ID_DoF(4))%Coord
!          Jac = |T| / |\hat T|
          Jac = Vol_Tetra_3D( Ve1, Ve2, Ve3, Ve4 ) * 6.0_Kr
          Jac_Inv = 1.0_Kr / Jac
          Nb_BF = 4
          ! Base Functions: 1-X-Y-Z, X, Y, Z

          OrderG : Select Case(Order)
          Case (1,2,3)
             ! 3rd ordeer cubature on a tetrahedron, cf
             ! J.E. Akin lecture notes / book

             Nb_Gauss = 5
             Allocate (GaussP(Nb_Gauss))

             GaussP%X = (/ 0.25_Kr, 0.5_Kr, InvOf6, InvOf6, InvOf6 /)  
             GaussP%Y = (/ 0.25_Kr, InvOf6, 0.5_Kr, InvOf6, InvOf6 /)
             GaussP%Z = (/ 0.25_Kr, InvOf6, InvOf6, 0.5_Kr, InvOf6 /)
             ! The coordinates of the Gauss Points in the reference element

             Elems_db(iE)%NB_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(Nb_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF,Nb_Gauss))
             Allocate(Elems_db(iE)%Grad_BF(NB_BF,Nb_Gauss))
             
             Elems_db(iE)%Gauss_C(1)     = -2.0_Kr / 15.0_Kr * Jac
             Elems_db(iE)%Gauss_C(2:5)   =  3.0_Kr / 40.0_Kr * Jac

             Elems_db(iE)%BF(1,:) = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(2,:) = GaussP%X
             Elems_db(iE)%BF(3,:) = GaussP%Y
             Elems_db(iE)%BF(4,:) = GaussP%Z

             F_Inv = reshape((/Ve2%X - Ve1%X, Ve2%Y - Ve1%Y, Ve2%Z - Ve1%Z,   &
               &               Ve3%X - Ve1%X, Ve3%Y - Ve1%Y, Ve3%Z - Ve1%Z,   &
               &               Ve4%X - Ve1%X, Ve4%Y - Ve1%Y, Ve4%Z - Ve1%Z/), &
               &             (/3,3/) )
             F = F_Inv
             Call GaussJordan_Inverse(F, GJStat)
             
             If ( .NOT. GJStat) Then
                Print*, 'Current element is degenerate, cannot compute'
                Print*, 'the linear transformation', iE
                STOP
             End If
             DeAllocate (GaussP)

          Case (4,5)
             ! 5th ordeer cubature on a tetrahedron, cf
             ! A.H. Stroud: "A fifth degree integration formula for the n-th 
             ! simplex SIAM J. Numer. Anal. Vol 6, No1, March 1969
             Nb_Gauss = 15
             r1 = (7.0_Kr - sqrt(15.0_Kr)) / 34.0_Kr
             s1 = 1.0_Kr - 3.0_Kr * r1
             r2 = (7.0_Kr + sqrt(15.0_Kr)) / 34.0_Kr
             s2 = 1.0_Kr - 3.0_Kr * r2
             
             u1 = (5.0_Kr + sqrt(15.0_Kr))/20.0_Kr
             v1 = (5.0_Kr - sqrt(15.0_Kr))/20.0_Kr 
             Allocate (GaussP(Nb_Gauss))
             GaussP%X = (/ 0.25_Kr, r1, s1, r1, r1, r2, s2, r2, r2,           &
                  &        u1, u1, v1, u1, v1, v1 /)
             GaussP%Y = (/ 0.25_Kr, r1, r1, s1, r1, r2, r2, s2, r2,           &
                  &        u1, v1, u1, v1, u1, v1 /)
             GaussP%Z = (/ 0.25_Kr, r1, r1, r1, s1, r2, r2, r2, s2,           &
                  &        v1, u1, u1, v1, v1, u1 /)
             ! The coordinates of the Gauss Points in the reference element

             Elems_db(iE)%NB_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(Nb_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF,Nb_Gauss))
             Allocate(Elems_db(iE)%Grad_BF(NB_BF,Nb_Gauss))
             
             Elems_db(iE)%Gauss_C(1)     = 16.0_Kr / 135.0_Kr * Jac * InvOf6
             Elems_db(iE)%Gauss_C(2:5)   = (2665.0_Kr +14.0_Kr*sqrt(15.0_Kr)) &
                  &                       / 37800.0_Kr * Jac * InvOf6
             Elems_db(iE)%Gauss_C(6:9)   = (2665.0_Kr -14.0_Kr*sqrt(15.0_Kr)) &
                  &                       / 37800.0_Kr * Jac * InvOf6
             Elems_db(iE)%Gauss_C(10:15) = 10.0_Kr / 189.0_Kr * Jac * InvOf6 


             Elems_db(iE)%BF(1,:) = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(2,:) = GaussP%X
             Elems_db(iE)%BF(3,:) = GaussP%Y
             Elems_db(iE)%BF(4,:) = GaussP%Z

             F_Inv = reshape((/Ve2%X - Ve1%X, Ve2%Y - Ve1%Y, Ve2%Z - Ve1%Z,   &
               &               Ve3%X - Ve1%X, Ve3%Y - Ve1%Y, Ve3%Z - Ve1%Z,   &
               &               Ve4%X - Ve1%X, Ve4%Y - Ve1%Y, Ve4%Z - Ve1%Z/), &
               &             (/3,3/) )
             F = F_Inv
             Call GaussJordan_Inverse(F, GJStat)
             If ( .NOT. GJStat) Then
                Print*, 'Current element is degenerate, cannot compute'
                Print*, 'the linear transformation', iE
                STOP
             End If
             DeAllocate (GaussP)
          Case Default
             Write(*,*) 'Order non implemented for TETRA4', Order
             STOP
          End Select OrderG

          Elems_db(iE)%Grad_BF(1,:)%X = (-F(1,1)-F(2,1)-F(3,1))
          Elems_db(iE)%Grad_BF(1,:)%Y = (-F(1,2)-F(2,2)-F(3,2))
          Elems_db(iE)%Grad_BF(1,:)%Z = (-F(1,3)-F(2,3)-F(3,3))
          
          Elems_db(iE)%Grad_BF(2,:)%X = F(1,1)
          Elems_db(iE)%Grad_BF(2,:)%Y = F(1,2)
          Elems_db(iE)%Grad_BF(2,:)%Z = F(1,3)
          
          Elems_db(iE)%Grad_BF(3,:)%X = F(2,1)
          Elems_db(iE)%Grad_BF(3,:)%Y = F(2,2)
          Elems_db(iE)%Grad_BF(3,:)%Z = F(2,3)
          
          Elems_db(iE)%Grad_BF(4,:)%X = F(3,1)
          Elems_db(iE)%Grad_BF(4,:)%Y = F(3,2)
          Elems_db(iE)%Grad_BF(4,:)%Z = F(3,3)
          
          DeAllocate (F)
          DeAllocate (F_Inv)          
       Case Default
          Write(*,*) 'Element type non implemented yet',                      &
               & Trim(Geom%Elem_Blk(iBlk)%Type)
       End Select Elem_Type
    End Do Do_eNum

    DeAllocate(Elements)
  End Subroutine Init_Gauss_EXO_3D_Scal

  Subroutine Init_Gauss_EXO_3D(Elems_db, Nodes_db, Geom, Order, Elem,         &
       & Elem_List)
    Type (Element3D), Dimension(:), Pointer             :: Elems_db
    Type (Node3D), Dimension(:), Pointer                :: Nodes_db
    Type (EXO_Geom_Info), intent(IN)                    :: Geom
    Integer                                             :: Order
    Integer, Intent(IN), Optional                       :: Elem
    Integer, Dimension(:), Pointer, Optional            :: Elem_List
    
    Integer                                             :: iE, iBlk, eNum
    Integer                                             :: NB_BF, NB_Gauss
    Real(Kind = Kr)                                     :: Jac, Jac_Inv

    Real(Kind = Kr)                                     :: r1, s1, r2, s2
    Real(Kind = Kr)                                     :: u1, v1

    Real(Kind = Kr), Dimension(:,:), Pointer            :: F, F_Inv
    ! F = linear transformation from the current element to the 3rd simplex
    Type(Vect3D)                                        :: Ve1, Ve2, Ve3, Ve4
    Logical                                             :: GJStat
    Type(Vect3D), Dimension(:), Pointer                 :: GaussP
    Integer, Dimension(:), Pointer                      :: Elements

    ! Make sure the databases are allocated
    If ((.NOT. Associated(Elems_db)) .OR. (.NOT. Associated(Nodes_db))) Then
       Print*, 'One of the required databases is not allocated...'
       Stop
    EndIf
!!$    If (Geom%Numbering /= Numbering_PerNodes) Then
!!$       Write(*,*) 'Error, incompatible numbering scheme for local DoF'
!!$    End If


    !!! Test args here, and allocate elem_list accordingly, if necessary
    If ( Present(Elem) .AND. Present(Elem_List) ) Then
       Print*, 'Error: Init_Gauss_EXO_3D, cannot specify Elem and',           &
            &  'Elem_List at the same time'
       STOP
    ElseIf ( Present(Elem) ) Then
       Allocate (Elements(1))
       Elements = Elem
    ElseIf ( .NOT. Present(Elem_List) ) Then
       Allocate(Elements(Geom%Num_Elems))
       Elements = (/ (iE, iE=1,Geom%Num_Elems) /)
    Else
       Allocate(Elements(Size(ELem_List)))
       Elements = Elem_List
    End If
    
    Do_eNum: Do eNum = 1, Size(Elements)
       iE = Elements(eNum)
       iBlk = Elems_db(iE)%ID_EL
       Elem_Type: Select case (Trim(Geom%Elem_Blk(iBlk)%Type))
       Case('TETRA4', 'TETRA')
          Allocate(F_Inv(3,3))
          Allocate(F(3, 3))
          Ve1 = Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord
          Ve2 = Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord
          Ve3 = Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord
          Ve4 = Nodes_db(Elems_db(iE)%ID_DoF(4))%Coord
          Jac = Vol_Tetra_3D( Ve1, Ve2, Ve3, Ve4 ) * 6.0_Kr
          Jac_Inv = 1.0_Kr / Jac
             
          Nb_BF = 12
          ! Base Functions: (1-X-Y-Z,0,0), (X,0,0), (Y,0,0), (Z,0,0)
          !                 (0,1-X-Y-Z,0), (0,X,0), (0,Y,0), (0,Z,0)
          !                 (0,0,1-X-Y-Z), (0,0,X), (0,0,Y), (0,0,Z)
          OrderG : Select Case(Order)
          Case (1,2,3)
             ! 3rd ordeer cubature on a tetrahedron, cf
             ! J.E. Akin lecture notes / book

             Nb_Gauss = 5
             Allocate (GaussP(Nb_Gauss))

             GaussP%X = (/ 0.25_Kr, 0.5_Kr, InvOf6, InvOf6, InvOf6 /)  
             GaussP%Y = (/ 0.25_Kr, InvOf6, 0.5_Kr, InvOf6, InvOf6 /)
             GaussP%Z = (/ 0.25_Kr, InvOf6, InvOf6, 0.5_Kr, InvOf6 /)
             ! The coordinates of the Gauss Points in the reference element

             Elems_db(iE)%NB_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(Nb_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF,Nb_Gauss))
             Allocate(Elems_db(iE)%Der_BF(NB_BF,Nb_Gauss))
             
             Elems_db(iE)%Gauss_C(1)     = -2.0_Kr / 15.0_Kr * Jac
             Elems_db(iE)%Gauss_C(2:5)   =  3.0_Kr / 40.0_Kr * Jac


             Elems_db(iE)%BF(:,:)%X  = 0.0_Kr
             Elems_db(iE)%BF(:,:)%Y  = 0.0_Kr
             Elems_db(iE)%BF(:,:)%Z  = 0.0_Kr

             Elems_db(iE)%BF(1,:)%X  = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(2,:)%X  = GaussP%X
             Elems_db(iE)%BF(3,:)%X  = GaussP%Y
             Elems_db(iE)%BF(4,:)%X  = GaussP%Z

             Elems_db(iE)%BF(5,:)%Y  = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(6,:)%Y  = GaussP%X
             Elems_db(iE)%BF(7,:)%Y  = GaussP%Y
             Elems_db(iE)%BF(8,:)%Y  = GaussP%Z

             Elems_db(iE)%BF(9,:)%Z  = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(10,:)%Z = GaussP%X
             Elems_db(iE)%BF(11,:)%Z = GaussP%Y
             Elems_db(iE)%BF(12,:)%Z = GaussP%Z

             F_Inv = reshape((/Ve2%X - Ve1%X, Ve2%Y - Ve1%Y, Ve2%Z - Ve1%Z,   &
               &               Ve3%X - Ve1%X, Ve3%Y - Ve1%Y, Ve3%Z - Ve1%Z,   &
               &               Ve4%X - Ve1%X, Ve4%Y - Ve1%Y, Ve4%Z - Ve1%Z/), &
               &             (/3,3/) )
             F = F_Inv
             Call GaussJordan_Inverse(F, GJStat)
             If ( .NOT. GJStat) Then
                Print*, 'Current element is degenerate, cannot compute'
                Print*, 'the linear transformation', iE
                STOP
             End If
             DeAllocate (GaussP)

          Case (4,5)
             
             ! 5th ordeer cubature on a tetrahedron, cf
             ! A.H. Strous: "A fifth degree integration formula for the n-th 
             ! simplex SIAM J. Numer. Anal. Vol 6, No1, March 1969
             Nb_Gauss = 15
             r1 = (7.0_Kr - sqrt(15.0_Kr)) / 34.0_Kr
             s1 = 1.0_Kr - 3.0_Kr * r1
             r2 = (7.0_Kr + sqrt(15.0_Kr)) / 34.0_Kr
             s2 = 1.0_Kr - 3.0_Kr * r2
             
             u1 = (5.0_Kr + sqrt(15.0_Kr))/20.0_Kr
             v1 = (5.0_Kr - sqrt(15.0_Kr))/20.0_Kr 
             Allocate (GaussP(Nb_Gauss))
             GaussP%X = (/ 0.25_Kr, r1, s1, r1, r1, r2, s2, r2, r2,           &
                  &        u1, u1, v1, u1, v1, v1 /)
             GaussP%Y = (/ 0.25_Kr, r1, r1, s1, r1, r2, r2, s2, r2,           &
                  &        u1, v1, u1, v1, u1, v1 /)
             GaussP%Z = (/ 0.25_Kr, r1, r1, r1, s1, r2, r2, r2, s2,           &
                  &        v1, u1, u1, v1, v1, u1 /)
             ! The coordinates of the Gauss Points in the reference element

             Elems_db(iE)%NB_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(Nb_Gauss))
             Allocate(Elems_db(iE)%BF(Nb_BF,Nb_Gauss))
             Allocate(Elems_db(iE)%Der_BF(NB_BF,Nb_Gauss))
             
             Elems_db(iE)%Gauss_C(1)     = 16.0_Kr / 135.0_Kr * Jac * InvOf6
             Elems_db(iE)%Gauss_C(2:5)   = (2665.0_Kr +14.0_Kr*sqrt(15.0_Kr)) &
                  &                       / 37800.0_Kr * Jac * InvOf6
             Elems_db(iE)%Gauss_C(6:9)   = (2665.0_Kr -14.0_Kr*sqrt(15.0_Kr)) &
                  &                       / 37800.0_Kr * Jac * InvOf6
             Elems_db(iE)%Gauss_C(10:15) = 10.0_Kr / 189.0_Kr * Jac * InvOf6


             Elems_db(iE)%BF(:,:)%X  = 0.0_Kr
             Elems_db(iE)%BF(:,:)%Y  = 0.0_Kr
             Elems_db(iE)%BF(:,:)%Z  = 0.0_Kr

             Elems_db(iE)%BF(1,:)%X  = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(2,:)%X  = GaussP%X
             Elems_db(iE)%BF(3,:)%X  = GaussP%Y
             Elems_db(iE)%BF(4,:)%X  = GaussP%Z

             Elems_db(iE)%BF(5,:)%Y  = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(6,:)%Y  = GaussP%X
             Elems_db(iE)%BF(7,:)%Y  = GaussP%Y
             Elems_db(iE)%BF(8,:)%Y  = GaussP%Z

             Elems_db(iE)%BF(9,:)%Z  = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(10,:)%Z = GaussP%X
             Elems_db(iE)%BF(11,:)%Z = GaussP%Y
             Elems_db(iE)%BF(12,:)%Z = GaussP%Z

             F_Inv = reshape((/Ve2%X - Ve1%X, Ve2%Y - Ve1%Y, Ve2%Z - Ve1%Z,   &
               &               Ve3%X - Ve1%X, Ve3%Y - Ve1%Y, Ve3%Z - Ve1%Z,   &
               &               Ve4%X - Ve1%X, Ve4%Y - Ve1%Y, Ve4%Z - Ve1%Z/), &
               &             (/3,3/) )
             F = F_Inv
             Call GaussJordan_Inverse(F, GJStat)
             If ( .NOT. GJStat) Then
                Print*, 'Current element is degenerate, cannot compute'
                Print*, 'the linear transformation', iE
                STOP
             End If
             DeAllocate (GaussP)
          Case Default
             Write(*,*) 'Order non implemented for TETRA4', Order
             STOP
          End Select OrderG
          Elems_db(iE)%Der_BF(:,:)%XX = 0.0_Kr
          Elems_db(iE)%Der_BF(:,:)%XY = 0.0_Kr
          Elems_db(iE)%Der_BF(:,:)%XZ = 0.0_Kr
          Elems_db(iE)%Der_BF(:,:)%YX = 0.0_Kr
          Elems_db(iE)%Der_BF(:,:)%YY = 0.0_Kr
          Elems_db(iE)%Der_BF(:,:)%YZ = 0.0_Kr
          Elems_db(iE)%Der_BF(:,:)%ZX = 0.0_Kr
          Elems_db(iE)%Der_BF(:,:)%ZY = 0.0_Kr
          Elems_db(iE)%Der_BF(:,:)%ZZ = 0.0_Kr

          Elems_db(iE)%Der_BF(1,:)%XX = (-F(1,1)-F(2,1)-F(3,1))
          Elems_db(iE)%Der_BF(1,:)%XY = (-F(1,2)-F(2,2)-F(3,2))
          Elems_db(iE)%Der_BF(1,:)%XZ = (-F(1,3)-F(2,3)-F(3,3))
          
          Elems_db(iE)%Der_BF(2,:)%XX = F(1,1)
          Elems_db(iE)%Der_BF(2,:)%XY = F(1,2)
          Elems_db(iE)%Der_BF(2,:)%XZ = F(1,3)
          
          Elems_db(iE)%Der_BF(3,:)%XX = F(2,1)
          Elems_db(iE)%Der_BF(3,:)%XY = F(2,2)
          Elems_db(iE)%Der_BF(3,:)%XZ = F(2,3)
          
          Elems_db(iE)%Der_BF(4,:)%XX = F(3,1)
          Elems_db(iE)%Der_BF(4,:)%XY = F(3,2)
          Elems_db(iE)%Der_BF(4,:)%XZ = F(3,3)

          Elems_db(iE)%Der_BF(5,:)%YX = (-F(1,1)-F(2,1)-F(3,1))
          Elems_db(iE)%Der_BF(5,:)%YY = (-F(1,2)-F(2,2)-F(3,2))
          Elems_db(iE)%Der_BF(5,:)%YZ = (-F(1,3)-F(2,3)-F(3,3))
          
          Elems_db(iE)%Der_BF(6,:)%YX = F(1,1)
          Elems_db(iE)%Der_BF(6,:)%YY = F(1,2)
          Elems_db(iE)%Der_BF(6,:)%YZ = F(1,3)
          
          Elems_db(iE)%Der_BF(7,:)%YX = F(2,1)
          Elems_db(iE)%Der_BF(7,:)%YY = F(2,2)
          Elems_db(iE)%Der_BF(7,:)%YZ = F(2,3)
          
          Elems_db(iE)%Der_BF(8,:)%YX = F(3,1)
          Elems_db(iE)%Der_BF(8,:)%YY = F(3,2)
          Elems_db(iE)%Der_BF(8,:)%YZ = F(3,3)

          Elems_db(iE)%Der_BF(9,:)%ZX = (-F(1,1)-F(2,1)-F(3,1))
          Elems_db(iE)%Der_BF(9,:)%ZY = (-F(1,2)-F(2,2)-F(3,2))
          Elems_db(iE)%Der_BF(9,:)%ZZ = (-F(1,3)-F(2,3)-F(3,3))
          
          Elems_db(iE)%Der_BF(10,:)%ZX = F(1,1)
          Elems_db(iE)%Der_BF(10,:)%ZY = F(1,2)
          Elems_db(iE)%Der_BF(10,:)%ZZ = F(1,3)
          
          Elems_db(iE)%Der_BF(11,:)%ZX = F(2,1)
          Elems_db(iE)%Der_BF(11,:)%ZY = F(2,2)
          Elems_db(iE)%Der_BF(11,:)%ZZ = F(2,3)
          
          Elems_db(iE)%Der_BF(12,:)%ZX = F(3,1)
          Elems_db(iE)%Der_BF(12,:)%ZY = F(3,2)
          Elems_db(iE)%Der_BF(12,:)%ZZ = F(3,3)
          DeAllocate (F)
          DeAllocate (F_Inv)          
       Case Default
          Write(*,*) 'Element type non implemented yet',                      &
               & Trim(Geom%Elem_Blk(iBlk)%Type)
       End Select Elem_Type
    End Do Do_eNum
    DeAllocate (Elements)
  End Subroutine Init_Gauss_EXO_3D


  Subroutine Init_Gauss_EXO_3D_Elast(Elems_db, Nodes_db, Geom, Order, Elem,   &
       & Elem_List)
    Type (Element3D_Elast), Dimension(:), Pointer       :: Elems_db
    Type (Node3D), Dimension(:), Pointer                :: Nodes_db
    Type (EXO_Geom_Info), intent(IN)                    :: Geom
    Integer                                             :: Order
    Integer, Intent(IN), Optional                       :: Elem
    Integer, Dimension(:), Pointer, Optional            :: Elem_List
    
    Integer                                             :: iE, iBlk
    Integer                                             :: eNum
    Integer                                             :: iBF, iG
    Integer                                             :: NB_BF, NB_Gauss
    Real(Kind = Kr)                                     :: Jac, Jac_Inv

    Real(Kind = Kr)                                     :: r1, s1, r2, s2
    Real(Kind = Kr)                                     :: u1, v1

    Real(Kind = Kr), Dimension(:,:), Pointer            :: F, F_Inv
    ! F = linear transformation from the current element to the 3rd simplex
    Type(Vect3D)                                        :: Ve1, Ve2, Ve3, Ve4
    Logical                                             :: GJStat
    Type(Vect3D), Dimension(:), Pointer                 :: GaussP
    Integer                                             :: Alloc_Stat
    Integer, Dimension(:), Pointer                      :: Elements

    ! Make sure the databases are allocated
    If ((.NOT. Associated(Elems_db)) .OR. (.NOT. Associated(Nodes_db))) Then
       Print*, 'One of the required databases is not allocated...'
       Stop
    EndIf
    Allocate(F_Inv(3,3))
    Allocate(F(3, 3))

    !!! Test args here, and allocate elem_list accordingly, if necessary
    If ( Present(Elem) .AND. Present(Elem_List) ) Then
       Print*, 'Error: Init_Gauss_EXO_3D_Elast, cannot specify Elem and',     &
            &  'Elem_List at the same time'
       STOP
    ElseIf ( Present(Elem) ) Then
       Allocate (Elements(1))
       Elements = Elem
    ElseIf ( .NOT. Present(Elem_List) ) Then
       Allocate(Elements(Geom%Num_Elems))
       Elements = (/ (iE, iE=1,Geom%Num_Elems) /)
    Else
       Allocate(Elements(Size(ELem_List)))
       Elements = Elem_List
    End If
    
    Do_eNum: Do eNum = 1, Size(Elements)
       iE = Elements(eNum)
       iBlk = Elems_db(iE)%ID_EL
       Elem_Type: Select case (Trim(Geom%Elem_Blk(iBlk)%Type))
       Case('TETRA4', 'TETRA')
          Ve1 = Nodes_db(Elems_db(iE)%ID_DoF(1))%Coord
          Ve2 = Nodes_db(Elems_db(iE)%ID_DoF(2))%Coord
          Ve3 = Nodes_db(Elems_db(iE)%ID_DoF(3))%Coord
          Ve4 = Nodes_db(Elems_db(iE)%ID_DoF(4))%Coord
          Jac = Vol_Tetra_3D( Ve1, Ve2, Ve3, Ve4 ) * 6.0_Kr
          Jac_Inv = 1.0_Kr / Jac

          Nb_BF = 12
          ! Base Functions: (1-X-Y-Z,0,0), (X,0,0), (Y,0,0), (Z,0,0)
          !                 (0,1-X-Y-Z,0), (0,X,0), (0,Y,0), (0,Z,0)
          !                 (0,0,1-X-Y-Z), (0,0,X), (0,0,Y), (0,0,Z)
          OrderG : Select Case(Order)
          Case (1,2,3)
             ! 3rd ordeer cubature on a tetrahedron, cf
             ! J.E. Akin lecture notes / book

             Nb_Gauss = 5
             Allocate (GaussP(Nb_Gauss), Stat = Alloc_Stat)
             If (Alloc_Stat /= 0) Then
                Write(*,*) 'Failed to allocate GaussP: iBlk, iE', iBlk, iE
                STOP
             End If

             GaussP%X = (/ 0.25_Kr, 0.5_Kr, InvOf6, InvOf6, InvOf6 /)  
             GaussP%Y = (/ 0.25_Kr, InvOf6, 0.5_Kr, InvOf6, InvOf6 /)
             GaussP%Z = (/ 0.25_Kr, InvOf6, InvOf6, 0.5_Kr, InvOf6 /)
             ! The coordinates of the Gauss Points in the reference element

             Elems_db(iE)%NB_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(Nb_Gauss), Stat = Alloc_Stat)
             If (Alloc_Stat /= 0) Then
                Write(*,*) 'Failed to allocate GaussP: iBlk, iE', iBlk, iE,   &
                     & Alloc_Stat
                STOP
             End If
             Allocate(Elems_db(iE)%BF(Nb_BF,Nb_Gauss), Stat = Alloc_Stat)
             If (Alloc_Stat /= 0) Then
                Write(*,*) 'Failed to allocate Elems_db(iE)%BF: ', iE,        &
                     & Alloc_Stat
                STOP
             End If
             Allocate(Elems_db(iE)%GradS_BF(NB_BF,Nb_Gauss), Stat = Alloc_Stat)
             If (Alloc_Stat /= 0) Then
                Write(*,*) 'Failed to allocate Elems_db(iE)%GradS_BF: iE',    &
                     & iE, Alloc_Stat
                STOP
             End If

             Elems_db(iE)%Gauss_C(1)     = -2.0_Kr / 15.0_Kr * Jac 
             Elems_db(iE)%Gauss_C(2:5)   =  3.0_Kr / 40.0_Kr * Jac 

             Elems_db(iE)%BF(:,:)%X  = 0.0_Kr
             Elems_db(iE)%BF(:,:)%Y  = 0.0_Kr
             Elems_db(iE)%BF(:,:)%Z  = 0.0_Kr

             Elems_db(iE)%BF(1,:)%X  = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(2,:)%X  = GaussP%X
             Elems_db(iE)%BF(3,:)%X  = GaussP%Y
             Elems_db(iE)%BF(4,:)%X  = GaussP%Z

             Elems_db(iE)%BF(5,:)%Y  = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(6,:)%Y  = GaussP%X
             Elems_db(iE)%BF(7,:)%Y  = GaussP%Y
             Elems_db(iE)%BF(8,:)%Y  = GaussP%Z

             Elems_db(iE)%BF(9,:)%Z  = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(10,:)%Z = GaussP%X
             Elems_db(iE)%BF(11,:)%Z = GaussP%Y
             Elems_db(iE)%BF(12,:)%Z = GaussP%Z

             F_Inv = reshape((/Ve2%X - Ve1%X, Ve2%Y - Ve1%Y, Ve2%Z - Ve1%Z,   &
               &               Ve3%X - Ve1%X, Ve3%Y - Ve1%Y, Ve3%Z - Ve1%Z,   &
               &               Ve4%X - Ve1%X, Ve4%Y - Ve1%Y, Ve4%Z - Ve1%Z/), &
               &             (/3,3/) )
             F = F_Inv
             Call GaussJordan_Inverse(F, GJStat)
             If ( .NOT. GJStat) Then
                Print*, 'Current element is degenerate, cannot compute'
                Print*, 'the linear transformation', iE
                STOP
             End If
             DeAllocate (GaussP)
          Case (4,5)
             
             ! 5th ordeer cubature on a tetrahedron, cf
             ! A.H. Strous: "A fifth degree integration formula for the n-th 
             ! simplex SIAM J. Numer. Anal. Vol 6, No1, March 1969
             Nb_Gauss = 15
             r1 = (7.0_Kr - sqrt(15.0_Kr)) / 34.0_Kr
             s1 = 1.0_Kr - 3.0_Kr * r1
             r2 = (7.0_Kr + sqrt(15.0_Kr)) / 34.0_Kr
             s2 = 1.0_Kr - 3.0_Kr * r2
             
             u1 = (5.0_Kr + sqrt(15.0_Kr))/20.0_Kr
             v1 = (5.0_Kr - sqrt(15.0_Kr))/20.0_Kr 
             Allocate (GaussP(Nb_Gauss), Stat = Alloc_Stat)
             If (Alloc_Stat /= 0) Then
                Write(*,*) 'Failed to allocate GaussP: iBlk, iE', iBlk, iE
                STOP
             End If

             GaussP%X = (/ 0.25_Kr, r1, s1, r1, r1, r2, s2, r2, r2,           &
                  &        u1, u1, v1, u1, v1, v1 /)
             GaussP%Y = (/ 0.25_Kr, r1, r1, s1, r1, r2, r2, s2, r2,           &
                  &        u1, v1, u1, v1, u1, v1 /)
             GaussP%Z = (/ 0.25_Kr, r1, r1, r1, s1, r2, r2, r2, s2,           &
                  &        v1, u1, u1, v1, v1, u1 /)
             ! The coordinates of the Gauss Points in the reference element

             Elems_db(iE)%NB_Gauss = Nb_Gauss
             Allocate(Elems_db(iE)%Gauss_C(Nb_Gauss), Stat = Alloc_Stat)
             If (Alloc_Stat /= 0) Then
                Write(*,*) 'Failed to allocate GaussP: iBlk, iE', iBlk, iE,   &
                     & Alloc_Stat
                STOP
             End If
             Allocate(Elems_db(iE)%BF(Nb_BF,Nb_Gauss), Stat = Alloc_Stat)
             If (Alloc_Stat /= 0) Then
                Write(*,*) 'Failed to allocate Elems_db(iE)%BF: ', iE,        &
                     & Alloc_Stat
                STOP
             End If
             Allocate(Elems_db(iE)%GradS_BF(NB_BF,Nb_Gauss), Stat = Alloc_Stat)
             If (Alloc_Stat /= 0) Then
                Write(*,*) 'Failed to allocate Elems_db(iE)%GradS_BF: iE',    &
                     & iE, Alloc_Stat
                STOP
             End If
             
             Elems_db(iE)%Gauss_C(1)     = 16.0_Kr / 135.0_Kr * Jac * InvOf6
             Elems_db(iE)%Gauss_C(2:5)   = (2665.0_Kr +14.0_Kr*sqrt(15.0_Kr)) &
                  &                       / 37800.0_Kr * Jac * InvOf6
             Elems_db(iE)%Gauss_C(6:9)   = (2665.0_Kr -14.0_Kr*sqrt(15.0_Kr)) &
                  &                       / 37800.0_Kr * Jac * InvOf6
             Elems_db(iE)%Gauss_C(10:15) = 10.0_Kr / 189.0_Kr * Jac * InvOf6


             Elems_db(iE)%BF(:,:)%X  = 0.0_Kr
             Elems_db(iE)%BF(:,:)%Y  = 0.0_Kr
             Elems_db(iE)%BF(:,:)%Z  = 0.0_Kr

             Elems_db(iE)%BF(1,:)%X  = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(2,:)%X  = GaussP%X
             Elems_db(iE)%BF(3,:)%X  = GaussP%Y
             Elems_db(iE)%BF(4,:)%X  = GaussP%Z

             Elems_db(iE)%BF(5,:)%Y  = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(6,:)%Y  = GaussP%X
             Elems_db(iE)%BF(7,:)%Y  = GaussP%Y
             Elems_db(iE)%BF(8,:)%Y  = GaussP%Z

             Elems_db(iE)%BF(9,:)%Z  = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
             Elems_db(iE)%BF(10,:)%Z = GaussP%X
             Elems_db(iE)%BF(11,:)%Z = GaussP%Y
             Elems_db(iE)%BF(12,:)%Z = GaussP%Z

             F_Inv = reshape((/Ve2%X - Ve1%X, Ve2%Y - Ve1%Y, Ve2%Z - Ve1%Z,   &
               &               Ve3%X - Ve1%X, Ve3%Y - Ve1%Y, Ve3%Z - Ve1%Z,   &
               &               Ve4%X - Ve1%X, Ve4%Y - Ve1%Y, Ve4%Z - Ve1%Z/), &
               &             (/3,3/) )
             F = F_Inv
             Call GaussJordan_Inverse(F, GJStat)
             If ( .NOT. GJStat) Then
                Print*, 'Current element is degenerate, cannot compute'
                Print*, 'the linear transformation', iE
                STOP
             End If
             DeAllocate (GaussP, Stat = Alloc_Stat)
             If (Alloc_Stat /= 0) Then
                Write(*,*) 'Failed to DeAllocate GaussP'
                STOP
             End If
          Case Default
             Write(*,*) 'Order non implemented for TETRA4', Order
             STOP
          End Select OrderG

          Elems_db(iE)%GradS_BF(:,:)%XX = 0.0_Kr
          Elems_db(iE)%GradS_BF(:,:)%YY = 0.0_Kr
          Elems_db(iE)%GradS_BF(:,:)%ZZ = 0.0_Kr
          Elems_db(iE)%GradS_BF(:,:)%YZ = 0.0_Kr
          Elems_db(iE)%GradS_BF(:,:)%XZ = 0.0_Kr
          Elems_db(iE)%GradS_BF(:,:)%XY = 0.0_Kr

          Elems_db(iE)%GradS_BF(1,:)%XX = (-F(1,1)-F(2,1)-F(3,1)) 
          Elems_db(iE)%GradS_BF(1,:)%XY = (-F(1,2)-F(2,2)-F(3,2)) * 0.5_Kr
          Elems_db(iE)%GradS_BF(1,:)%XZ = (-F(1,3)-F(2,3)-F(3,3)) * 0.5_Kr
          
          Elems_db(iE)%GradS_BF(2,:)%XX = F(1,1) 
          Elems_db(iE)%GradS_BF(2,:)%XY = F(1,2) * 0.5_Kr
          Elems_db(iE)%GradS_BF(2,:)%XZ = F(1,3) * 0.5_Kr
          
          Elems_db(iE)%GradS_BF(3,:)%XX = F(2,1)
          Elems_db(iE)%GradS_BF(3,:)%XY = F(2,2) * 0.5_Kr
          Elems_db(iE)%GradS_BF(3,:)%XZ = F(2,3) * 0.5_Kr
          
          Elems_db(iE)%GradS_BF(4,:)%XX = F(3,1)
          Elems_db(iE)%GradS_BF(4,:)%XY = F(3,2) * 0.5_Kr
          Elems_db(iE)%GradS_BF(4,:)%XZ = F(3,3) * 0.5_Kr

          Elems_db(iE)%GradS_BF(5,:)%XY = (-F(1,1)-F(2,1)-F(3,1)) * 0.5_Kr
          Elems_db(iE)%GradS_BF(5,:)%YY = (-F(1,2)-F(2,2)-F(3,2))
          Elems_db(iE)%GradS_BF(5,:)%YZ = (-F(1,3)-F(2,3)-F(3,3)) * 0.5_Kr
          
          Elems_db(iE)%GradS_BF(6,:)%XY = F(1,1) * 0.5_Kr
          Elems_db(iE)%GradS_BF(6,:)%YY = F(1,2)
          Elems_db(iE)%GradS_BF(6,:)%YZ = F(1,3) * 0.5_Kr
          
          Elems_db(iE)%GradS_BF(7,:)%XY = F(2,1) * 0.5_Kr
          Elems_db(iE)%GradS_BF(7,:)%YY = F(2,2)
          Elems_db(iE)%GradS_BF(7,:)%YZ = F(2,3) * 0.5_Kr
          
          Elems_db(iE)%GradS_BF(8,:)%XY = F(3,1) * 0.5_Kr
          Elems_db(iE)%GradS_BF(8,:)%YY = F(3,2)
          Elems_db(iE)%GradS_BF(8,:)%YZ = F(3,3) * 0.5_Kr

          Elems_db(iE)%GradS_BF(9,:)%XZ = (-F(1,1)-F(2,1)-F(3,1)) * 0.5_Kr
          Elems_db(iE)%GradS_BF(9,:)%YZ = (-F(1,2)-F(2,2)-F(3,2)) * 0.5_Kr
          Elems_db(iE)%GradS_BF(9,:)%ZZ = (-F(1,3)-F(2,3)-F(3,3))
          
          Elems_db(iE)%GradS_BF(10,:)%XZ = F(1,1) * 0.5_Kr
          Elems_db(iE)%GradS_BF(10,:)%YZ = F(1,2) * 0.5_Kr
          Elems_db(iE)%GradS_BF(10,:)%ZZ = F(1,3)
          
          Elems_db(iE)%GradS_BF(11,:)%XZ = F(2,1) * 0.5_Kr
          Elems_db(iE)%GradS_BF(11,:)%YZ = F(2,2) * 0.5_Kr
          Elems_db(iE)%GradS_BF(11,:)%ZZ = F(2,3)
          
          Elems_db(iE)%GradS_BF(12,:)%XZ = F(3,1) * 0.5_Kr
          Elems_db(iE)%GradS_BF(12,:)%YZ = F(3,2) * 0.5_Kr
          Elems_db(iE)%GradS_BF(12,:)%ZZ = F(3,3)
       Case Default
          Write(*,*) 'Element type non implemented yet',                      &
               & Trim(Geom%Elem_Blk(iBlk)%Type)
       End Select Elem_Type
    End Do Do_eNum
    DeAllocate (F)
    DeAllocate (F_Inv)     
    DeAllocate (Elements)
  End Subroutine Init_Gauss_EXO_3D_Elast

  Subroutine InitGauss_3D_Scal_P1(Elem_bd, Nodes_bd, Order)
    Type(Element3D_Scal), Dimension(:), Pointer         :: Elem_bd
    Type(Node3D), Dimension(:), Pointer                 :: Nodes_bd
    Integer                                             :: Order

    Integer(Kind = Ki)                                  :: NE, iE
    Integer(Kind = Ki)                                  :: iGauss, Nb_Gauss
    Integer(Kind = Ki)                                  :: NB_BF
!    Type (Vect2D), Dimension(:),Pointer                 :: GradBF
    Real(Kind = Kr)                                     :: Vol_iE, InvOfVol
    
    Real(Kind = Kr)                                     :: r1, s1, r2, s2
    Real(Kind = Kr)                                     :: u1, v1
    Type(Vect3D), Dimension(:), Pointer                 :: GaussP
    Real(Kind = Kr), Dimension(:,:), Pointer            :: F, F_Inv
    ! F = linear transformation from the current element to the 3rd simplex
    Type(Vect3D)                                        :: Ve1, Ve2, Ve3, Ve4
    Logical                                             :: GJStat


   If ((.NOT. Associated(Elem_bd)) .OR. (.NOT. Associated(Nodes_bd))) Then
      Print*, 'One of the required database is not allocated.'
      Stop
   EndIf
   
   NE = Size(Elem_bd)
   Allocate(F_Inv(3,3))
   Allocate(F(3, 3))

   OrderG : Select Case(Order)
   Case (1,2,3,4,5)
       Nb_BF = 4

      ! 5th ordeer cubature on a tetrahedron, cf
      ! A.H. Strous: "A fifth degree integration formula for the n-th simplex
      ! SIAM J. Numer. Anal. Vol 6, No1, March 1969
       Nb_Gauss = 15
       r1 = (7.0_Kr - sqrt(15.0_Kr)) / 34.0_Kr
       s1 = 1.0_Kr - 3.0_Kr * r1
       r2 = (7.0_Kr + sqrt(15.0_Kr)) / 34.0_Kr
       s2 = 1.0_Kr - 3.0_Kr * r2
       
       u1 = (5.0_Kr + sqrt(15.0_Kr))/20.0_Kr
       v1 = (5.0_Kr - sqrt(15.0_Kr))/20.0_Kr 
       Allocate (GaussP(Nb_Gauss))
       GaussP%X = (/ 0.25_Kr, r1, s1, r1, r1, r2, s2, r2, r2,                 &
            &        u1, u1, v1, u1, v1, v1 /)
       GaussP%Y = (/ 0.25_Kr, r1, r1, s1, r1, r2, r2, s2, r2,                 &
            &        u1, v1, u1, v1, u1, v1 /)
       GaussP%Z = (/ 0.25_Kr, r1, r1, r1, s1, r2, r2, r2, s2,                 &
            &        v1, u1, u1, v1, v1, u1 /)
       Do iE = 1, NE
          Elem_bd(iE)%NB_Gauss = Nb_Gauss
          Allocate(Elem_bd(iE)%Gauss_C(Nb_Gauss))
          Allocate(Elem_bd(iE)%BF(Nb_BF,Nb_Gauss))
          Allocate(Elem_bd(iE)%Grad_BF(NB_BF,Nb_Gauss))

          Ve1 = Nodes_bd(Elem_bd(iE)%ID_DoF(1))%Coord
          Ve2 = Nodes_bd(Elem_bd(iE)%ID_DoF(2))%Coord
          Ve3 = Nodes_bd(Elem_bd(iE)%ID_DoF(3))%Coord
          Ve4 = Nodes_bd(Elem_bd(iE)%ID_DoF(4))%Coord
          Vol_iE = Vol_Tetra_3D( Ve1, Ve2, Ve3, Ve4 )

          Elem_bd(iE)%Gauss_C(1)     = 16.0_Kr / 135.0_Kr * Vol_iE
          Elem_bd(iE)%Gauss_C(2:5)   = (2665.0_Kr + 14.0_Kr * sqrt(15.0_Kr))  &
               &                       / 37800.0_Kr * Vol_iE
          Elem_bd(iE)%Gauss_C(6:9)   = (2665.0_Kr - 14.0_Kr * sqrt(15.0_Kr))  &
               &                       / 37800.0_Kr * Vol_iE
          Elem_bd(iE)%Gauss_C(10:15) = 10.0_Kr / 189.0_Kr * Vol_iE

          InvOfVol = 1.0_Kr / Vol_iE

          Elem_bd(iE)%BF(1,:) = 1.0_Kr - GaussP%X - GaussP%Y - GaussP%Z
          Elem_bd(iE)%BF(2,:) = GaussP%X
          Elem_bd(iE)%BF(3,:) = GaussP%Y
          Elem_bd(iE)%BF(4,:) = GaussP%Z

          F_Inv = reshape ( (/ Ve2%X - Ve1%X, Ve2%Y - Ve1%Y, Ve2%Z - Ve1%Z,   &
               &               Ve3%X - Ve1%X, Ve3%Y - Ve1%Y, Ve3%Z - Ve1%Z,   &
               &               Ve4%X - Ve1%X, Ve4%Y - Ve1%Y, Ve4%Z - Ve1%Z/), &
               &               (/3,3/) )
          F = F_Inv
          Call GaussJordan_Inverse(F, GJStat)
          If ( .NOT. GJStat) Then
             Print*, 'Current element is degenerate, cannot compute'
             Print*, 'the linear transformation', iE
             STOP
          End If

          Elem_bd(iE)%Grad_BF(1,:)%X = (-F(1,1)-F(2,1)-F(3,1)) 
          Elem_bd(iE)%Grad_BF(1,:)%Y = (-F(1,2)-F(2,2)-F(3,2)) 
          Elem_bd(iE)%Grad_BF(1,:)%Z = (-F(1,3)-F(2,3)-F(3,3)) 

          Elem_bd(iE)%Grad_BF(2,:)%X = F(1,1)
          Elem_bd(iE)%Grad_BF(2,:)%Y = F(1,2)
          Elem_bd(iE)%Grad_BF(2,:)%Z = F(1,3)

          Elem_bd(iE)%Grad_BF(3,:)%X = F(2,1)
          Elem_bd(iE)%Grad_BF(3,:)%Y = F(2,2)
          Elem_bd(iE)%Grad_BF(3,:)%Z = F(2,3)

          Elem_bd(iE)%Grad_BF(4,:)%X = F(3,1)
          Elem_bd(iE)%Grad_BF(4,:)%Y = F(3,2)
          Elem_bd(iE)%Grad_BF(4,:)%Z = F(3,3)
       End Do
       
       
       DeAllocate(GaussP)
    Case Default
       Print*, 'This order is not implemented yet...', Order
       Stop
       
    End Select OrderG

    DeAllocate(F, F_Inv)
  End Subroutine InitGauss_3D_Scal_P1


  Subroutine Destroy_Gauss_EXO_3D(Elem_db, Elem, Elem_List)
    Type(Element3D), Dimension(:), Pointer           :: Elem_db
    Integer, Intent(IN), Optional                    :: Elem
    Integer, Dimension(:), Pointer, Optional         :: Elem_List

    Integer, Dimension(:), Pointer                   :: Elements
    Integer                                          :: iE, eNum
    

    If ( Present(Elem) .AND. Present(Elem_List) ) Then
       Print*, 'Error: Destroy_Gauss_EXO_3D, cannot specify Elem and',        &
            &  'Elem_List at the same time'
       STOP
    ElseIf ( Present(Elem) ) Then
       Allocate (Elements(1))
       Elements = Elem
    ElseIf ( .NOT. Present(Elem_List) ) Then
       Allocate(Elements(Size(Elem_db)))
       Elements = (/ (iE, iE=1, Size(Elements)) /)
    Else
       Allocate(Elements(Size(Elem_List)))
       Elements = Elem_List
    End If
    
    Do_eNum: Do eNum = 1, Size(Elements)
       iE = Elements(eNum)
       DeAllocate (Elem_db(iE)%BF)
       DeAllocate (Elem_db(iE)%Der_BF)
       DeAllocate (Elem_db(iE)%Gauss_C)
    End Do Do_eNum

    DeAllocate (Elements)
  End Subroutine Destroy_Gauss_EXO_3D

  Subroutine Destroy_Gauss_EXO_3D_Scal(Elem_db, Elem, Elem_List)
    Type(Element3D_Scal), Dimension(:), Pointer      :: Elem_db
    Integer, Intent(IN), Optional                    :: Elem
    Integer, Dimension(:), Pointer, Optional         :: Elem_List

    Integer, Dimension(:), Pointer                   :: Elements
    Integer                                          :: iE, eNum
    

    If ( Present(Elem) .AND. Present(Elem_List) ) Then
       Print*, 'Error: Destroy_Gauss_EXO_3D_Scal, cannot specify Elem and',   &
            &  'Elem_List at the same time'
       STOP
    ElseIf ( Present(Elem) ) Then
       Allocate (Elements(1))
       Elements = Elem
    ElseIf ( .NOT. Present(Elem_List) ) Then
       Allocate(Elements(Size(Elem_db)))
       Elements = (/ (iE, iE=1, Size(Elements)) /)
    Else
       Allocate(Elements(Size(Elem_List)))
       Elements = Elem_List
    End If
    
    Do_eNum: Do eNum = 1, Size(Elements)
       iE = Elements(eNum)
       DeAllocate (Elem_db(iE)%BF)
       DeAllocate (Elem_db(iE)%Grad_BF)
       DeAllocate (Elem_db(iE)%Gauss_C)
    End Do Do_eNum

    DeAllocate (Elements)
  End Subroutine Destroy_Gauss_EXO_3D_Scal

  Subroutine Destroy_Gauss_EXO_3D_Elast(Elem_db, Elem, Elem_List)
    Type(Element3D_Elast), Dimension(:), Pointer     :: Elem_db
    Integer, Intent(IN), Optional                    :: Elem
    Integer, Dimension(:), Pointer, Optional         :: Elem_List

    Integer, Dimension(:), Pointer                   :: Elements
    Integer                                          :: iE, eNum
    

    If ( Present(Elem) .AND. Present(Elem_List) ) Then
       Print*, 'Error: Destroy_Gauss_EXO_3D_Elast, cannot specify Elem and',  &
            &  'Elem_List at the same time'
       STOP
    ElseIf ( Present(Elem) ) Then
       Allocate (Elements(1))
       Elements = Elem
    ElseIf ( .NOT. Present(Elem_List) ) Then
       Allocate(Elements(Size(Elem_db)))
       Elements = (/ (iE, iE=1, Size(Elements)) /)
    Else
       Allocate(Elements(Size(Elem_List)))
       Elements = Elem_List
    End If
    
    Do_eNum: Do eNum = 1, Size(Elements)
       iE = Elements(eNum)
       DeAllocate (Elem_db(iE)%BF)
       DeAllocate (Elem_db(iE)%GradS_BF)
       DeAllocate (Elem_db(iE)%Gauss_C)
    End Do Do_eNum

    DeAllocate (Elements)
  End Subroutine Destroy_Gauss_EXO_3D_Elast



  Subroutine Destroy_Gauss_EXO_2D(Elem_db, Elem, Elem_List)
    Type(Element2D), Dimension(:), Pointer           :: Elem_db
    Integer, Intent(IN), Optional                    :: Elem
    Integer, Dimension(:), Pointer, Optional         :: Elem_List

    Integer, Dimension(:), Pointer                   :: Elements
    Integer                                          :: iE, eNum
    

    If ( Present(Elem) .AND. Present(Elem_List) ) Then
       Print*, 'Error: Destroy_Gauss_EXO_2D, cannot specify Elem and',      &
            &  'Elem_List at the same time'
       STOP
    ElseIf ( Present(Elem) ) Then
       Allocate (Elements(1))
       Elements = Elem
    ElseIf ( .NOT. Present(Elem_List) ) Then
       Allocate(Elements(Size(Elem_db)))
       Elements = (/ (iE, iE=1, Size(Elements)) /)
    Else
       Allocate(Elements(Size(Elem_List)))
       Elements = Elem_List
    End If
    
    Do_eNum: Do eNum = 1, Size(Elements)
       iE = Elements(eNum)
       DeAllocate (Elem_db(iE)%BF)
       DeAllocate (Elem_db(iE)%Der_BF)
       DeAllocate (Elem_db(iE)%Gauss_C)
    End Do Do_eNum

    DeAllocate (Elements)
  End Subroutine Destroy_Gauss_EXO_2D

  Subroutine Destroy_Gauss_EXO_2D_Scal(Elem_db, Elem, Elem_List)
    Type(Element2D_Scal), Dimension(:), Pointer      :: Elem_db
    Integer, Intent(IN), Optional                    :: Elem
    Integer, Dimension(:), Pointer, Optional         :: Elem_List

    Integer, Dimension(:), Pointer                   :: Elements
    Integer                                          :: iE, eNum
    

    If ( Present(Elem) .AND. Present(Elem_List) ) Then
       Print*, 'Error: Destroy_Gauss_EXO_2D_Scal, cannot specify Elem and',   &
            &  'Elem_List at the same time'
       STOP
    ElseIf ( Present(Elem) ) Then
       Allocate (Elements(1))
       Elements = Elem
    ElseIf ( .NOT. Present(Elem_List) ) Then
       Allocate(Elements(Size(Elem_db)))
       Elements = (/ (iE, iE=1, Size(Elements)) /)
    Else
       Allocate(Elements(Size(Elem_List)))
       Elements = Elem_List
    End If
    
    Do_eNum: Do eNum = 1, Size(Elements)
       iE = Elements(eNum)
       DeAllocate (Elem_db(iE)%BF)
       DeAllocate (Elem_db(iE)%Grad_BF)
       DeAllocate (Elem_db(iE)%Gauss_C)
    End Do Do_eNum

    DeAllocate (Elements)
  End Subroutine Destroy_Gauss_EXO_2D_Scal

  Subroutine Destroy_Gauss_EXO_2D_Elast(Elem_db, Elem, Elem_List)
    Type(Element2D_Elast), Dimension(:), Pointer     :: Elem_db
    Integer, Intent(IN), Optional                    :: Elem
    Integer, Dimension(:), Pointer, Optional         :: Elem_List

    Integer, Dimension(:), Pointer                   :: Elements
    Integer                                          :: iE, eNum
    

    If ( Present(Elem) .AND. Present(Elem_List) ) Then
       Print*, 'Error: Destroy_Gauss_EXO_2D_Elast, cannot specify Elem and',  &
            &  'Elem_List at the same time'
       STOP
    ElseIf ( Present(Elem) ) Then
       Allocate (Elements(1))
       Elements = Elem
    ElseIf ( .NOT. Present(Elem_List) ) Then
       Allocate(Elements(Size(Elem_db)))
       Elements = (/ (iE, iE=1, Size(Elements)) /)
    Else
       Allocate(Elements(Size(Elem_List)))
       Elements = Elem_List
    End If
    
    Do_eNum: Do eNum = 1, Size(Elements)
       iE = Elements(eNum)
       DeAllocate (Elem_db(iE)%BF)
       DeAllocate (Elem_db(iE)%GradS_BF)
       DeAllocate (Elem_db(iE)%Gauss_C)
    End Do Do_eNum

    DeAllocate (Elements)
  End Subroutine Destroy_Gauss_EXO_2D_Elast

End Module m_MEF_Gauss
