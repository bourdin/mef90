Program TestSieve

   Use m_MEF90
   Use m_TestSieve
   Implicit NONE
    
#include "finclude/petsc.h"
#include "finclude/petscvec.h"
#include "finclude/petscviewer.h"
#include "finclude/petscviewer.h90"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscmat.h90"
#include "finclude/petscksp.h"
#include "finclude/petscpc.h"
#include "finclude/petscao.h"
#include "finclude/petscmesh.h"
#include "finclude/petscmesh.h90"

   Type (MeshTopology_Info)                     :: MeshTopology
   Type (EXO_Info)                              :: EXO
   Type (Element2D_Scal), Dimension(:), Pointer :: Elem2DA
   Type (Vect3D), Dimension(:), Pointer         :: Coords
   PetscReal, Dimension(:,:), Pointer           :: Vertices
   
   PetscReal                                    :: ObjectiveFunction
   Vec                                          :: U, F, Gradient
   PetscReal, Dimension(:), Pointer             :: U_Ptr, F_Ptr
   Mat                                          :: Hessian
   PetscTruth                                   :: HasF
   Integer                                      :: iErr
   Integer                                      :: iBlk, iELoc, iE, iS
   Character(len=256)                           :: CharBuffer
   Integer, Dimension(:,:), Pointer             :: Tmp_Connect
   
   Integer                                      :: vers
   Integer, Parameter                           :: exo_cpu_ws = 8
   Integer, Parameter                           :: exo_io_ws = 8

     
   Call MEF90_Initialize()
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-f', EXO%filename, HasF, iErr)    
   If (.NOT. HasF) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file given\n", iErr)
      Call MEF90_Finalize()
      STOP
   End If



   EXO%Comm = PETSC_COMM_WORLD
   
   Call Read_MeshTopology_Info_EXO(MeshTopology, Coords, Elem2DA, EXO)
   
   MeshTopology%Elem_Blk%Elem_Type    = MEF90_P1_Lagrange
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Call Init_Elem_Blk_Info(MeshTopology%Elem_Blk(iBlk), MeshTopology%num_dim)
   End Do
   Call Show_MeshTopology_Info(MeshTopology)

   Allocate (Vertices(2,3))

   !!! Initialize the element   
   Do iBlk = 1, MeshTopology%Num_Elem_Blks
      Do iELoc = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
         iE = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iELoc)
         Vertices(1,:) = Coords(Elem2DA(iE)%ID_DoF(:))%X
         Vertices(2,:) = Coords(Elem2DA(iE)%ID_DoF(:))%Y
         Call Init_Element(Elem2DA(iE), Vertices, 4, MeshTopology%Elem_Blk(iBlk)%Elem_Type)
      End Do
   End Do

!   Call Show_Elem2D_Scal(Elem2DA)

   Call VecCreateSeq(PETSC_COMM_WORLD, MeshTopology%Num_Vert, U, iErr)
   Call VecCreateSeq(PETSC_COMM_WORLD, MeshTopology%Num_Vert, F, iErr)

   Call VecSet(U, 1.0_Kr, iErr)
   Call VecSet(F, 0.0_Kr, iErr)

   Call VecGetArrayF90(F, F_Ptr, iErr)
   Call VecGetArrayF90(U, U_Ptr, iErr)

   Do iS = 1, Size(U_Ptr)
      U_Ptr(iS) = Coords(iS)%X
!      F_Ptr(iS) = Coords(iS)%Y
   End Do

   Call VecRestoreArrayF90(U, U_Ptr, iErr)
   Call VecRestoreArrayF90(F, F_Ptr, iErr)

   Call FormObjectiveFunction(ObjectiveFunction, MeshTopology, Elem2DA, U, F)

   Write(*,*) 'Objective Function: ', ObjectiveFunction
   Call MEF90_Finalize()

 Contains
   Subroutine Show_Elem2D_Scal(dElems, Unit)
      Type (Element2D_Scal), DImension(:), Pointer   :: dElems
      Integer, Optional                              :: Unit
      
      Integer                                        :: iE, Nb_Gauss, Nb_DoF, iG, iDoF
      Integer                                        :: IO_Unit
    
      If ( Present(Unit) ) Then
         IO_Unit = Unit
      Else
         IO_Unit = 6 !!! STDOUT
      End If 
 
      
      Write(IO_Unit, 100) Size(dElems)
      
      Do iE = 1, Size(dElems)
         Write(IO_Unit, 101) iE
         Nb_DoF   = size(dElems(iE)%BF,1)
         Nb_Gauss = size(dElems(iE)%BF,2)
         Write(IO_Unit, 102) Nb_DoF
         Write(IO_Unit, 103) Nb_Gauss

         Do iDoF = 1, Nb_DoF
            Write(IO_Unit, 201) iDoF
            Write(IO_Unit, 200, advance = 'no') 'BF'
            Write(IO_Unit, *) dElems(iE)%BF(iDoF, :)
            Write(IO_Unit, 200, advance = 'no') 'Grad_BF%X'
            Write(IO_Unit, *) dElems(iE)%Grad_BF(iDoF, :)%X
            Write(IO_Unit, 200, advance = 'no') 'Grad_BF%Y'
            Write(IO_Unit, *) dElems(iE)%Grad_BF(iDoF, :)%Y
         End Do
      End Do
100 Format('    Number of elements =================== ', I9)
101 Format('*** Element  ', I9)
102 Format('    Nb_DoF   ', I9)
103 Format('    Nb_Gauss ', I9)
200 Format(A)
201 Format('    *** DoF  ', I9)
   End Subroutine Show_Elem2D_Scal
   
   
   Subroutine Read_MeshTopology_Info_EXO(dMeshTopology, Coords, Elem2DA, dEXO)
      Type (MeshTopology_Info)                     :: dMeshTopology
      Type (Vect3D), Dimension(:), Pointer         :: Coords
      Type (Element2D_Scal), Dimension(:), Pointer :: Elem2DA
      Type (EXO_Info)                              :: dEXO
      Integer                                      :: iErr, iBlk, iSet, Offset, i
      Character(len=256)                           :: CharBuffer
      
      PetscReal, Dimension(:,:), Pointer           :: array
      PetscInt, Dimension(:,:), Pointer            :: arrayCon
      Integer                                      :: EXO_DummyInteger
      Character(len=MXSTLN)                        :: EXO_DummyStr   
      Integer                                      :: embedDim
      Integer                                      :: iNode, iElem
      Mesh                                         :: mesh

      ! Open File
      call MeshCreateExodus(PETSC_COMM_WORLD, dEXO%filename, mesh, ierr)

      ! Read Global Geometric Parameters
      call MeshExodusGetInfo(mesh, dMeshTopology%Num_Dim, dMeshTopology%Num_Vert, dMeshTopology%Num_Elems, &
           & dMeshTopology%Num_Elem_Blks, dMeshTopology%Num_Node_Sets, iErr)
      ! Read Elem blocks informations
      Allocate(dMeshTopology%Elem_blk(dMeshTopology%Num_Elem_blks))
      CharBuffer = 'CellBlocks'
      If (dMeshTopology%Num_Elem_blks > 0) Then
         Do iBlk = 1, dMeshTopology%Num_Elem_Blks
            call MeshGetStratumSize(mesh, CharBuffer, iBlk-1, dMeshTopology%elem_blk(iBlk)%Num_Elems, ierr)
            write(6,*) 'Number of elements in block',iBlk-1,dMeshTopology%Elem_blk(iBlk)%Num_Elems
            Allocate(dMeshTopology%Elem_blk(iBlk)%Elem_ID(dMeshTopology%elem_blk(iBlk)%Num_Elems))
            call MeshGetStratum(mesh, CharBuffer, iBlk-1, dMeshTopology%Elem_blk(iBlk)%Elem_ID, ierr)
            dMeshTopology%Elem_blk(iBlk)%Elem_ID = dMeshTopology%Elem_blk(iBlk)%Elem_ID + 1
         End Do
      End If
      
      ! Read Node sets informations
      Allocate (dMeshTopology%Node_Set(dMeshTopology%Num_Node_Sets))
      CharBuffer = 'VertexSets'
      If (dMeshTopology%Num_Node_Sets > 0) Then
         Do iSet = 1, dMeshTopology%Num_node_sets
            call MeshGetStratumSize(mesh, CharBuffer, iSet-1, dMeshTopology%Node_Set(iSet)%Num_Nodes, ierr)
            write(6,*) 'Number of nodes in set',iSet-1,dMeshTopology%Node_set(iSet)%Num_Nodes
            Allocate(dMeshTopology%Node_Set(iSet)%Node_ID(dMeshTopology%Node_Set(iSet)%Num_Nodes))
            call MeshGetStratum(mesh, CharBuffer, iSet-1, dMeshTopology%Node_Set(iSet)%Node_ID, ierr)
            dMeshTopology%Node_Set(iSet)%Node_ID = dMeshTopology%Node_Set(iSet)%Node_ID - MeshTopology%Num_Elems
         End Do
      End If

      ! Read the vertices coordinates
      Allocate(Coords(MeshTopology%Num_Vert))
      call MeshGetCoordinatesF90(mesh, array, iErr)
      embedDim = size(array,2)
      Do iNode = 1, dMeshTopology%Num_Vert
         Coords(iNode)%X = array(iNode,1)
         If (embedDim > 1) Then
            Coords(iNode)%Y = array(iNode,2)
            If (embedDim > 2) Then
               Coords(iNode)%Z = array(iNode,3)
            Else
               Coords(iNode)%Z = 0.0
            End If
         Else
            Coords(iNode)%Y = 0.0
            Coords(iNode)%Z = 0.0
         End If
         write(6,*) 'Node',iNode,Coords(iNode)
      End Do
      call MeshRestoreCoordinatesF90(mesh, array, iErr)

      ! Read the connectivity table
      Allocate(Elem2DA(MeshTopology%Num_Elems))
      call MeshGetElementsF90(mesh, arrayCon, iErr)
      Do iBlk = 1, MeshTopology%Num_Elem_Blks
         Do iE = 1, MeshTopology%Elem_Blk(iBlk)%Num_Elems
            iElem = MeshTopology%Elem_Blk(iBlk)%Elem_ID(iE)
            Allocate(Elem2DA(iElem)%ID_DoF(3))
            !Elem2DA(iElem)%ID_DoF(:) = arrayCon(iE,:)
            Do i = 1, 3
               write(6,*) 'arrayCon',iElem,i,arrayCon(iElem,i)
               Elem2DA(iElem)%ID_DoF(i) = arrayCon(iElem,i)
            End Do
            write(6,*) 'Elem',iElem,Elem2DA(iElem)%ID_DoF
         End Do
      End Do
      call MeshRestoreElementsF90(mesh, arrayCon, iErr)

      dEXO%exoid = 0
    End Subroutine Read_MeshTopology_Info_EXO

   Subroutine Show_MeshTopology_Info(dMeshTopology, Unit)
      Type (MeshTopology_Info), Intent(IN)           :: dMeshTopology
      Integer, Optional                              :: Unit
      
      Integer                                        :: IO_Unit, i
    
      If ( Present(Unit) ) Then
         IO_Unit = Unit
      Else
         IO_Unit = 6 !!! STDOUT
      End If  


      Write(IO_Unit, 103) MeshTopology%num_dim
      Write(IO_Unit, 104) MeshTopology%num_Vert
      Write(IO_Unit, 105) MeshTopology%num_elems
      Write(IO_Unit, 106) MeshTopology%Num_Elem_blks
      Write(IO_Unit, 107) MeshTopology%Num_Node_Sets
      Write(IO_Unit, 108) MeshTopology%Num_Side_Sets
      
      Write(IO_Unit, *)
      Write(IO_Unit, 200)
      Write(IO_Unit, 201) MeshTopology%num_elem_blks
      Do i = 1, dMeshTopology%Num_Elem_blks
         Write(IO_Unit, 203) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%Num_Elems
         Write(IO_Unit, 204) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%Elem_Type
         Write(IO_Unit, 205) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%DoF_Location
!         Write(IO_Unit, 206) dMeshTopology%Elem_Blk(i)%ID, dMeshTopology%Elem_blk(i)%Nb_Gauss
         Write(IO_Unit, 207, advance = 'no') dMeshTopology%Elem_Blk(i)%ID
         Write(IO_Unit, *) dMeshTopology%Elem_blk(i)%Elem_ID
      End Do
      
      
      Write(IO_Unit, *)
      Write(IO_Unit, 300)
      Write(IO_Unit, 301) dMeshTopology%num_node_sets
      Do i = 1, dMeshTopology%num_node_sets
         Write(IO_Unit, 302) dMeshTopology%Node_Set(i)%ID, dMeshTopology%Node_Set(i)%Num_Nodes
         Write(IO_Unit, 303, advance = 'no') dMeshTopology%Node_Set(i)%ID
         Write(IO_Unit, *) dMeshTopology%Node_Set(i)%Node_ID
      End Do
      
      
      Write(IO_Unit, *)
      Write(IO_Unit, 400)
      Write(IO_Unit, 401) dMeshTopology%num_side_sets
      
100 Format('*** GLOBAL INFORMATIONS ***')
103 Format('    Number of dimensions ============ ', I6)
104 Format('    Number of vertices ============== ', I6)
105 Format('    Number of elements ============== ', I6)
106 Format('    Number of elements blocks ======= ', I6)
107 Format('    Number of node sets ============= ', I6)
108 Format('    Number of side sets ============= ', I6)

200 Format('*** ELEMENT BLOCKS ***')
201 Format('    Number of blocks ================ ', I4)
202 Format('    Block ', I3, ' Elements type ========= ', A)
203 Format('    Block ', I3, ' Number of elements ==== ', I4)
204 Format('    Block ', I3, ' Element type ========== ', I4)
205 Format('    Block ', I3, ' DoF location ========== ', 4(I4, ' '))
!206 Format('    Block ', I3, ' Nb_Gauss==== ========== ', I4)
207 Format('    Block ', I3, ' IDs: ')

300 Format('*** NODE SETS ***')
301 Format('    Number of sets ================== ', I4)
302 Format('    Set ', I3, ' Number of nodes ========= ', I4)
303 Format('    Set ', I3, ' IDs: ')
304 Format('    Set ', I3, ' Number of dist. factors = ', I4)

400 Format('*** SIDE SETS ***')
401 Format('    Number of side sets ============= ', I4)
    
   End Subroutine Show_MeshTopology_Info

End Program TestSieve
