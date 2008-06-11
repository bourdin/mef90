Program TestNewLayout

   Use m_MEF90
   Implicit NONE
    
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscviewer.h"
#include "include/finclude/petscviewer.h90"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscao.h"

   Type (EXO_Geom_Info)                          :: Geom, MyGeom
   Type (Layout_Info)                            :: Layout
   Integer, Dimension(:,:), Pointer              :: MyConnect
   Type(Vect2D), Dimension(:), Pointer           :: MyCoord
   Type (Element2D_Scal), Dimension(:), Pointer  :: Elem
   Type (Node2D), Dimension(:), Pointer          :: Node
   Vec                                           :: VN_Dist, VN_IO, VNG_Dist, VNG_Local
   Vec                                           :: VE_Dist, VE_IO
   Real(Kind = Kr), Dimension(:), Pointer        :: VN_Ptr
   
   Character(len=128)                            :: CharBuffer
   Integer                                       :: iErr, i
   PetscTruth                                    :: HasF
   PetscReal                                     :: Val
   
   PetscViewer                                   :: viewer
   
   
   Call MEF90_Initialize()
   Call PetscOptionsGetString(PETSC_NULL_CHARACTER, '-f', Geom%filename, HasF, iErr)    
   If (.NOT. HasF) Then
      Call PetscPrintf(PETSC_COMM_WORLD, "No input file given\n"c, iErr)
      Call MEF90_Finalize()
      STOP
   End If

   Geom%Numbering = Numbering_PerNodes
   Geom%Comm = PETSC_COMM_WORLD
   Call Read_EXO_Geom_Info(Geom, MyGeom)
   
   Write(CharBuffer,101) MEF90_MyRank
 101 Format('.',I3.3)  
   Call Init_Layout_TRI3(Geom, MyGeom, Layout, Node, Elem)
   MyGeom%filename = trim(Geom%filename)//adjustL(CharBuffer)


   Write(CharBuffer,*) 'Total number of dof:        ', Geom%Num_Nodes, '\n'c
   Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
   Write(CharBuffer,*) 'Local number of dof/ghosts: ', Layout%Num_Local_dof, Layout%Num_Ghost_dof, '\n'c
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr)
   
   Write(CharBuffer,*) 'Total number of elements:   ', Geom%Num_Elems, '\n'c
   Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
   Write(CharBuffer,*) 'Local number of elements:    ', Layout%Num_Local_elems, '\n'c
   Call PetscSynchronizedPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
   Call PetscSynchronizedFlush(PETSC_COMM_WORLD, iErr)


   Call Show_EXO_Geom_Info(MyGeom, 200+MEF90_MyRank)
   Call Write_EXO_Geom_Info(MyGeom)
!   Write(*,*) "MyGeom%filename: ", MyGeom%filename
   MyGeom%Numbering = Numbering_PerNodes
   Call Write_EXO_Node_Coord(MyGeom, Node, 1)
   Call Write_EXO_Connect(MyGeom, Elem)
   
   Call VecCreateMPI(Geom%Comm, Layout%num_local_dof, Geom%num_nodes, VN_Dist, iErr)
   Call VecCreateMPI(Geom%Comm, Layout%num_local_elems, Geom%num_elems, VE_Dist, iErr)
   
   If (MEF90_MyRank == 0) Then
      Call VecCreateMPI(Geom%Comm, Geom%num_nodes, Geom%num_nodes, VN_IO, iErr)
      Call VecCreateMPI(Geom%Comm, Geom%num_elems, Geom%num_elems, VE_IO, iErr)
   Else
      Call VecCreateMPI(Geom%Comm, 0, Geom%num_nodes, VN_IO, iErr)
      Call VecCreateMPI(Geom%Comm, 0, Geom%num_elems, VE_IO, iErr)
   End If   
   
   Call VecCreateGhost(Geom%Comm, Layout%num_local_dof, Geom%num_nodes, Layout%num_ghost_dof, Layout%ghost_dof-1, VNG_Dist, iErr)
	Call VecGhostGetLocalForm(VNG_Dist, VNG_Local, iErr)
	
	
	
   Call VecGetArrayF90(VNG_Local, VN_Ptr, iErr)
   VN_Ptr = MEF90_MyRank
   Do i = 1, Layout%num_local_dof
      VN_Ptr(i) = i!Node(i)%Coord%X**2 + Node(i)%Coord%Y**2
   End Do
   Write(700+MEF90_MyRank, *) 'Before Ghost Update'
   Write(700+MEF90_MyRank, *) VN_Ptr
   Call VecRestoreArrayF90(VNG_Local, VN_Ptr, iErr)
   
   Call VecGhostUpdateBegin(VNG_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
   Call VecGhostUpdateEnd  (VNG_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)

   Call VecGetArrayF90(VNG_Local, VN_Ptr, iErr)
   Write(700+MEF90_MyRank, *) 'After Ghost Update'
   Write(700+MEF90_MyRank, *) VN_Ptr
   Call VecRestoreArrayF90(VNG_Local, VN_Ptr, iErr)

   Call Write_EXO_Result_Nodes(Geom, Layout, 4, 2, VNG_Dist)
   Call Write_EXO_Result_Nodes(MyGeom, Layout, 4, 2, VNG_Local)
   
   Write(700+MEF90_MyRank, *) "===== COORDINATES"
   Do i = 1, MyGeom%Num_Nodes
      Write(700+MEF90_MyRank, *) i, Node(i)%Coord%X, Node(i)%Coord%Y
   End Do
   Call MEF90_Finalize()
   STOP

	
	
	val = MEF90_MyRank + 1.0_Kr
   Call VecSet(VNG_Local, val, iErr)
   Call PetscPrintf(PETSC_COMM_WORLD, "VNG_Dist\n"c, iErr)
   Call VecView(VNG_Dist, PETSC_VIEWER_STDOUT_WORLD, iErr)
   
   Call VecGetArrayF90    (VNG_Local, VN_Ptr, iErr)
   Write(500+MEF90_MyRank, *) "Before Ghost Update"
   Do i = 1, Size(VN_Ptr)
      Write(500+MEF90_MyRank, *) i, VN_Ptr(i)
   End Do
   Call VecRestoreArrayF90(VNG_Local, VN_Ptr, iErr)
   
   Call VecGhostUpdateBegin(VNG_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)
   Call VecGhostUpdateEnd  (VNG_Dist, INSERT_VALUES, SCATTER_FORWARD, iErr)

   Call VecGetArrayF90    (VNG_Local, VN_Ptr, iErr)
   Write(500+MEF90_MyRank, *) "After Ghost Update"
   Do i = 1, Size(VN_Ptr)
      Write(500+MEF90_MyRank, *) i, VN_Ptr(i)
   End Do
   Call VecRestoreArrayF90(VNG_Local, VN_Ptr, iErr)
   Write(500+MEF90_MyRank, * ) "Layout%Ghost_dof", Layout%Ghost_dof 
   
   Call PetscFinalize(iErr)
   STOP


	Call VecSetLocalToGlobalMapping(VE_Dist, Layout%Mapping_E, iErr)
	Do i = 1, MyGeom%num_elems
	  Val = 100.0_Kr * (MEF90_MyRank + 1.0_Kr) + i - 1.0_Kr
      Call VecSetValueLocal(VE_Dist, i-1, Val, INSERT_VALUES, iErr)
	End Do
	Call VecAssemblyBegin(VE_Dist, iErr)
	Call VecAssemblyEnd  (VE_Dist, iErr)
!	Call PetscPrintf(PETSC_COMM_WORLD, "================================ VE_Dist\n"c, iErr) 
!	Call VecView(VE_Dist, PETSC_VIEWER_STDOUT_WORLD, iErr)
	
   Call VecScatterBegin(Layout%ToIOSeq_E, VE_Dist, VE_IO, INSERT_VALUES, SCATTER_FORWARD, iErr) 
   Call VecScatterEnd  (Layout%ToIOSeq_E, VE_Dist, VE_IO, INSERT_VALUES, SCATTER_FORWARD, iErr) 
!   Call PetscPrintf(PETSC_COMM_WORLD, "================================ VE_IO\n"c, iErr) 
!	Call VecView(VE_IO, PETSC_VIEWER_STDOUT_WORLD, iErr)

	Call VecSetLocalToGlobalMapping(VN_Dist, Layout%Mapping_N, iErr)
	Do i = 1, MyGeom%num_nodes
	  Val = 100.0_Kr * (MEF90_MyRank + 1.0_Kr) + i - 1.0_Kr
	  Call VecSetValueLocal(VN_Dist, i-1, Val, Add_VALUES, iErr)
	End Do
	Call VecAssemblyBegin(VN_Dist, iErr)
	Call VecAssemblyEnd  (VN_Dist, iErr)

	Call VecAssemblyBegin(VN_Dist, iErr)
	Call VecAssemblyEnd  (VN_Dist, iErr)
!	Call PetscPrintf(PETSC_COMM_WORLD, "================================ VN_Dist\n"c, iErr) 
!	Call VecView(VN_Dist, PETSC_VIEWER_STDOUT_WORLD, iErr)

   Call VecScatterBegin(Layout%ToIOSeq_N, VN_Dist, VN_IO, ADD_VALUES, SCATTER_FORWARD, iErr) 
   Call VecScatterEnd  (Layout%ToIOSeq_N, VN_Dist, VN_IO, ADD_VALUES, SCATTER_FORWARD, iErr) 
!	Call PetscPrintf(PETSC_COMM_WORLD, "================================ VN_IO\n"c, iErr) 
!	Call VecView(VN_IO, PETSC_VIEWER_STDOUT_WORLD, iErr)
	
   Call VecGetArrayF90(VN_IO, VN_Ptr, iErr)
   If (MEF90_MyRank == 0) Then
      Do i = 1, Geom%Num_Nodes
         VN_Ptr(i) = 1.0_Kr + i - 1.0_Kr
      End Do
   End If
   Call VecRestoreArrayF90(VN_IO, VN_Ptr, iErr)
!	Call PetscPrintf(PETSC_COMM_WORLD, "================================ VN_IO\n"c, iErr) 
!	Call VecView(VN_IO, PETSC_VIEWER_STDOUT_WORLD, iErr)
   Call VecScatterBegin(Layout%ToIOSeq_N, VN_IO, VN_Dist, INSERT_VALUES, SCATTER_REVERSE, iErr) 
   Call VecScatterEnd  (Layout%ToIOSeq_N, VN_IO, VN_Dist, INSERT_VALUES, SCATTER_REVERSE, iErr) 
!	Call PetscPrintf(PETSC_COMM_WORLD, "================================ VN_IO\n"c, iErr) 
!	Call VecView(VN_IO, PETSC_VIEWER_STDOUT_WORLD, iErr)
!	Call PetscPrintf(PETSC_COMM_WORLD, "================================ VN_Dist\n"c, iErr) 
!	Call VecView(VN_Dist, PETSC_VIEWER_STDOUT_WORLD, iErr)

   Call MEF90_Finalize()
 !  STOP 

End Program TestNewLayout