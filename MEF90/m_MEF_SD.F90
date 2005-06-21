Module m_MEF_SD

  Use m_MEF_Types
!  Use m_MEF_Communication
  IMPLICIT NONE
  Private

#include "include/finclude/petsc.h"
#include "include/finclude/petscsys.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscis.h"
#include "include/finclude/petscis.h90"
#include "include/finclude/petscao.h"

  Public :: SD_Info
  Public :: Init_SD_NoOvlp

  Interface Init_SD_NoOvlp
     Module Procedure Init_SD_NoOvlp_2DScal, Init_SD_NoOvlp_2DElast,          &
          & Init_SD_NoOvlp_2D, Init_SD_NoOvlp_3DScal, Init_SD_NoOvlp_3DElast, &
          & Init_SD_NoOvlp_3D
  End Interface

  Type SD_Info
     Integer                                         :: Num_Elems
     Integer                                         :: Num_GhostElems
     Integer                                         :: Num_Nodes
     Integer                                         :: Num_GhostNodes

     Integer, Dimension(:), Pointer                  :: Elem
     Integer, Dimension(:), Pointer                  :: GhostElem
     Integer, Dimension(:), Pointer                  :: Node
     Integer, Dimension(:), Pointer                  :: GhostNode

     Logical, Dimension(:), Pointer                  :: IsLocal_Elem
     Logical, Dimension(:), Pointer                  :: IsLocal_GhostElem
     Logical, Dimension(:), Pointer                  :: IsLocal_Node
     Logical, Dimension(:), Pointer                  :: IsLocal_GhostNode

     VecScatter                                      :: ToMaster
 
     IS                                              :: EXO_IS
     IS                                              :: PETSc_IS

     AO                                              :: EXO_AO
     AO                                              :: Loc_AO
  End Type SD_Info


  
  Integer, Parameter                              :: METIS_etype_tri  = 1
  Integer, Parameter                              :: METIS_etype_tet  = 2
  Integer, Parameter                              :: METIS_etype_hex  = 3
  Integer, Parameter                              :: METIS_etype_quad = 4

  Integer, Parameter                              :: METIS_numflag = 1

  Contains 
    Subroutine Init_SD_NoOvlp_2DScal(Geom, Elem_db, SD_db, Vect_Dist,         &
         & Vect_Local, Vect_Master)
      Type (EXO_Geom_Info), Intent(IN)                 :: Geom
      Type (Element2D_Scal), Dimension(:), Pointer     :: Elem_db
      Type (SD_Info), Intent(OUT)                      :: SD_db
      Vec, Intent(INOUT)                               :: Vect_Dist
      Vec, Intent(INOUT)                               :: Vect_Local
      Vec, Intent(INOUT)                               :: Vect_Master


      Integer, Dimension(:), Pointer                   :: METIS_epart
      Integer, Dimension(:), Pointer                   :: METIS_npart
      Integer, Dimension(:), Pointer                   :: METIS_elmnts
      Integer                                          :: METIS_edgecut
      
      Integer                                          :: MyRank
      Integer                                          :: NumProcs
      Integer                                          :: Nb_DoF
      Integer                                          :: iE, iELoc
      Integer                                          :: iN, iNLoc
      Integer                                          :: iErr

      Integer, Dimension(:), Pointer                   :: Tmp_Idx, Tmp_Idx2
      Integer                                          :: MyIdxMin, MyIdxMax
      Integer                                          :: i, iG, inG

      
      ! No overlap => no ghost elements
      Allocate(SD_db%IsLocal_Node(Geom%Num_Nodes))
      Allocate(SD_db%IsLocal_GhostNode(Geom%Num_Nodes))
      Allocate(SD_db%IsLocal_Elem(Geom%Num_Elems))
      Allocate(SD_db%IsLocal_GhostElem(Geom%Num_Elems))
      SD_db%IsLocal_Node      = .FALSE.
      SD_db%IsLocal_GhostNode = .FALSE.
      SD_db%IsLocal_Elem      = .FALSE.
      SD_db%IsLocal_GhostElem = .FALSE.

      Nb_DoF = 3

      Call MPI_COMM_RANK(MPI_COMM_WORLD, MyRank, iErr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, NumProcs, iErr)
 
      ! Step 1: Call METIS and generate the subdomains
      Allocate (METIS_elmnts(Geom%Num_Elems*Nb_DoF))
      Allocate(METIS_epart(Geom%Num_Elems))
      Allocate(METIS_npart(Geom%Num_Nodes))

      ! This assumes that all elements in the mesh are triangular.
      ! Unfortunately, METIS doesn't seem to be able to deal with 
      ! mixed elements meshes.
      
      Allocate (METIS_elmnts(Geom%Num_Elems*Nb_DoF))
      Do iE = 1, Geom%Num_Elems
         METIS_elmnts((iE-1)*Nb_DoF +1:iE*Nb_DoF) = Elem_db(iE)%ID_DoF
      End Do

!      Write(1000, *) METIS_elmnts
!
      If (NumProcs == 1) Then
         METIS_epart = 1
         METIS_npart = 1
      Else
         Call METIS_PartMeshNodal(Geom%Num_Elems, Geom%Num_Nodes,             &
              & METIS_elmnts, METIS_etype_tri, METIS_numflag, NumProcs,       &
              & METIS_edgecut, METIS_epart, METIS_npart)
      End If

      ! Step 2: Allocate the fields in SD_Info
      SD_db%Num_Nodes = Count(METIS_npart == MyRank+1)
      Allocate(SD_db%Node(SD_db%Num_Nodes))
      SD_db%Num_Elems = Count(METIS_epart == MyRank+1)
      Allocate(SD_db%Elem(SD_db%Num_Elems))


      iELoc = 0
      Do iE = 1, Geom%Num_Elems
         If (METIS_epart(iE) == MyRank + 1) Then
            SD_db%IsLocal_Elem(iE) = .TRUE.
            iELoc = iELoc + 1
            SD_db%Elem(iELoc) = iE
         End If
      End Do

      iNLoc = 0
      Do iN = 1, Geom%Num_Nodes
         If (METIS_npart(iN) == MyRank + 1) Then
            SD_db%IsLocal_Node(iN) = .TRUE.
            iNLoc = iNLoc + 1
            SD_db%Node(iNLoc) = iN
         End If
      End Do

      ! Step 3: Find the Ghosts Nodes and Elems
      SD_db%Num_GhostElems = 0
      Allocate(SD_db%GhostElem(0))

      Do iELoc = 1,SD_db%Num_Elems
         iE = SD_db%Elem(iELoc)
         ! Here we don't assume that the number of nodes per element
         ! is the same everywhere. Hopefully, when this issue is fixed 
         ! with METIS, the following will then be more general.
         Do iNLoc = 1, Elem_db(iE)%NB_DoF
            iN = Elem_db(iE)%ID_DoF(iNLoc)
            If ( .NOT. SD_db%IsLocal_Node(iN)) Then
               ! iN belongs to a local element, but is not a local node
               !! it is therefore a ghost node
               SD_db%IsLocal_GhostNode(iN) = .TRUE.
            End If
         End Do
      End Do
      SD_db%Num_GhostNodes = Count(SD_db%IsLocal_GhostNode)
      Allocate (SD_db%GhostNode(SD_db%Num_GhostNodes))

      iNLoc = 0
      Do iN = 1, Geom%Num_Nodes
         If (SD_db%IsLocal_GhostNode(iN)) Then
            iNLoc = iNLoc + 1
            SD_db%GhostNode(iNLoc) = iN
         End If
      End Do

      ! Step 4:IS, AO's
      ! This is silly!
      ! I should instead manually compute MyIdxMin, and MyIdxMax
      ! OTH, this way I am sure that this will not break, even if PETSc
      ! change its way of distributing values accross proocessors

      Call VecCreateMPI(PETSC_COMM_WORLD, SD_db%Num_Nodes, Geom%Num_Nodes,    &
           & Vect_Dist, iErr)
      Call VecGetOwnerShipRange(Vect_Dist, MyIdxMin, MyIdxMax, iErr)
      Call VecDestroy(Vect_Dist, iErr)
      Allocate(Tmp_Idx(SD_db%Num_Nodes))
      Tmp_Idx = (/ ( i, i = MyIdxMin, MyIdxMax-1) /)
      Call ISCreateGeneral(PETSC_COMM_WORLD, SD_db%Num_Nodes, Tmp_Idx,        &
           & SD_db%EXO_IS, iErr)
      Tmp_Idx = SD_db%Node - 1
      Call ISCreateGeneral(PETSC_COMM_WORLD, SD_db%Num_Nodes, Tmp_Idx,        &
           & SD_db%PETSc_IS, iErr)

      Call AOCreateBasicIS(SD_db%PETSc_IS, SD_db%EXO_IS, SD_db%EXO_AO, iErr)
      DeAllocate(Tmp_Idx)

      Allocate(Tmp_Idx(Geom%Num_Nodes))
      Allocate(Tmp_Idx2(Geom%Num_Nodes))
      
      Tmp_Idx2 = (/ (i, i = 0, Geom%Num_Nodes-1) /)
      Tmp_Idx(1:SD_db%Num_Nodes) = SD_db%Node - 1
      Tmp_Idx(SD_db%Num_Nodes + 1: SD_db%Num_Nodes + SD_db%Num_GhostNodes) =  &
           & SD_db%GhostNode - 1
      inG = SD_db%Num_Nodes + SD_db%Num_GhostNodes


      Do i = 1, Geom%Num_Nodes
         If (.NOT. (SD_db%IsLocal_Node(i) .OR. SD_db%IsLocal_GhostNode(i) ))  &
              & Then
            inG = inG + 1
            Tmp_Idx(inG) = i - 1
         End If
      End Do
      Call AOCreateBasic(PETSC_COMM_SELF, Geom%Num_Nodes, Tmp_Idx, Tmp_Idx2,  &
           & SD_db%Loc_AO, iErr)

     DeAllocate (Tmp_Idx)
     DeAllocate (Tmp_Idx2)

      ! Step 5: Scatter and PETSc ghost points
      Allocate(Tmp_Idx(SD_db%Num_GhostNodes))
      Tmp_Idx = SD_db%GhostNode-1
      Call AOApplicationToPetsc(SD_db%EXO_AO, SD_db%Num_GhostNodes, Tmp_Idx,  &
           & iErr)
      Call VecCreateGhost(PETSC_Comm_World, SD_db%Num_Nodes, Geom%Num_Nodes,  &
           & SD_db%Num_GhostNodes, Tmp_Idx, Vect_Dist, iErr)
      DeAllocate (Tmp_Idx)

      Call VecGhostGetLocalForm(Vect_Dist, Vect_Local, iErr)
      If(MyRank == 0) Then
         Call VecCreateMPI(PETSC_COMM_WORLD, Geom%Num_Nodes, Geom%Num_Nodes,  &
              & Vect_Master, iErr)
      Else
         Call VecCreateMPI(PETSC_COMM_WORLD, 0, Geom%Num_Nodes, Vect_Master,  &
              & iErr)
      End If

      Call VecScatterCreate(Vect_Dist, SD_db%EXO_IS, Vect_Master,             &
           & SD_db%PETSc_IS, SD_db%ToMaster, iErr)

    End Subroutine Init_SD_NoOvlp_2DScal
    
    Subroutine Init_SD_NoOvlp_2DElast(Geom, Elem_db, SD_db, Vect_Dist,        &
         & Vect_Local, Vect_Master)
      Type (EXO_Geom_Info), Intent(IN)                 :: Geom
      Type (Element2D_Elast), Dimension(:), Pointer    :: Elem_db
      Type (SD_Info), Intent(OUT)                      :: SD_db
      Vec, Intent(INOUT)                               :: Vect_Dist
      Vec, Intent(INOUT)                               :: Vect_Local
      Vec, Intent(INOUT)                               :: Vect_Master


      Integer, Dimension(:), Pointer                   :: METIS_epart
      Integer, Dimension(:), Pointer                   :: METIS_npart
      Integer, Dimension(:), Pointer                   :: METIS_elmnts
      Integer                                          :: METIS_edgecut
      
      Integer                                          :: MyRank
      Integer                                          :: NumProcs
      Integer                                          :: Nb_DoF
      Integer                                          :: iE, iELoc
      Integer                                          :: iN, iNLoc
      Integer                                          :: iErr

      Integer, Dimension(:), Pointer                   :: Tmp_Idx, Tmp_Idx2
      Integer                                          :: MyIdxMin, MyIdxMax
      Integer                                          :: i, iG, inG
      Integer                                          :: nDim = 2

      
      ! No overlap => no ghost elements
      Allocate(SD_db%IsLocal_Node(Geom%Num_Nodes*nDim))
      Allocate(SD_db%IsLocal_GhostNode(Geom%Num_Nodes*nDim))
      Allocate(SD_db%IsLocal_Elem(Geom%Num_Elems))
      Allocate(SD_db%IsLocal_GhostElem(Geom%Num_Elems))
      SD_db%IsLocal_Node      = .FALSE.
      SD_db%IsLocal_GhostNode = .FALSE.
      SD_db%IsLocal_Elem      = .FALSE.
      SD_db%IsLocal_GhostElem = .FALSE.

      Nb_DoF = 3
      Call MPI_COMM_RANK(MPI_COMM_WORLD, MyRank, iErr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, NumProcs, iErr)
 
      ! Step 1: Call METIS and generate the subdomains
      Allocate(METIS_elmnts(Geom%Num_Elems*Nb_DoF))
      Allocate(METIS_epart(Geom%Num_Elems))
      Allocate(METIS_npart(Geom%Num_Nodes))

      ! This assumes that all elements in the mesh are triangular.
      ! Unfortunately, METIS doesn't seem to be able to deal with 
      ! mixed elements meshes.
      
      Allocate (METIS_elmnts(Geom%Num_Elems*Nb_DoF))
      Do iE = 1, Geom%Num_Elems
!!$         METIS_elmnts((iE-1)*Nb_DoF +1:iE*Nb_DoF) =                           &
!!$              & (Elem_db(iE)%ID_DoF(1:Elem_db(iE)%Nb_DoF:nDim)-1)/nDim+1
         METIS_elmnts((iE-1)*Nb_DoF +1:iE*Nb_DoF) =                           &
              & (Elem_db(iE)%ID_DoF(1:Elem_db(iE)%Nb_DoF/nDim)+1)/nDim
      End Do

      If (NumProcs == 1) Then
         METIS_epart = 1
         METIS_npart = 1
      Else
         Call METIS_PartMeshNodal(Geom%Num_Elems, Geom%Num_Nodes,             &
              & METIS_elmnts, METIS_etype_tri, METIS_numflag, NumProcs,       &
              & METIS_edgecut, METIS_epart, METIS_npart)
      End If
!      Call METIS_PartMeshNodal(Geom%Num_Elems, Geom%Num_Nodes, METIS_elmnts,  &
!           & METIS_etype_tri, METIS_numflag, NumProcs, METIS_edgecut,         &
!           & METIS_epart, METIS_npart)

      ! Step 2: Allocate the fields in SD_Info
      SD_db%Num_Nodes = Count(METIS_npart == MyRank+1)*nDim
      Allocate(SD_db%Node(SD_db%Num_Nodes))
      SD_db%Num_Elems = Count(METIS_epart == MyRank+1)
      Allocate(SD_db%Elem(SD_db%Num_Elems))


      iELoc = 0
      Do iE = 1, Geom%Num_Elems
         If (METIS_epart(iE) == MyRank + 1) Then
            SD_db%IsLocal_Elem(iE) = .TRUE.
            iELoc = iELoc + 1
            SD_db%Elem(iELoc) = iE
         End If
      End Do

      iNLoc = 0
      Do iN = 1, Geom%Num_Nodes
         If (METIS_npart(iN) == MyRank + 1) Then
            SD_db%IsLocal_Node(nDim*iN-1) = .TRUE.
            SD_db%IsLocal_Node(nDim*iN)   = .TRUE.
            SD_db%Node(iNLoc+1) = nDim*iN-1
            SD_db%Node(iNLoc+2) = nDim*iN
            iNLoc = iNLoc + nDim
!!$            SD_db%IsLocal_Node(2*iN-1) = .TRUE.
!!$            SD_db%IsLocal_Node(2*iN)   = .TRUE.
!!$            SD_db%Node(iNLoc+1) = 2*iN-1
!!$            SD_db%Node(iNLoc+2)   = 2*iN
!!$            iNLoc = iNLoc + 2
         End If
      End Do

      ! Step 3: Find the Ghosts Nodes and Elems
      SD_db%Num_GhostElems = 0
      Allocate(SD_db%GhostElem(0))

      Do iELoc = 1,SD_db%Num_Elems
         iE = SD_db%Elem(iELoc)
         ! Here we don't assume that the number of nodes per element
         ! is the same everywhere. Hopefully, when this issue is fixed 
         ! with METIS, the following will then be more general.
         Do iNLoc = 1, Elem_db(iE)%NB_DoF
            iN = Elem_db(iE)%ID_DoF(iNLoc)
            If ( .NOT. SD_db%IsLocal_Node(iN)) Then
               ! iN belongs to a local element, but is not a local node
               !! it is therefore a ghost node
               SD_db%IsLocal_GhostNode(iN) = .TRUE.
            End If
         End Do
      End Do
      SD_db%Num_GhostNodes = Count(SD_db%IsLocal_GhostNode)
      Allocate (SD_db%GhostNode(SD_db%Num_GhostNodes))

      iNLoc = 0
      Do iN = 1, Geom%Num_Nodes*nDim
         If (SD_db%IsLocal_GhostNode(iN)) Then
            iNLoc = iNLoc + 1
            SD_db%GhostNode(iNLoc) = iN
         End If
      End Do

      ! Step 4:IS, AO's
      ! This is silly!
      ! I should instead manually compute MyIdxMin, and MyIdxMax
      ! OTH, this way I am sure that this will not break, even if PETSc
      ! change its way of distributing values accross proocessors
      Call VecCreateMPI(PETSC_COMM_WORLD, SD_db%Num_Nodes,                    &
           & Geom%Num_Nodes*nDim, Vect_Dist, iErr)
      Call VecGetOwnerShipRange(Vect_Dist, MyIdxMin, MyIdxMax, iErr)
      Call VecDestroy(Vect_Dist, iErr)
      Allocate(Tmp_Idx(SD_db%Num_Nodes))
      Tmp_Idx = (/ ( i, i = MyIdxMin, MyIdxMax-1) /)
      Call ISCreateGeneral(PETSC_COMM_WORLD, SD_db%Num_Nodes, Tmp_Idx,        &
           & SD_db%EXO_IS, iErr)
      Tmp_Idx = SD_db%Node - 1
      Call ISCreateGeneral(PETSC_COMM_WORLD, SD_db%Num_Nodes, Tmp_Idx,        &
           & SD_db%PETSc_IS, iErr)

      Call AOCreateBasicIS(SD_db%PETSc_IS, SD_db%EXO_IS, SD_db%EXO_AO, iErr)
      DeAllocate(Tmp_Idx)

      Allocate(Tmp_Idx(Geom%Num_Nodes*nDim))
      Allocate(Tmp_Idx2(Geom%Num_Nodes*nDim))
      
      Tmp_Idx2 = (/ (i, i = 0, Geom%Num_Nodes*nDim-1) /)
      Tmp_Idx(1:SD_db%Num_Nodes) = SD_db%Node - 1
      Tmp_Idx(SD_db%Num_Nodes + 1: SD_db%Num_Nodes + SD_db%Num_GhostNodes) =  &
           & SD_db%GhostNode - 1
      inG = SD_db%Num_Nodes + SD_db%Num_GhostNodes

      Do i = 1, Geom%Num_Nodes * nDim 
         If (.NOT. (SD_db%IsLocal_Node(i) .OR. SD_db%IsLocal_GhostNode(i) ))  &
              & Then
            inG = inG + 1
            Tmp_Idx(inG) = i - 1
         End If
      End Do
      Call AOCreateBasic(PETSC_COMM_SELF, Geom%Num_Nodes*nDim, Tmp_Idx,       &
           & Tmp_Idx2, SD_db%Loc_AO, iErr)
     DeAllocate (Tmp_Idx)
     DeAllocate (Tmp_Idx2)


      ! Step 5: Scatter and PETSc ghost points
      Allocate(Tmp_Idx(SD_db%Num_GhostNodes))
      Tmp_Idx = SD_db%GhostNode-1
      Call AOApplicationToPetsc(SD_db%EXO_AO, SD_db%Num_GhostNodes, Tmp_Idx,  &
           & iErr)
      Call VecCreateGhost(PETSC_Comm_World, SD_db%Num_Nodes,                  &
           & Geom%Num_Nodes * nDim, SD_db%Num_GhostNodes, Tmp_Idx, Vect_Dist, &
           & iErr)
      DeAllocate (Tmp_Idx)

      Call VecGhostGetLocalForm(Vect_Dist, Vect_Local, iErr)
      If(MyRank == 0) Then
         Call VecCreateMPI(PETSC_COMM_WORLD, Geom%Num_Nodes * nDim,           &
              & Geom%Num_Nodes * nDim, Vect_Master, iErr)
      Else
         Call VecCreateMPI(PETSC_COMM_WORLD, 0, Geom%Num_Nodes * nDim,        &
              & Vect_Master, iErr)
      End If

      Call VecScatterCreate(Vect_Dist, SD_db%EXO_IS, Vect_Master,             &
           & SD_db%PETSc_IS, SD_db%ToMaster, iErr)

    End Subroutine Init_SD_NoOvlp_2DElast

    Subroutine Init_SD_NoOvlp_2D(Geom, Elem_db, SD_db, Vect_Dist, Vect_Local, &
         & Vect_Master)
      Type (EXO_Geom_Info), Intent(IN)                 :: Geom
      Type (Element2D), Dimension(:), Pointer          :: Elem_db
      Type (SD_Info), Intent(OUT)                      :: SD_db
      Vec, Intent(INOUT)                               :: Vect_Dist
      Vec, Intent(INOUT)                               :: Vect_Local
      Vec, Intent(INOUT)                               :: Vect_Master


      Integer, Dimension(:), Pointer                   :: METIS_epart
      Integer, Dimension(:), Pointer                   :: METIS_npart
      Integer, Dimension(:), Pointer                   :: METIS_elmnts
      Integer                                          :: METIS_edgecut
      
      Integer                                          :: MyRank
      Integer                                          :: NumProcs
      Integer                                          :: Nb_DoF
      Integer                                          :: iE, iELoc
      Integer                                          :: iN, iNLoc
      Integer                                          :: iErr

      Integer, Dimension(:), Pointer                   :: Tmp_Idx, Tmp_Idx2
      Integer                                          :: MyIdxMin, MyIdxMax
      Integer                                          :: i, iG, inG
      Integer                                          :: nDim = 2

      
      ! No overlap => no ghost elements
      Allocate(SD_db%IsLocal_Node(Geom%Num_Nodes*nDim))
      Allocate(SD_db%IsLocal_GhostNode(Geom%Num_Nodes*nDim))
      Allocate(SD_db%IsLocal_Elem(Geom%Num_Elems))
      Allocate(SD_db%IsLocal_GhostElem(Geom%Num_Elems))
      SD_db%IsLocal_Node      = .FALSE.
      SD_db%IsLocal_GhostNode = .FALSE.
      SD_db%IsLocal_Elem      = .FALSE.
      SD_db%IsLocal_GhostElem = .FALSE.

      Nb_DoF = 3
      Call MPI_COMM_RANK(MPI_COMM_WORLD, MyRank, iErr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, NumProcs, iErr)
 
      ! Step 1: Call METIS and generate the subdomains
      Allocate(METIS_elmnts(Geom%Num_Elems*Nb_DoF))
      Allocate(METIS_epart(Geom%Num_Elems))
      Allocate(METIS_npart(Geom%Num_Nodes))

      ! This assumes that all elements in the mesh are triangular.
      ! Unfortunately, METIS doesn't seem to be able to deal with 
      ! mixed elements meshes.
      
      Allocate (METIS_elmnts(Geom%Num_Elems*Nb_DoF))
      Do iE = 1, Geom%Num_Elems
!!$         METIS_elmnts((iE-1)*Nb_DoF +1:iE*Nb_DoF) =                           &
!!$              & (Elem_db(iE)%ID_DoF(1:Elem_db(iE)%Nb_DoF:nDim)-1)/nDim+1
         METIS_elmnts((iE-1)*Nb_DoF +1:iE*Nb_DoF) =                           &
              & (Elem_db(iE)%ID_DoF(1:Elem_db(iE)%Nb_DoF/nDim)+1)/nDim
      End Do

      If (NumProcs == 1) Then
         METIS_epart = 1
         METIS_npart = 1
      Else
         Call METIS_PartMeshNodal(Geom%Num_Elems, Geom%Num_Nodes,             &
              & METIS_elmnts, METIS_etype_tri, METIS_numflag, NumProcs,       &
              & METIS_edgecut, METIS_epart, METIS_npart)
      End If
!      Call METIS_PartMeshNodal(Geom%Num_Elems, Geom%Num_Nodes, METIS_elmnts,  &
!           & METIS_etype_tri, METIS_numflag, NumProcs, METIS_edgecut,         &
!           & METIS_epart, METIS_npart)

      ! Step 2: Allocate the fields in SD_Info
      SD_db%Num_Nodes = Count(METIS_npart == MyRank+1)*nDim
      Allocate(SD_db%Node(SD_db%Num_Nodes))
      SD_db%Num_Elems = Count(METIS_epart == MyRank+1)
      Allocate(SD_db%Elem(SD_db%Num_Elems))


      iELoc = 0
      Do iE = 1, Geom%Num_Elems
         If (METIS_epart(iE) == MyRank + 1) Then
            SD_db%IsLocal_Elem(iE) = .TRUE.
            iELoc = iELoc + 1
            SD_db%Elem(iELoc) = iE
         End If
      End Do

      iNLoc = 0
      Do iN = 1, Geom%Num_Nodes
         If (METIS_npart(iN) == MyRank + 1) Then
            SD_db%IsLocal_Node(nDim*iN-1) = .TRUE.
            SD_db%IsLocal_Node(nDim*iN)   = .TRUE.
            SD_db%Node(iNLoc+1) = nDim*iN-1
            SD_db%Node(iNLoc+2) = nDim*iN
            iNLoc = iNLoc + nDim
!!$            SD_db%IsLocal_Node(2*iN-1) = .TRUE.
!!$            SD_db%IsLocal_Node(2*iN)   = .TRUE.
!!$            SD_db%Node(iNLoc+1) = 2*iN-1
!!$            SD_db%Node(iNLoc+2)   = 2*iN
!!$            iNLoc = iNLoc + 2
         End If
      End Do

      ! Step 3: Find the Ghosts Nodes and Elems
      SD_db%Num_GhostElems = 0
      Allocate(SD_db%GhostElem(0))

      Do iELoc = 1,SD_db%Num_Elems
         iE = SD_db%Elem(iELoc)
         ! Here we don't assume that the number of nodes per element
         ! is the same everywhere. Hopefully, when this issue is fixed 
         ! with METIS, the following will then be more general.
         Do iNLoc = 1, Elem_db(iE)%NB_DoF
            iN = Elem_db(iE)%ID_DoF(iNLoc)
            If ( .NOT. SD_db%IsLocal_Node(iN)) Then
               ! iN belongs to a local element, but is not a local node
               !! it is therefore a ghost node
               SD_db%IsLocal_GhostNode(iN) = .TRUE.
            End If
         End Do
      End Do
      SD_db%Num_GhostNodes = Count(SD_db%IsLocal_GhostNode)
      Allocate (SD_db%GhostNode(SD_db%Num_GhostNodes))

      iNLoc = 0
      Do iN = 1, Geom%Num_Nodes*nDim
         If (SD_db%IsLocal_GhostNode(iN)) Then
            iNLoc = iNLoc + 1
            SD_db%GhostNode(iNLoc) = iN
         End If
      End Do

      ! Step 4:IS, AO's
      ! This is silly!
      ! I should instead manually compute MyIdxMin, and MyIdxMax
      ! OTH, this way I am sure that this will not break, even if PETSc
      ! change its way of distributing values accross proocessors
      Call VecCreateMPI(PETSC_COMM_WORLD, SD_db%Num_Nodes,                    &
           & Geom%Num_Nodes*nDim, Vect_Dist, iErr)
      Call VecGetOwnerShipRange(Vect_Dist, MyIdxMin, MyIdxMax, iErr)
      Call VecDestroy(Vect_Dist, iErr)
      Allocate(Tmp_Idx(SD_db%Num_Nodes))
      Tmp_Idx = (/ ( i, i = MyIdxMin, MyIdxMax-1) /)
      Call ISCreateGeneral(PETSC_COMM_WORLD, SD_db%Num_Nodes, Tmp_Idx,        &
           & SD_db%EXO_IS, iErr)
      Tmp_Idx = SD_db%Node - 1
      Call ISCreateGeneral(PETSC_COMM_WORLD, SD_db%Num_Nodes, Tmp_Idx,        &
           & SD_db%PETSc_IS, iErr)

      Call AOCreateBasicIS(SD_db%PETSc_IS, SD_db%EXO_IS, SD_db%EXO_AO, iErr)
      DeAllocate(Tmp_Idx)

      Allocate(Tmp_Idx(Geom%Num_Nodes*nDim))
      Allocate(Tmp_Idx2(Geom%Num_Nodes*nDim))
      
      Tmp_Idx2 = (/ (i, i = 0, Geom%Num_Nodes*nDim-1) /)
      Tmp_Idx(1:SD_db%Num_Nodes) = SD_db%Node - 1
      Tmp_Idx(SD_db%Num_Nodes + 1: SD_db%Num_Nodes + SD_db%Num_GhostNodes) =  &
           & SD_db%GhostNode - 1
      inG = SD_db%Num_Nodes + SD_db%Num_GhostNodes

      Do i = 1, Geom%Num_Nodes * nDim 
         If (.NOT. (SD_db%IsLocal_Node(i) .OR. SD_db%IsLocal_GhostNode(i) ))  &
              & Then
            inG = inG + 1
            Tmp_Idx(inG) = i - 1
         End If
      End Do
      Call AOCreateBasic(PETSC_COMM_SELF, Geom%Num_Nodes*nDim, Tmp_Idx,       &
           & Tmp_Idx2, SD_db%Loc_AO, iErr)
     DeAllocate (Tmp_Idx)
     DeAllocate (Tmp_Idx2)


      ! Step 5: Scatter and PETSc ghost points
      Allocate(Tmp_Idx(SD_db%Num_GhostNodes))
      Tmp_Idx = SD_db%GhostNode-1
      Call AOApplicationToPetsc(SD_db%EXO_AO, SD_db%Num_GhostNodes, Tmp_Idx,  &
           & iErr)
      Call VecCreateGhost(PETSC_Comm_World, SD_db%Num_Nodes,                  &
           & Geom%Num_Nodes * nDim, SD_db%Num_GhostNodes, Tmp_Idx, Vect_Dist, &
           & iErr)
      DeAllocate (Tmp_Idx)

      Call VecGhostGetLocalForm(Vect_Dist, Vect_Local, iErr)
      If(MyRank == 0) Then
         Call VecCreateMPI(PETSC_COMM_WORLD, Geom%Num_Nodes * nDim,           &
              & Geom%Num_Nodes * nDim, Vect_Master, iErr)
      Else
         Call VecCreateMPI(PETSC_COMM_WORLD, 0, Geom%Num_Nodes * nDim,        &
              & Vect_Master, iErr)
      End If

      Call VecScatterCreate(Vect_Dist, SD_db%EXO_IS, Vect_Master,             &
           & SD_db%PETSc_IS, SD_db%ToMaster, iErr)

    End Subroutine Init_SD_NoOvlp_2D

    Subroutine Init_SD_NoOvlp_3DScal(Geom, Elem_db, SD_db, Vect_Dist,         &
         & Vect_Local, Vect_Master)
      Type (EXO_Geom_Info), Intent(IN)                 :: Geom
      Type (Element3D_Scal), Dimension(:), Pointer     :: Elem_db
      Type (SD_Info), Intent(OUT)                      :: SD_db
      Vec, Intent(INOUT)                               :: Vect_Dist
      Vec, Intent(INOUT)                               :: Vect_Local
      Vec, Intent(INOUT)                               :: Vect_Master


      Integer, Dimension(:), Pointer                   :: METIS_epart
      Integer, Dimension(:), Pointer                   :: METIS_npart
      Integer, Dimension(:), Pointer                   :: METIS_elmnts
      Integer                                          :: METIS_edgecut
      
      Integer                                          :: MyRank
      Integer                                          :: NumProcs
      Integer                                          :: Nb_DoF
      Integer                                          :: iE, iELoc
      Integer                                          :: iN, iNLoc
      Integer                                          :: iErr

      Integer, Dimension(:), Pointer                   :: Tmp_Idx, Tmp_Idx2
      Integer                                          :: MyIdxMin, MyIdxMax
      Integer                                          :: i, iG, inG

      
      ! No overlap => no ghost elements
      Allocate(SD_db%IsLocal_Node(Geom%Num_Nodes))
      Allocate(SD_db%IsLocal_GhostNode(Geom%Num_Nodes))
      Allocate(SD_db%IsLocal_Elem(Geom%Num_Elems))
      Allocate(SD_db%IsLocal_GhostElem(Geom%Num_Elems))
      SD_db%IsLocal_Node      = .FALSE.
      SD_db%IsLocal_GhostNode = .FALSE.
      SD_db%IsLocal_Elem      = .FALSE.
      SD_db%IsLocal_GhostElem = .FALSE.

      Nb_DoF = 4

      Call MPI_COMM_RANK(MPI_COMM_WORLD, MyRank, iErr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, NumProcs, iErr)
 
      ! Step 1: Call METIS and generate the subdomains
      Allocate(METIS_elmnts(Geom%Num_Elems*Nb_DoF))
      Allocate(METIS_epart(Geom%Num_Elems))
      Allocate(METIS_npart(Geom%Num_Nodes))

      ! This assumes that all elements in the mesh are tets.
      ! Unfortunately, METIS doesn't seem to be able to deal with 
      ! mixed elements meshes.
      
      Allocate (METIS_elmnts(Geom%Num_Elems*Nb_DoF))
      Do iE = 1, Geom%Num_Elems
         METIS_elmnts((iE-1)*Nb_DoF +1:iE*Nb_DoF) = Elem_db(iE)%ID_DoF
      End Do

      If (NumProcs == 1) Then
         METIS_epart = 1
         METIS_npart = 1
      Else
         Call METIS_PartMeshNodal(Geom%Num_Elems, Geom%Num_Nodes,             &
              & METIS_elmnts, METIS_etype_tet, METIS_numflag, NumProcs,       &
              & METIS_edgecut, METIS_epart, METIS_npart)
      End If
!      Call METIS_PartMeshNodal(Geom%Num_Elems, Geom%Num_Nodes, METIS_elmnts,  &
!           & METIS_etype_tet, METIS_numflag, NumProcs, METIS_edgecut,         &
!           & METIS_epart, METIS_npart)

      ! Step 2: Allocate the fields in SD_Info
      SD_db%Num_Nodes = Count(METIS_npart == MyRank+1)
      Allocate(SD_db%Node(SD_db%Num_Nodes))
      SD_db%Num_Elems = Count(METIS_epart == MyRank+1)
      Allocate(SD_db%Elem(SD_db%Num_Elems))


      iELoc = 0
      Do iE = 1, Geom%Num_Elems
         If (METIS_epart(iE) == MyRank + 1) Then
            SD_db%IsLocal_Elem(iE) = .TRUE.
            iELoc = iELoc + 1
            SD_db%Elem(iELoc) = iE
         End If
      End Do

      iNLoc = 0
      Do iN = 1, Geom%Num_Nodes
         If (METIS_npart(iN) == MyRank + 1) Then
            SD_db%IsLocal_Node(iN) = .TRUE.
            iNLoc = iNLoc + 1
            SD_db%Node(iNLoc) = iN
         End If
      End Do

      ! Step 3: Find the Ghosts Nodes and Elems
      SD_db%Num_GhostElems = 0
      Allocate(SD_db%GhostElem(0))

      Do iELoc = 1,SD_db%Num_Elems
         iE = SD_db%Elem(iELoc)
         ! Here we don't assume that the number of nodes per element
         ! is the same everywhere. Hopefully, when this issue is fixed 
         ! with METIS, the following will then be more general.
         Do iNLoc = 1, Elem_db(iE)%NB_DoF
            iN = Elem_db(iE)%ID_DoF(iNLoc)
            If ( .NOT. SD_db%IsLocal_Node(iN)) Then
               ! iN belongs to a local element, but is not a local node
               !! it is therefore a ghost node
               SD_db%IsLocal_GhostNode(iN) = .TRUE.
            End If
         End Do
      End Do
      SD_db%Num_GhostNodes = Count(SD_db%IsLocal_GhostNode)
      Allocate (SD_db%GhostNode(SD_db%Num_GhostNodes))

      iNLoc = 0
      Do iN = 1, Geom%Num_Nodes
         If (SD_db%IsLocal_GhostNode(iN)) Then
            iNLoc = iNLoc + 1
            SD_db%GhostNode(iNLoc) = iN
         End If
      End Do

      ! Step 4:IS, AO's
      ! This is silly!
      ! I should instead manually compute MyIdxMin, and MyIdxMax
      ! OTH, this way I am sure that this will not break, even if PETSc
      ! change its way of distributing values accross proocessors

      Call VecCreateMPI(PETSC_COMM_WORLD, SD_db%Num_Nodes, Geom%Num_Nodes,    &
           & Vect_Dist, iErr)
      Call VecGetOwnerShipRange(Vect_Dist, MyIdxMin, MyIdxMax, iErr)
      Call VecDestroy(Vect_Dist, iErr)
      Allocate(Tmp_Idx(SD_db%Num_Nodes))
      Tmp_Idx = (/ ( i, i = MyIdxMin, MyIdxMax-1) /)
      Call ISCreateGeneral(PETSC_COMM_WORLD, SD_db%Num_Nodes, Tmp_Idx,        &
           & SD_db%EXO_IS, iErr)
      Tmp_Idx = SD_db%Node - 1
      Call ISCreateGeneral(PETSC_COMM_WORLD, SD_db%Num_Nodes, Tmp_Idx,        &
           & SD_db%PETSc_IS, iErr)

      Call AOCreateBasicIS(SD_db%PETSc_IS, SD_db%EXO_IS, SD_db%EXO_AO, iErr)
      DeAllocate(Tmp_Idx)

      Allocate(Tmp_Idx(Geom%Num_Nodes))
      Allocate(Tmp_Idx2(Geom%Num_Nodes))
      
      Tmp_Idx2 = (/ (i, i = 0, Geom%Num_Nodes-1) /)
      Tmp_Idx(1:SD_db%Num_Nodes) = SD_db%Node - 1
      Tmp_Idx(SD_db%Num_Nodes + 1: SD_db%Num_Nodes + SD_db%Num_GhostNodes) =  &
           & SD_db%GhostNode - 1
      inG = SD_db%Num_Nodes + SD_db%Num_GhostNodes


      Do i = 1, Geom%Num_Nodes
         If (.NOT. (SD_db%IsLocal_Node(i) .OR. SD_db%IsLocal_GhostNode(i) ))  &
              & Then
            inG = inG + 1
            Tmp_Idx(inG) = i - 1
         End If
      End Do
      Call AOCreateBasic(PETSC_COMM_SELF, Geom%Num_Nodes, Tmp_Idx, Tmp_Idx2,  &
           & SD_db%Loc_AO, iErr)

     DeAllocate (Tmp_Idx)
     DeAllocate (Tmp_Idx2)

      ! Step 5: Scatter and PETSc ghost points
      Allocate(Tmp_Idx(SD_db%Num_GhostNodes))
      Tmp_Idx = SD_db%GhostNode-1
      Call AOApplicationToPetsc(SD_db%EXO_AO, SD_db%Num_GhostNodes, Tmp_Idx,  &
           & iErr)
      Call VecCreateGhost(PETSC_Comm_World, SD_db%Num_Nodes, Geom%Num_Nodes,  &
           & SD_db%Num_GhostNodes, Tmp_Idx, Vect_Dist, iErr)
      DeAllocate (Tmp_Idx)

      Call VecGhostGetLocalForm(Vect_Dist, Vect_Local, iErr)
      If(MyRank == 0) Then
         Call VecCreateMPI(PETSC_COMM_WORLD, Geom%Num_Nodes, Geom%Num_Nodes,  &
              & Vect_Master, iErr)
      Else
         Call VecCreateMPI(PETSC_COMM_WORLD, 0, Geom%Num_Nodes, Vect_Master,  &
              & iErr)
      End If

      Call VecScatterCreate(Vect_Dist, SD_db%EXO_IS, Vect_Master,             &
           & SD_db%PETSc_IS, SD_db%ToMaster, iErr)

    End Subroutine Init_SD_NoOvlp_3DScal
    
    Subroutine Init_SD_NoOvlp_3DElast(Geom, Elem_db, SD_db, Vect_Dist,        &
         & Vect_Local, Vect_Master)
      Type (EXO_Geom_Info), Intent(IN)                 :: Geom
      Type (Element3D_Elast), Dimension(:), Pointer    :: Elem_db
      Type (SD_Info), Intent(OUT)                      :: SD_db
      Vec, Intent(INOUT)                               :: Vect_Dist
      Vec, Intent(INOUT)                               :: Vect_Local
      Vec, Intent(INOUT)                               :: Vect_Master


      Integer, Dimension(:), Pointer                   :: METIS_epart
      Integer, Dimension(:), Pointer                   :: METIS_npart
      Integer, Dimension(:), Pointer                   :: METIS_elmnts
      Integer                                          :: METIS_edgecut
      
      Integer                                          :: MyRank
      Integer                                          :: NumProcs
      Integer                                          :: Nb_DoF
      Integer                                          :: iE, iELoc
      Integer                                          :: iN, iNLoc
      Integer                                          :: iErr

      Integer, Dimension(:), Pointer                   :: Tmp_Idx, Tmp_Idx2
      Integer                                          :: MyIdxMin, MyIdxMax
      Integer                                          :: i, iG, inG
      Integer                                          :: nDim = 3

      
      ! No overlap => no ghost elements
      Allocate(SD_db%IsLocal_Node(Geom%Num_Nodes*nDim))
      Allocate(SD_db%IsLocal_GhostNode(Geom%Num_Nodes*nDim))
      Allocate(SD_db%IsLocal_Elem(Geom%Num_Elems))
      Allocate(SD_db%IsLocal_GhostElem(Geom%Num_Elems))
      SD_db%IsLocal_Node      = .FALSE.
      SD_db%IsLocal_GhostNode = .FALSE.
      SD_db%IsLocal_Elem      = .FALSE.
      SD_db%IsLocal_GhostElem = .FALSE.

      Nb_DoF = 4
      Call MPI_COMM_RANK(MPI_COMM_WORLD, MyRank, iErr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, NumProcs, iErr)
 
      ! Step 1: Call METIS and generate the subdomains
      Allocate(METIS_elmnts(Geom%Num_Elems*Nb_DoF))
      Allocate(METIS_epart(Geom%Num_Elems))
      Allocate(METIS_npart(Geom%Num_Nodes))

      ! This assumes that all elements in the mesh are tets.
      ! Unfortunately, METIS doesn't seem to be able to deal with 
      ! mixed elements meshes.
      
      Allocate (METIS_elmnts(Geom%Num_Elems*Nb_DoF))
      Do iE = 1, Geom%Num_Elems
!!$         METIS_elmnts((iE-1)*Nb_DoF +1:iE*Nb_DoF) =                           &
!!$              & (Elem_db(iE)%ID_DoF(1:Elem_db(iE)%Nb_DoF:nDim)-1)/nDim+1
         METIS_elmnts((iE-1)*Nb_DoF +1:iE*Nb_DoF) =                           &
              & (Elem_db(iE)%ID_DoF(1:Elem_db(iE)%Nb_DoF/nDim)+2)/nDim
      End Do

      If (NumProcs == 1) Then
         METIS_epart = 1
         METIS_npart = 1
      Else
         Call METIS_PartMeshNodal(Geom%Num_Elems, Geom%Num_Nodes,             &
              & METIS_elmnts, METIS_etype_tet, METIS_numflag, NumProcs,       &
              & METIS_edgecut, METIS_epart, METIS_npart)
      End If
!      If (NumProcs > 1) Then
!      Call METIS_PartMeshNodal(Geom%Num_Elems, Geom%Num_Nodes, METIS_elmnts,  &
!           & METIS_etype_tet, METIS_numflag, NumProcs, METIS_edgecut,         &
!           & METIS_epart, METIS_npart)!
!      Else
!         METIS_epart = 1
!         METIS_npart = 1
!      End If

      ! Step 2: Allocate the fields in SD_Info
      SD_db%Num_Nodes = Count(METIS_npart == MyRank+1)*nDim
      Allocate(SD_db%Node(SD_db%Num_Nodes))
      SD_db%Num_Elems = Count(METIS_epart == MyRank+1)
      Allocate(SD_db%Elem(SD_db%Num_Elems))


      iELoc = 0
      Do iE = 1, Geom%Num_Elems
         If (METIS_epart(iE) == MyRank + 1) Then
            SD_db%IsLocal_Elem(iE) = .TRUE.
            iELoc = iELoc + 1
            SD_db%Elem(iELoc) = iE
         End If
      End Do

      iNLoc = 0
      Do iN = 1, Geom%Num_Nodes
         If (METIS_npart(iN) == MyRank + 1) Then
            SD_db%IsLocal_Node(nDim*iN-2) = .TRUE.
            SD_db%IsLocal_Node(nDim*iN-1) = .TRUE.
            SD_db%IsLocal_Node(nDim*iN)   = .TRUE.
            SD_db%Node(iNLoc+1) = nDim*iN-2
            SD_db%Node(iNLoc+2) = nDim*iN-1
            SD_db%Node(iNLoc+3) = nDim*iN
            iNLoc = iNLoc + nDim
         End If
      End Do

      ! Step 3: Find the Ghosts Nodes and Elems
      SD_db%Num_GhostElems = 0
      Allocate(SD_db%GhostElem(0))

      Do iELoc = 1,SD_db%Num_Elems
         iE = SD_db%Elem(iELoc)
         ! Here we don't assume that the number of nodes per element
         ! is the same everywhere. Hopefully, when this issue is fixed 
         ! with METIS, the following will then be more general.
         Do iNLoc = 1, Elem_db(iE)%NB_DoF
            iN = Elem_db(iE)%ID_DoF(iNLoc)
            If ( .NOT. SD_db%IsLocal_Node(iN)) Then
               ! iN belongs to a local element, but is not a local node
               !! it is therefore a ghost node
               SD_db%IsLocal_GhostNode(iN) = .TRUE.
            End If
         End Do
      End Do
      SD_db%Num_GhostNodes = Count(SD_db%IsLocal_GhostNode)
      Allocate (SD_db%GhostNode(SD_db%Num_GhostNodes))

      iNLoc = 0
      Do iN = 1, Geom%Num_Nodes*nDim
         If (SD_db%IsLocal_GhostNode(iN)) Then
            iNLoc = iNLoc + 1
            SD_db%GhostNode(iNLoc) = iN
         End If
      End Do

      ! Step 4:IS, AO's
      ! This is silly!
      ! I should instead manually compute MyIdxMin, and MyIdxMax
      ! OTH, this way I am sure that this will not break, even if PETSc
      ! change its way of distributing values accross proocessors
      Call VecCreateMPI(PETSC_COMM_WORLD, SD_db%Num_Nodes,                    &
           & Geom%Num_Nodes*nDim, Vect_Dist, iErr)
      Call VecGetOwnerShipRange(Vect_Dist, MyIdxMin, MyIdxMax, iErr)
      Call VecDestroy(Vect_Dist, iErr)
      Allocate(Tmp_Idx(SD_db%Num_Nodes))
      Tmp_Idx = (/ ( i, i = MyIdxMin, MyIdxMax-1) /)
      Call ISCreateGeneral(PETSC_COMM_WORLD, SD_db%Num_Nodes, Tmp_Idx,        &
           & SD_db%EXO_IS, iErr)

      Tmp_Idx = SD_db%Node - 1
      Call ISCreateGeneral(PETSC_COMM_WORLD, SD_db%Num_Nodes, Tmp_Idx,        &
           & SD_db%PETSc_IS, iErr)

      Call AOCreateBasicIS(SD_db%PETSc_IS, SD_db%EXO_IS, SD_db%EXO_AO, iErr)
      DeAllocate(Tmp_Idx)

      Allocate(Tmp_Idx(Geom%Num_Nodes*nDim))
      Allocate(Tmp_Idx2(Geom%Num_Nodes*nDim))
      
      Tmp_Idx2 = (/ (i, i = 0, Geom%Num_Nodes*nDim-1) /)
      Tmp_Idx(1:SD_db%Num_Nodes) = SD_db%Node - 1
      Tmp_Idx(SD_db%Num_Nodes + 1: SD_db%Num_Nodes + SD_db%Num_GhostNodes) =  &
           & SD_db%GhostNode - 1
      inG = SD_db%Num_Nodes + SD_db%Num_GhostNodes

      Do i = 1, Geom%Num_Nodes * nDim 
         If (.NOT. (SD_db%IsLocal_Node(i) .OR. SD_db%IsLocal_GhostNode(i) ))  &
              & Then
            inG = inG + 1
            Tmp_Idx(inG) = i - 1
         End If
      End Do
      Call AOCreateBasic(PETSC_COMM_SELF, Geom%Num_Nodes*nDim, Tmp_Idx,       &
           & Tmp_Idx2, SD_db%Loc_AO, iErr)
     DeAllocate (Tmp_Idx)
     DeAllocate (Tmp_Idx2)


      ! Step 5: Scatter and PETSc ghost points
      Allocate(Tmp_Idx(SD_db%Num_GhostNodes))
      Tmp_Idx = SD_db%GhostNode-1
      Call AOApplicationToPetsc(SD_db%EXO_AO, SD_db%Num_GhostNodes, Tmp_Idx,  &
           & iErr)
      Call VecCreateGhost(PETSC_Comm_World, SD_db%Num_Nodes,                  &
           & Geom%Num_Nodes * nDim, SD_db%Num_GhostNodes, Tmp_Idx, Vect_Dist, &
           & iErr)
      DeAllocate (Tmp_Idx)

      Call VecGhostGetLocalForm(Vect_Dist, Vect_Local, iErr)
      If(MyRank == 0) Then
         Call VecCreateMPI(PETSC_COMM_WORLD, Geom%Num_Nodes * nDim,           &
              & Geom%Num_Nodes * nDim, Vect_Master, iErr)
      Else
         Call VecCreateMPI(PETSC_COMM_WORLD, 0, Geom%Num_Nodes * nDim,        &
              & Vect_Master, iErr)
      End If

      Call VecScatterCreate(Vect_Dist, SD_db%EXO_IS, Vect_Master,             &
           & SD_db%PETSc_IS, SD_db%ToMaster, iErr)

    End Subroutine Init_SD_NoOvlp_3DElast

    Subroutine Init_SD_NoOvlp_3D(Geom, Elem_db, SD_db, Vect_Dist, Vect_Local, &
         & Vect_Master)
      Type (EXO_Geom_Info), Intent(IN)                 :: Geom
      Type (Element3D), Dimension(:), Pointer          :: Elem_db
      Type (SD_Info), Intent(OUT)                      :: SD_db
      Vec, Intent(INOUT)                               :: Vect_Dist
      Vec, Intent(INOUT)                               :: Vect_Local
      Vec, Intent(INOUT)                               :: Vect_Master


      Integer, Dimension(:), Pointer                   :: METIS_epart
      Integer, Dimension(:), Pointer                   :: METIS_npart
      Integer, Dimension(:), Pointer                   :: METIS_elmnts
      Integer                                          :: METIS_edgecut
      
      Integer                                          :: MyRank
      Integer                                          :: NumProcs
      Integer                                          :: Nb_DoF
      Integer                                          :: iE, iELoc
      Integer                                          :: iN, iNLoc
      Integer                                          :: iErr

      Integer, Dimension(:), Pointer                   :: Tmp_Idx, Tmp_Idx2
      Integer                                          :: MyIdxMin, MyIdxMax
      Integer                                          :: i, iG, inG
      Integer                                          :: nDim = 3

      
      ! No overlap => no ghost elements
      Allocate(SD_db%IsLocal_Node(Geom%Num_Nodes*nDim))
      Allocate(SD_db%IsLocal_GhostNode(Geom%Num_Nodes*nDim))
      Allocate(SD_db%IsLocal_Elem(Geom%Num_Elems))
      Allocate(SD_db%IsLocal_GhostElem(Geom%Num_Elems))
      SD_db%IsLocal_Node      = .FALSE.
      SD_db%IsLocal_GhostNode = .FALSE.
      SD_db%IsLocal_Elem      = .FALSE.
      SD_db%IsLocal_GhostElem = .FALSE.

      Nb_DoF = 4
      Call MPI_COMM_RANK(MPI_COMM_WORLD, MyRank, iErr)
      Call MPI_COMM_SIZE(MPI_COMM_WORLD, NumProcs, iErr)
 
      ! Step 1: Call METIS and generate the subdomains
      Allocate(METIS_elmnts(Geom%Num_Elems*Nb_DoF))
      Allocate(METIS_epart(Geom%Num_Elems))
      Allocate(METIS_npart(Geom%Num_Nodes))

      ! This assumes that all elements in the mesh are tets.
      ! Unfortunately, METIS doesn't seem to be able to deal with 
      ! mixed elements meshes.
      
      Allocate (METIS_elmnts(Geom%Num_Elems*Nb_DoF))
      Do iE = 1, Geom%Num_Elems
!!$         METIS_elmnts((iE-1)*Nb_DoF +1:iE*Nb_DoF) =                           &
!!$              & (Elem_db(iE)%ID_DoF(1:Elem_db(iE)%Nb_DoF:nDim)-1)/nDim+1
         METIS_elmnts((iE-1)*Nb_DoF +1:iE*Nb_DoF) =                           &
              & (Elem_db(iE)%ID_DoF(1:Elem_db(iE)%Nb_DoF/nDim)+2)/nDim
      End Do

      If (NumProcs == 1) Then
         METIS_epart = 1
         METIS_npart = 1
      Else
         Call METIS_PartMeshNodal(Geom%Num_Elems, Geom%Num_Nodes,             &
              & METIS_elmnts, METIS_etype_tet, METIS_numflag, NumProcs,       &
              & METIS_edgecut, METIS_epart, METIS_npart)
      End If
!      Call METIS_PartMeshNodal(Geom%Num_Elems, Geom%Num_Nodes, METIS_elmnts,  &
!           & METIS_etype_tet, METIS_numflag, NumProcs, METIS_edgecut,         &
!           & METIS_epart, METIS_npart)

      ! Step 2: Allocate the fields in SD_Info
      SD_db%Num_Nodes = Count(METIS_npart == MyRank+1)*nDim
      Allocate(SD_db%Node(SD_db%Num_Nodes))
      SD_db%Num_Elems = Count(METIS_epart == MyRank+1)
      Allocate(SD_db%Elem(SD_db%Num_Elems))


      iELoc = 0
      Do iE = 1, Geom%Num_Elems
         If (METIS_epart(iE) == MyRank + 1) Then
            SD_db%IsLocal_Elem(iE) = .TRUE.
            iELoc = iELoc + 1
            SD_db%Elem(iELoc) = iE
         End If
      End Do

      iNLoc = 0
      Do iN = 1, Geom%Num_Nodes
         If (METIS_npart(iN) == MyRank + 1) Then
            SD_db%IsLocal_Node(nDim*iN-2) = .TRUE.
            SD_db%IsLocal_Node(nDim*iN-1) = .TRUE.
            SD_db%IsLocal_Node(nDim*iN)   = .TRUE.
            SD_db%Node(iNLoc+1) = nDim*iN-2
            SD_db%Node(iNLoc+2) = nDim*iN-1
            SD_db%Node(iNLoc+3) = nDim*iN
            iNLoc = iNLoc + nDim
         End If
      End Do

      ! Step 3: Find the Ghosts Nodes and Elems
      SD_db%Num_GhostElems = 0
      Allocate(SD_db%GhostElem(0))

      Do iELoc = 1,SD_db%Num_Elems
         iE = SD_db%Elem(iELoc)
         ! Here we don't assume that the number of nodes per element
         ! is the same everywhere. Hopefully, when this issue is fixed 
         ! with METIS, the following will then be more general.
         Do iNLoc = 1, Elem_db(iE)%NB_DoF
            iN = Elem_db(iE)%ID_DoF(iNLoc)
            If ( .NOT. SD_db%IsLocal_Node(iN)) Then
               ! iN belongs to a local element, but is not a local node
               !! it is therefore a ghost node
               SD_db%IsLocal_GhostNode(iN) = .TRUE.
            End If
         End Do
      End Do
      SD_db%Num_GhostNodes = Count(SD_db%IsLocal_GhostNode)
      Allocate (SD_db%GhostNode(SD_db%Num_GhostNodes))

      iNLoc = 0
      Do iN = 1, Geom%Num_Nodes*nDim
         If (SD_db%IsLocal_GhostNode(iN)) Then
            iNLoc = iNLoc + 1
            SD_db%GhostNode(iNLoc) = iN
         End If
      End Do

      ! Step 4:IS, AO's
      ! This is silly!
      ! I should instead manually compute MyIdxMin, and MyIdxMax
      ! OTH, this way I am sure that this will not break, even if PETSc
      ! change its way of distributing values accross proocessors
      Call VecCreateMPI(PETSC_COMM_WORLD, SD_db%Num_Nodes,                    &
           & Geom%Num_Nodes*nDim, Vect_Dist, iErr)
      Call VecGetOwnerShipRange(Vect_Dist, MyIdxMin, MyIdxMax, iErr)
      Call VecDestroy(Vect_Dist, iErr)
      Allocate(Tmp_Idx(SD_db%Num_Nodes))
      Tmp_Idx = (/ ( i, i = MyIdxMin, MyIdxMax-1) /)
      Call ISCreateGeneral(PETSC_COMM_WORLD, SD_db%Num_Nodes, Tmp_Idx,        &
           & SD_db%EXO_IS, iErr)
      Tmp_Idx = SD_db%Node - 1
      Call ISCreateGeneral(PETSC_COMM_WORLD, SD_db%Num_Nodes, Tmp_Idx,        &
           & SD_db%PETSc_IS, iErr)

      Call AOCreateBasicIS(SD_db%PETSc_IS, SD_db%EXO_IS, SD_db%EXO_AO, iErr)
      DeAllocate(Tmp_Idx)

      Allocate(Tmp_Idx(Geom%Num_Nodes*nDim))
      Allocate(Tmp_Idx2(Geom%Num_Nodes*nDim))
      
      Tmp_Idx2 = (/ (i, i = 0, Geom%Num_Nodes*nDim-1) /)
      Tmp_Idx(1:SD_db%Num_Nodes) = SD_db%Node - 1
      Tmp_Idx(SD_db%Num_Nodes + 1: SD_db%Num_Nodes + SD_db%Num_GhostNodes) =  &
           & SD_db%GhostNode - 1
      inG = SD_db%Num_Nodes + SD_db%Num_GhostNodes

      Do i = 1, Geom%Num_Nodes * nDim 
         If (.NOT. (SD_db%IsLocal_Node(i) .OR. SD_db%IsLocal_GhostNode(i) ))  &
              & Then
            inG = inG + 1
            Tmp_Idx(inG) = i - 1
         End If
      End Do
      Call AOCreateBasic(PETSC_COMM_SELF, Geom%Num_Nodes*nDim, Tmp_Idx,       &
           & Tmp_Idx2, SD_db%Loc_AO, iErr)
     DeAllocate (Tmp_Idx)
     DeAllocate (Tmp_Idx2)


      ! Step 5: Scatter and PETSc ghost points
      Allocate(Tmp_Idx(SD_db%Num_GhostNodes))
      Tmp_Idx = SD_db%GhostNode-1
      Call AOApplicationToPetsc(SD_db%EXO_AO, SD_db%Num_GhostNodes, Tmp_Idx,  &
           & iErr)
      Call VecCreateGhost(PETSC_Comm_World, SD_db%Num_Nodes,                  &
           & Geom%Num_Nodes * nDim, SD_db%Num_GhostNodes, Tmp_Idx, Vect_Dist, &
           & iErr)
      DeAllocate (Tmp_Idx)

      Call VecGhostGetLocalForm(Vect_Dist, Vect_Local, iErr)
      If(MyRank == 0) Then
         Call VecCreateMPI(PETSC_COMM_WORLD, Geom%Num_Nodes * nDim,           &
              & Geom%Num_Nodes * nDim, Vect_Master, iErr)
      Else
         Call VecCreateMPI(PETSC_COMM_WORLD, 0, Geom%Num_Nodes * nDim,        &
              & Vect_Master, iErr)
      End If

      Call VecScatterCreate(Vect_Dist, SD_db%EXO_IS, Vect_Master,             &
           & SD_db%PETSc_IS, SD_db%ToMaster, iErr)

    End Subroutine Init_SD_NoOvlp_3D

End Module m_MEF_SD
