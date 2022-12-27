Module m_MEF90_EXO
#include "petsc/finclude/petsc.h"
   Use m_MEF90_DMPlex
   Use m_MEF90_Parameters
   Use m_MEF90_Ctx
   Use m_MEF90_Utils
   Use m_MEF90_Elements
   Use petsc
   IMPLICIT NONE
#include "../mef90version.h"

   Private 
   PetscInt,Public                                 :: exo_ver

   
   Public :: MEF90CtxOpenEXO
   Public :: MEF90CtxCloseEXO
   Public :: MEF90EXOFormat
   Public :: MEF90EXODMView
   Public :: MEF90EXOVecView
   Public :: MEF90EXOVecLoad

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90CtxOpenEXO"
!!!
!!!  
!!!  MEF90CtxOpenEXO:
!!!  
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca
!!!      2022 Blaise Bourdin  bourdin@mcmaster.ca
!!!
   Subroutine MEF90CtxOpenEXO(MEF90Ctx,Viewer,mode,ierr)
      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx
      Type(tPetscViewer), Intent(INOUT)               :: Viewer
      PetscEnum,Intent(IN)                            :: mode
      PetscErrorCode,Intent(INOUT)                    :: ierr

      Integer                                         :: opts

#ifdef PETSC_USE_DEBUG
      opts = EXVRBS + EXDEBG
#else
      opts = 0
#endif
      Call exopts(opts,ierr)
      PetscCall(PetscViewerExodusIIOpen(MEF90Ctx%Comm,MEF90Ctx%resultFile,mode,Viewer,ierr))
   End Subroutine MEF90CtxOpenEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxCloseEXO"
!!!
!!!  
!!!  MEF90CtxCloseEXO:
!!!  
!!!  (c) 2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!      2022      Blaise Bourdin  bourdin@mcmaster.ca
!!!
   Subroutine MEF90CtxCloseEXO(Viewer,ierr)
      Type(tPetscViewer),Intent(INOUT)                :: Viewer
      PetscErrorCode,Intent(INOUT)                    :: ierr

      PetscCall(PetscViewerDestroy(Viewer,ierr))
   End Subroutine MEF90CtxCloseEXO
      
#undef __FUNCT__
#define __FUNCT__ "MEF90EXOFormat"
!!!
!!!  
!!!  MEF90EXOFormat:
!!!  
!!!  (c) 2012-2022 Blaise Bourdin bourdin@lsu.edu
!!!           2022 Alexis Marboeuf marboeua@mcmaster.ca    
!!!
Subroutine MEF90EXOFormat(Viewer,nameG,nameC,nameV,nameS,time,ierr)
   Type(tPetscViewer),Intent(IN)                         :: Viewer
   Character(len=*),Dimension(:),Intent(IN)              :: nameG,nameC,nameV,nameS
   PetscReal,Dimension(:),Pointer                        :: time
   PetscErrorCode,Intent(INOUT)                          :: ierr
   
   PetscInt                                              :: numCS,numSS,numG,numC,numV,numS
   Integer                                               :: exoid
   PetscInt                                              :: step
   Character(len=MXSTLN)                                 :: sJunk
   PetscReal                                             :: rJunk
   Logical,Dimension(:,:),Pointer                        :: truthtable


   If (.NOT. associated(time)) Then
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_POINTER,"Time value must be allocated prior to calling MEF90EXOFormat")
      STOP
   End If
   PetscCall(PetscViewerExodusIIGetId(Viewer,exoid,ierr))

   If (exoid > 0) Then
      Call exinq(exoid,EX_INQ_SIDE_SETS,numSS,rJunk,sjunk,ierr)
      Call exinq(exoid, EX_INQ_ELEM_BLK,numCS,rJunk,sjunk,ierr)
      !!! Write variable names
      numG = size(nameG)
      If (numG > 0) Then
         Call expvp(exoid,"g",numG,ierr)
         Call expvan(exoid,"g",numG,nameG,ierr)
      End If
      numC = size(nameC)
       If (numC > 0) Then
          Call expvp(exoid,"e",numC,ierr)
          Call expvan(exoid,"e",numC,nameC,ierr)
       End If
       numV = size(nameV)
       If (numV > 0) Then
          Call expvp(exoid,"n",numV,ierr)
          Call expvan(exoid,"n",numV,nameV,ierr)
       End If
       numS = size(nameS)
       If (numS > 0) Then
          PetscCall(expvp(exoid,"s",numS,ierr))
          PetscCall(expvan(exoid,"s",numS,nameS,ierr))
       End If

       !!! Write truth tables
       If (numS > 0) Then
          Allocate(truthtable(numSS,numS))
          truthtable = .true.
          PetscCall(expsstt(exoid, numSS, numS, truthtable, ierr))
          DeAllocate(truthtable)
       End If

       If (numC > 0) Then
          Allocate(truthtable(numCS,numC))
          truthtable = .true.
          Call expvtt(exoid, numCS, numC, truthtable, ierr)
          DeAllocate(truthtable)
       End If

      Do step = 1,size(time)
          Call exptim(exoid,step,time(step),ierr)
       End Do
   End If
End Subroutine MEF90EXOFormat

#undef __FUNCT__
#define __FUNCT__ "MEF90EXODMView"
!!!
!!!  
!!!  MEF90EXODMView:
!!!
!!!  
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca    
!!!
   Subroutine MEF90EXODMView(dm,Viewer,order,ierr)
      Type(tPetscViewer),Intent(IN)                      :: Viewer
      Type(tDM),Intent(IN)                               :: dm
      PetscInt,Intent(IN)                                :: order
      PetscErrorCode,Intent(INOUT)                       :: ierr
      
      Character(len=PETSC_MAX_PATH_LEN)                  :: IOBuffer

      if ((order > 2) .or. (order < 1)) then
          write(IOBuffer,'("Unsupported polynomial order ", I2, " not in [1,2]")') order
          SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,IOBuffer)
      end if
      PetscCall(PetscViewerExodusIISetOrder(Viewer,order,ierr))
      PetscCall(DMView(dm,Viewer,ierr))
   End Subroutine MEF90EXODMView

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecView"
!!!
!!!  
!!!  MEF90EXOVecView:
!!!
!!!  
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca    
!!!
   Subroutine MEF90EXOVecView(v,sf,invSF,Viewer,step,ierr)
      Type(tVec),Intent(IN)                              :: v
      Type(tPetscSF),Intent(IN)                          :: sf,invSF
      Type(tPetscViewer),Intent(IN)                      :: Viewer
      PetscInt,Intent(IN)                                :: step
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Integer                                            :: exoid
      PetscInt                                           :: offsetN = -1,offsetZ = -1, offsetS = -1,bs
      Type(tVec)                                         :: iov
      Character(len=PETSC_MAX_PATH_LEN)                  :: vecname,IOBuffer

      PetscCall(PetscViewerExodusIIGetId(Viewer,exoid,ierr))
      PetscCall(PetscObjectGetName(v,vecname,ierr))

      PetscCall(VecGetBlockSize(v,bs,ierr))
      PetscCall(MEF90VecCreateIO(iov,bs,sf,ierr))
      PetscCall(PetscObjectSetName(iov,vecname,ierr))
      PetscCall(MEF90VecCopySF(v,iov,sf,ierr))

      PetscCall(MEF90EXOGetVarIndex_Private(exoid,"n",vecname,offsetN,ierr))
      PetscCall(MEF90EXOGetVarIndex_Private(exoid,"e",vecname,offsetZ,ierr))
      PetscCall(MEF90EXOGetVarIndex_Private(exoid,"s",vecname,offsetS,ierr))
      If (offsetN > 0) Then
         PetscCall(MEF90EXOVecViewNodal_Private(iov,exoid,step,offsetN,ierr))
      Else If (offsetZ > 0) Then
         PetscCall(MEF90EXOVecViewZonal_Private(iov,exoid,step,offsetZ,ierr))
      Else If (offsetS > 0) Then
         PetscCall(MEF90EXOVecViewSide_Private(iov,exoid,step,offsetS,ierr))
      Else
         write(IOBuffer,'("Could not find nodal or zonal variable ", A, " in exodus file. ")') trim(vecname)
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_UNEXPECTED,IOBuffer)
      End If
      PetscCall(VecDestroy(iov,ierr))
   End Subroutine MEF90EXOVecView

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecLoad"
!!!
!!!  
!!!  MEF90EXOVecLoad:
!!!
!!!  
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca    
!!!
   Subroutine MEF90EXOVecLoad(v,sf,invSF,Viewer,step,ierr)
      Type(tVec),Intent(INOUT)                           :: v
      Type(tPetscSF),Intent(IN)                          :: sf,invSF
      Type(tPetscViewer),Intent(IN)                      :: Viewer
      PetscInt,Intent(IN)                                :: step
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Integer                                            :: exoid 
      PetscInt                                           :: offsetN,offsetZ,offsetS,bs
      Type(tVec)                                         :: iov
      Character(len=PETSC_MAX_PATH_LEN)                  :: vecname,IOBuffer

      offsetN = -1
      offsetZ = -1
      offsetS = -1
      PetscCall(PetscViewerExodusIIGetId(Viewer,exoid,ierr))
      PetscCall(PetscObjectGetName(v,vecname,ierr))

      PetscCall(VecGetBlockSize(v,bs,ierr))
      PetscCall(MEF90VecCreateIO(iov,bs,sf,ierr))
      PetscCall(PetscObjectSetName(iov,vecname,ierr))

      PetscCall(MEF90EXOGetVarIndex_Private(exoid,"n",vecname,offsetN,ierr))
      PetscCall(MEF90EXOGetVarIndex_Private(exoid,"e",vecname,offsetZ,ierr))
      PetscCall(MEF90EXOGetVarIndex_Private(exoid,"s",vecname,offsetS,ierr))
      If (offsetN > 0) Then
         PetscCall(MEF90EXOVecLoadNodal_Private(iov,exoid,step,offsetN,ierr))
      Else If (offsetZ > 0) Then
         PetscCall(MEF90EXOVecLoadZonal_Private(iov,exoid,step,offsetZ,ierr))
      Else If (offsetS > 0) Then
         PetscCall(MEF90EXOVecLoadSide_Private(iov,exoid,step,offsetS,ierr))
      Else
         write(IOBuffer,'("Could not find nodal or zonal variable ", A, " in exodus file. ")') trim(vecname)
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_UNEXPECTED,IOBuffer)
      End If
      PetscCall(MEF90VecCopySF(iov,v,invSF,ierr))
      PetscCall(VecDestroy(iov,ierr))
   End Subroutine MEF90EXOVecLoad

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOGetVarIndex_Private"
   Subroutine MEF90EXOGetVarIndex_Private(exoid,obj_type,name,varIndex,ierr)
      Integer,Intent(IN)               :: exoid
      Character(len=MXSTLN),Intent(IN) :: name
      Character(len=1),Intent(IN)      :: obj_type
      PetscInt,Intent(OUT)             :: varIndex
      PetscErrorCode,Intent(INOUT)     :: ierr
   
      Integer                          :: i,j,num_suffix,num_vars
      Character(len=MXSTLN)            :: var_name,ext_name,suffix(5)
   
      suffix = ["   ","_X ","_XX","_1 ","_11"]
      num_suffix = size(suffix)

      varIndex = -1
      Call exgvp(exoid,obj_type,num_vars,ierr)
      Do i=1,num_vars
         Call exgvnm(exoid,obj_type,i,var_name,ierr)
         Do j=1,num_suffix
            ext_name = trim(name)//trim(suffix(j))
            If (ext_name == var_name) Then
               varIndex = i
            End If
         End Do
      End Do
   End Subroutine MEF90EXOGetVarIndex_Private
   
#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecViewNodal_Private"
   Subroutine MEF90EXOVecViewNodal_Private(v,exoid,step,offset,ierr)
      Integer,Intent(IN)               :: exoid
      PetscInt,Intent(IN)              :: step,offset
      Type(tVec),Intent(IN)            :: v
      PetscErrorCode,Intent(INOUT)     :: ierr
   
      PetscInt                         :: xs,xe,bs,c
      PetscScalar,Dimension(:),Pointer :: varray
      Type(tVec)                       :: vComp
      Type(tIS)                        :: compIS
   
      PetscCall(VecGetOwnershipRange(v,xs,xe,ierr))
      PetscCall(VecGetBlockSize(v,bs,ierr))

      If (bs == 1) Then
         PetscCall(VecGetArrayReadF90(v,varray,ierr))
         Call expnvs(exoid,step,offset,xs+1,xe-xs,varray,ierr)
         PetscCall(VecRestoreArrayReadF90(v,varray,ierr))
      Else 
         PetscCall(ISCreateStride(PETSC_COMM_WORLD,(xe-xs)/bs,xs,bs,compIS,ierr))
         Do c = 0,bs-1
            PetscCall(ISStrideSetStride(compIS,(xe-xs)/bs,xs+c,bs,ierr))
            PetscCall(VecGetSubVector(v,compIS,vComp,ierr))
            PetscCall(VecGetArrayReadF90(vComp,varray,ierr))
            Call expnvs(exoid,step,offset+c,xs/bs+1,(xe-xs)/bs,varray,ierr)
            PetscCall(VecRestoreArrayReadF90(vComp,varray,ierr))
            PetscCall(VecRestoreSubVector(v,compIS,vComp,ierr))
         End Do ! c
         PetscCall(ISDestroy(compIS,ierr))
      End If ! bs
   End Subroutine MEF90EXOVecViewNodal_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecLoadNodal_Private"
   Subroutine MEF90EXOVecLoadNodal_Private(v,exoid,step,offset,ierr)
      Integer,Intent(IN)               :: exoid
      PetscInt,Intent(IN)              :: step, offset
      Type(tVec),Intent(INOUT)         :: v
      PetscErrorCode,Intent(INOUT)     :: ierr
   
      PetscInt                         :: xs,xe,bs,c
      PetscScalar,Dimension(:),Pointer :: varray
      Type(tVec)                       :: vComp
      Type(tIS)                        :: compIS
   
      PetscCall(VecGetOwnershipRange(v,xs,xe,ierr))
      PetscCall(VecGetBlockSize(v,bs,ierr))
      If (bs == 1) Then
         PetscCall(VecGetArrayReadF90(v,varray,ierr));
         Call exgnnv(exoid,step,offset,xs+1,xe-xs,varray,ierr)
         PetscCall(VecRestoreArrayReadF90(v,varray,ierr))
      Else 
         PetscCall(ISCreateStride(PETSC_COMM_WORLD,(xe-xs)/bs,xs,bs,compIS,ierr))
         Do c = 0,bs-1
            PetscCall(ISStrideSetStride(compIS,(xe-xs)/bs,xs+c,bs,ierr))
            PetscCall(VecGetSubVector(v,compIS,vComp,ierr))
            PetscCall(VecGetArrayReadF90(vComp,varray,ierr))
            Call exgnnv(exoid,step,offset+c,xs/bs+1,(xe-xs)/bs,varray,ierr)
            PetscCall(VecRestoreArrayReadF90(vComp,varray,ierr))
            PetscCall(VecISCopy(v,compIS,SCATTER_FORWARD,vComp,ierr))
            PetscCall(VecRestoreSubVector(v,compIS,vComp,ierr))
         End Do
         PetscCall(ISDestroy(compIS,ierr))
      End If
   End Subroutine MEF90EXOVecLoadNodal_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecViewZonal_Private"
   Subroutine MEF90EXOVecViewZonal_Private(v,exoid,step,offset,ierr)
      Integer,Intent(IN)               :: exoid
      PetscInt,Intent(IN)              :: step, offset
      Type(tVec),Intent(IN)            :: v
      PetscErrorCode,Intent(INOUT)     :: ierr
   
      PetscInt                         :: xs,xe,bs,c,numCS,set,csLocalSize,csxs
      PetscScalar,Dimension(:),Pointer :: varray
      PetscInt,Dimension(:),Pointer    :: csID,csSize
      Type(tVec)                       :: vComp
      Type(tIS)                        :: compIS
      Character(len=MXSTLN)            :: elemType
      PetscMPIInt                      :: rank
   
      csxs = 0_Ki
      PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      numCS = exinqi(exoid,EX_INQ_ELEM_BLK)
      Allocate(csID(numCS))
      Allocate(csSize(numCS))
      Call exgebi(exoid,csID,ierr)
      Do set=1,numCS
         Call exgelb(exoid,csID(set),elemType,csSize(set),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)
      End Do
      PetscCall(VecGetOwnershipRange(v,xs,xe,ierr))
      PetscCall(VecGetBlockSize(v,bs,ierr))
      If (bs > 1) Then
         PetscCall(ISCreateStride(PETSC_COMM_WORLD,(xe-xs)/bs,xs,bs,compIS,ierr))
      End If
      Do set=1,numCS
         !  range of indices for set setID[set]: csxs:csxs + csSize[set]-1
         !  local slice of zonal values:         xs/bs,xm/bs-1
         !  intersection:                        max(xs/bs,csxs),min(xm/bs-1,csxs + csSize[set]-1)
         csLocalSize = max(0, min(xe/bs, csxs+csSize(set)) - max(xs/bs, csxs))
         If (bs == 1) Then
            PetscCall(VecGetArrayReadF90(v,varray,ierr))
            Call expevs(exoid,step,offset,csID(set),max(xs-csxs,0)+1,csLocalSize,varray,ierr)
            PetscCall(VecRestoreArrayReadF90(v,varray,ierr))
         Else
            Do c = 0,bs-1
               PetscCall(ISStrideSetStride(compIS,(xe-xs)/bs,xs+c,bs,ierr))
               PetscCall(VecGetSubVector(v,compIS,vComp,ierr))
               PetscCall(VecGetArrayReadF90(vComp,varray,ierr))
               Call expevs(exoid,step,offset+c,csID(set),max(xs/bs-csxs,0)+1,csLocalSize,varray(max(0,csxs-xs/bs)+1:max(0,csxs-xs/bs)+csLocalSize),ierr)
               PetscCall(VecRestoreArrayReadF90(vComp,varray,ierr))
               PetscCall(VecRestoreSubVector(v,compIS,vComp,ierr))
            End Do
         End If
         csxs = csxs + csSize(set)
      End Do
      If (bs > 1) Then 
         PetscCall(ISDestroy(compIS,ierr))
      End If
      DeAllocate(csID)
      DeAllocate(csSize)
   End Subroutine MEF90EXOVecViewZonal_Private
   
#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecLoadZonal_Private"
   Subroutine MEF90EXOVecLoadZonal_Private(v,exoid,step,offset,ierr)
      Integer,Intent(IN)               :: exoid
      PetscInt,Intent(IN)              :: step, offset
      Type(tVec),Intent(INOUT)         :: v
      PetscErrorCode,Intent(INOUT)     :: ierr
   
      PetscInt                         :: xs,xe,bs,c,numCS,set,csLocalSize,csxs
      PetscScalar,Dimension(:),Pointer :: varray
      PetscInt,Dimension(:),Pointer    :: csID,csSize
      Type(tVec)                       :: vComp
      Type(tIS)                        :: compIS
      Character(len=MXSTLN)            :: elemType
      PetscMPIInt                      :: rank
   
      csxs = 0_Ki
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      numCS = exinqi(exoid,EX_INQ_ELEM_BLK)
      Allocate(csID(numCS))
      Allocate(csSize(numCS))
      Call exgebi(exoid,csID,ierr)
      Do set=1,numCS
         Call exgelb(exoid,csID(set),elemType,csSize(set),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)
      End Do
      PetscCall(VecGetOwnershipRange(v,xs,xe,ierr))
      PetscCall(VecGetBlockSize(v,bs,ierr))
      If (bs > 1) Then
         PetscCall(ISCreateStride(PETSC_COMM_WORLD,(xe-xs)/bs,xs,bs,compIS,ierr))
      End If
      Do set=1,numCS
         !  range of indices for set setID[set]: csxs:csxs + csSize[set]-1
         !  local slice of zonal values:         xs/bs,xm/bs-1
         !  intersection:                        max(xs/bs,csxs),min(xm/bs-1,csxs + csSize[set]-1)
         csLocalSize = max(0, min(xe/bs, csxs+csSize(set)) - max(xs/bs, csxs))
         If (bs == 1) Then
            PetscCall(VecGetArrayF90(v,varray,ierr))
            Call exgnev(exoid,step,offset,csID(set),csSize(set),max(xs-csxs,0)+1,csLocalSize,varray,ierr)
            PetscCall(VecRestoreArrayF90(v,varray,ierr))
         Else
            Do c = 0,bs-1
               PetscCall(ISStrideSetStride(compIS,(xe-xs)/bs,xs+c,bs,ierr))
               PetscCall(VecGetSubVector(v,compIS,vComp,ierr))
               PetscCall(VecGetArrayF90(vComp,varray,ierr))
               Call exgnev(exoid,step,offset+c,csID(set),0_Ki,max(xs/bs-csxs,0)+1,csLocalSize,varray(max(0,csxs-xs/bs)+1:max(0,csxs-xs/bs)+csLocalSize),ierr)
               ! the 5th argument of exgnev is unused
               PetscCall(VecRestoreArrayF90(vComp,varray,ierr))
               PetscCall(VecISCopy(v,compIS,SCATTER_FORWARD,vComp,ierr))
               PetscCall(VecRestoreSubVector(v,compIS,vComp,ierr))
            End Do
         End If
         csxs = csxs + csSize(set)
      End Do
      If (bs > 1) Then 
         PetscCall(ISDestroy(compIS,ierr))
      End If
      DeAllocate(csID)
      DeAllocate(csSize)
   End Subroutine MEF90EXOVecLoadZonal_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecViewSide_Private"
Subroutine MEF90EXOVecViewSide_Private(v,exoid,step,offset,ierr)
   Integer,Intent(IN)               :: exoid
   PetscInt,Intent(IN)              :: step, offset
   Type(tVec),Intent(IN)            :: v
   PetscErrorCode,Intent(INOUT)     :: ierr

   PetscInt                         :: xs,xe,bs,c,numSS,set,ssLocalSize,ssxs,sscs
   PetscScalar,Dimension(:),Pointer :: varray
   PetscInt,Dimension(:),Pointer    :: ssID,ssSize
   Type(tVec)                       :: vComp
   Type(tIS)                        :: compIS
   PetscMPIInt                      :: rank

   ssxs = 0_Ki
   sscs = 0_Ki
   PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
   numSS = exinqi(exoid,EX_INQ_SIDE_SETS)
   Allocate(ssID(numSS))
   Allocate(ssSize(numSS))
   Call exgssi(exoid,ssID,ierr)
   Do set=1,numSS
      Call exgsp(exoid,ssID(set),ssSize(set),PETSC_NULL_INTEGER,ierr)
   End Do
   PetscCall(VecGetOwnershipRange(v,xs,xe,ierr))
   PetscCall(VecGetBlockSize(v,bs,ierr))
   Do set=1,numSS
      !  range of indices for set setID[set]: csxs:csxs + csSize[set]-1
      !  local slice of zonal values:         xs/bs,xm/bs-1
      !  intersection:                        max(xs/bs,csxs),min(xm/bs-1,csxs + csSize[set]-1)
      ssLocalSize = max(0, min(xe/bs, ssxs+ssSize(set)) - max(xs/bs, ssxs))
      If (bs == 1) Then
         PetscCall(VecGetArrayReadF90(v,varray,ierr))
         Call expssv(exoid,step,offset,ssID(set),ssLocalSize,varray,ierr)
         PetscCall(VecRestoreArrayReadF90(v,varray,ierr))
      Else
         PetscCall(ISCreateStride(PETSC_COMM_WORLD,ssLocalSize,xs+sscs,bs,compIS,ierr))
         Do c = 0,bs-1
            PetscCall(ISStrideSetStride(compIS,ssLocalSize,xs+sscs+c,bs,ierr))
            PetscCall(VecGetSubVector(v,compIS,vComp,ierr))
            PetscCall(VecGetArrayReadF90(vComp,varray,ierr))
            Call exppv(exoid,step,EX_SIDE_SET,offset+c,ssID(set),max(xs/bs-ssxs,0)+1,ssLocalSize,varray,ierr)
            PetscCall(VecRestoreArrayReadF90(vComp,varray,ierr))
            PetscCall(VecRestoreSubVector(v,compIS,vComp,ierr))
         End Do
         PetscCall(ISDestroy(compIS,ierr))
      End If
      ssxs = ssxs + ssSize(set)
      sscs = sscs + bs*ssLocalSize
   End Do
   DeAllocate(ssID)
   DeAllocate(ssSize)
End Subroutine MEF90EXOVecViewSide_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecLoadSide_Private"
Subroutine MEF90EXOVecLoadSide_Private(v,exoid,step,offset,ierr)
   Integer,Intent(IN)               :: exoid
   PetscInt,Intent(IN)              :: step, offset
   Type(tVec),Intent(IN)            :: v
   PetscErrorCode,Intent(INOUT)     :: ierr

   PetscInt                         :: xs,xe,bs,c,numSS,set,ssLocalSize,ssxs,sscs
   PetscScalar,Dimension(:),Pointer :: varray
   PetscInt,Dimension(:),Pointer    :: ssID,ssSize
   Type(tVec)                       :: vComp
   Type(tIS)                        :: compIS
   PetscMPIInt                      :: rank

   ssxs = 0_Ki
   sscs = 0_Ki
   PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
   numSS = exinqi(exoid,EX_INQ_SIDE_SETS)
   Allocate(ssID(numSS))
   Allocate(ssSize(numSS))
   Call exgssi(exoid,ssID,ierr)
   Do set=1,numSS
      Call exgsp(exoid,ssID(set),ssSize(set),PETSC_NULL_INTEGER,ierr)
   End Do
   PetscCall(VecGetOwnershipRange(v,xs,xe,ierr))
   PetscCall(VecGetBlockSize(v,bs,ierr))
   Do set=1,numSS
      !  range of indices for set setID[set]: csxs:csxs + csSize[set]-1
      !  local slice of zonal values:         xs/bs,xm/bs-1
      !  intersection:                        max(xs/bs,csxs),min(xm/bs-1,csxs + csSize[set]-1)
      ssLocalSize = max(0, min(xe/bs, ssxs+ssSize(set)) - max(xs/bs, ssxs))
      If (bs == 1) Then
         PetscCall(VecGetArrayReadF90(v,varray,ierr))
         PetscCall(exgssv(exoid,step,offset,ssID(set),ssLocalSize,varray,ierr))
         PetscCall(VecRestoreArrayReadF90(v,varray,ierr))
      Else
         PetscCall(ISCreateStride(PETSC_COMM_WORLD,ssLocalSize,xs+sscs,bs,compIS,ierr))
         Do c = 0,bs-1
            PetscCall(ISStrideSetStride(compIS,ssLocalSize,xs+sscs+c,bs,ierr))
            PetscCall(VecGetSubVector(v,compIS,vComp,ierr))
            PetscCall(VecGetArrayReadF90(vComp,varray,ierr))
            Call exgpv(exoid,step,EX_SIDE_SET,offset+c,ssID(set),max(xs/bs-ssxs,0)+1,ssLocalSize,varray,ierr)
            PetscCall(VecRestoreArrayReadF90(vComp,varray,ierr))
            PetscCall(VecISCopy(v,compIS,SCATTER_FORWARD,vComp,ierr))
            PetscCall(VecRestoreSubVector(v,compIS,vComp,ierr))
         End Do
         PetscCall(ISDestroy(compIS,ierr))
      End If
      ssxs = ssxs + ssSize(set)
      sscs = sscs + bs*ssLocalSize
   End Do
   DeAllocate(ssID)
   DeAllocate(ssSize)
End Subroutine MEF90EXOVecLoadSide_Private
End Module m_MEF90_EXO
