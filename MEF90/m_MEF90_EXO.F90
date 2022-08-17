Module m_MEF90_EXO
#include "petsc/finclude/petsc.h"
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
!!!
   Subroutine MEF90CtxOpenEXO(MEF90Ctx,Viewer,ierr)
      Type(MEF90Ctx_Type),Intent(IN)                  :: MEF90Ctx
      Type(tPetscViewer), Intent(INOUT)               :: Viewer
      PetscErrorCode,Intent(INOUT)                    :: ierr

      PetscCallA(exopts(EXVRBS+EXDEBG,ierr))
      PetscCallA(PetscViewerExodusIIOpen(MEF90Ctx%Comm,MEF90Ctx%resultFile,FILE_MODE_WRITE,Viewer,ierr))

   End Subroutine MEF90CtxOpenEXO

#undef __FUNCT__
#define __FUNCT__ "MEF90CtxCloseEXO"
!!!
!!!  
!!!  MEF90CtxCloseEXO:
!!!  
!!! 
!!!      2022      Alexis Marboeuf marboeua@mcmaster.ca
!!!
   Subroutine MEF90CtxCloseEXO(Viewer,ierr)
      Type(tPetscViewer),Intent(INOUT)                :: Viewer
      PetscErrorCode,Intent(INOUT)                    :: ierr

      PetscCallA(PetscViewerDestroy(Viewer,ierr))

   End Subroutine MEF90CtxCloseEXO
      
#undef __FUNCT__
#define __FUNCT__ "MEF90EXOFormat"
!!!
!!!  
!!!  MEF90EXOFormat:
!!!  
!!!  (c) 2012-2014 Blaise Bourdin bourdin@lsu.edu
!!!           2022 Alexis Marboeuf marboeua@mcmaster.ca    
!!!
   Subroutine MEF90EXOFormat(Viewer,NameG,NameC,NameV,numstep,ierr)
      Type(tPetscViewer),Intent(IN)                         :: Viewer
      Character(len=*),Dimension(:),Pointer,Intent(IN)      :: nameG,nameC,nameV
      PetscInt,Intent(IN)                                   :: numstep
      PetscErrorCode,Intent(INOUT)                          :: ierr
      
      PetscInt                                              :: numCS
      Integer                                               :: exoid
      PetscInt                                              :: step
      Character(len=MXSTLN)                                 :: sJunk
      Logical,Dimension(:,:),Pointer                        :: truthtable

      PetscCallA(PetscViewerExodusIIGetId(Viewer,exoid,ierr))
      If (exoid > 0) Then
         !!! Write variable names
         If (size(nameG) > 0) Then
            PetscCallA(expvp(exoid,"g",size(nameG),ierr))
            PetscCallA(expvan(exoid,"g",size(nameG),nameG,ierr))
         End If
         If (size(nameC) > 0) Then
            PetscCallA(expvp(exoid,"e",size(nameC),ierr))
            PetscCallA(expvan(exoid,"e",size(nameC),nameC,ierr))
         End If
         If (size(nameV) > 0) Then
            PetscCallA(expvp(exoid,"n",size(nameV),ierr))
            PetscCallA(expvan(exoid,"n",size(nameV),nameV,ierr))
         End If

         !!! Write truth tables
         PetscCallA(exinq(exoid, EX_INQ_ELEM_BLK,numCS,PETSC_NULL_REAL,sjunk,ierr))
         If (size(nameC) > 0) Then
            Allocate(truthtable(numCS,size(nameC)))
            truthtable = .true.
            PetscCallA(expvtt(exoid, numCS, size(nameC), truthtable, ierr))
            DeAllocate(truthtable)
         End If

         Do step = 1,numstep
            PetscCallA(exptim(exoid,step,Real(step,kind=Kr),ierr))
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
   Subroutine MEF90EXODMView(dm,Viewer,ierr)
      Type(tPetscViewer),Intent(IN)                      :: Viewer
      Type(tDM),Intent(IN)                               :: dm
      PetscErrorCode,Intent(INOUT)                       :: ierr
      
      PetscInt                                           :: order = 1
      PetscBool                                          :: flg
      Character(len=PETSC_MAX_PATH_LEN)                  :: IOBuffer

      PetscCallA(PetscOptionsGetInt(PETSC_NULL_OPTIONS,PETSC_NULL_CHARACTER,"-order",order,flg,ierr))
      if ((order > 2) .or. (order < 1)) then
          write(IOBuffer,'("Unsupported polynomial order ", I2, " not in [1,2]")') order
          SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,IOBuffer)
      end if
      PetscCallA(PetscViewerExodusIISetOrder(Viewer,order,ierr))
      PetscCallA(DMView(dm,Viewer,ierr))

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
   Subroutine MEF90EXOVecView(v,Viewer,ierr)
      Type(tVec),Intent(IN)                              :: v
      Type(tPetscViewer),Intent(IN)                      :: Viewer
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(tDM)                                          :: dm
      Integer                                            :: exoid
      PetscInt                                           :: offsetN = -1, offsetZ = -1, step
      PetscReal                                          :: time
      Character(len=PETSC_MAX_PATH_LEN)                  :: vecname, IOBuffer

      PetscCallA(PetscViewerExodusIIGetId(Viewer,exoid,ierr))
      PetscCallA(PetscObjectGetName(v, vecname,ierr))
      PetscCallA(VecGetDM(v,dm,ierr))
      PetscCallA(DMGetOutputSequenceNumber(dm,step,time,ierr))

      PetscCallA(MEF90EXOGetVarIndex_Private(exoid,"n",vecname,offsetN,ierr))
      PetscCallA(MEF90EXOGetVarIndex_Private(exoid,"e",vecname,offsetZ,ierr))
      If (offsetN > 0) Then
         PetscCallA(MEF90EXOVecViewNodal_Private(v,exoid,step+1,offsetN,ierr))
      Else If (offsetZ > 0) Then
         PetscCallA(MEF90EXOVecViewZonal_Private(v,exoid,step+1,offsetZ,ierr))
      Else
         write(IOBuffer,'("Could not find nodal or zonal variable ", A5, " in exodus file. ")') vecname
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_UNEXPECTED,IOBuffer)
      End If
   End Subroutine MEF90EXOVecView

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecView"
!!!
!!!  
!!!  MEF90EXOVecLoad:
!!!
!!!  
!!!  (c) 2022 Alexis Marboeuf marboeua@mcmaster.ca    
!!!
   Subroutine MEF90EXOVecLoad(v,Viewer,ierr)
      Type(tVec),Intent(IN)                              :: v
      Type(tPetscViewer),Intent(IN)                      :: Viewer
      PetscErrorCode,Intent(INOUT)                       :: ierr

      Type(tDM)                                          :: dm
      Integer                                            :: exoid
      PetscInt                                           :: offsetN = -1, offsetZ = -1, step
      PetscReal                                          :: time
      Character(len=PETSC_MAX_PATH_LEN)                  :: vecname, IOBuffer

      PetscCallA(PetscViewerExodusIIGetId(Viewer,exoid,ierr))
      PetscCallA(PetscObjectGetName(v, vecname,ierr))
      PetscCallA(VecGetDM(v,dm,ierr))
      PetscCallA(DMGetOutputSequenceNumber(dm,step,time,ierr))

      PetscCallA(MEF90EXOGetVarIndex_Private(exoid,"n",vecname,offsetN,ierr))
      PetscCallA(MEF90EXOGetVarIndex_Private(exoid,"e",vecname,offsetZ,ierr))
      If (offsetN > 0) Then
         PetscCallA(MEF90EXOVecLoadNodal_Private(v,exoid,step+1,offsetN,ierr))
      Else If (offsetZ > 0) Then
         PetscCallA(MEF90EXOVecLoadZonal_Private(v,exoid,step+1,offsetZ,ierr))
      Else
         write(IOBuffer,'("Could not find nodal or zonal variable ", A5, " in exodus file. ")') vecname
         SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_UNEXPECTED,IOBuffer)
      End If
   End Subroutine MEF90EXOVecLoad

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOGetVarIndex_Private"
   Subroutine MEF90EXOGetVarIndex_Private(exoid,obj_type,name,varIndex,ierr)
      Integer,Intent(IN)               :: exoid
      Character(len=MXSTLN),Intent(IN) :: name
      Character(len=1),Intent(IN)      :: obj_type
      PetscInt,Intent(OUT)             :: varIndex
      PetscErrorCode,Intent(INOUT)     :: ierr
   
      Integer                         :: i, j, num_suffix = 5, num_vars
      Character(len=MXSTLN)            :: var_name,ext_name, suffix(5)
   
      suffix(1:5) = ["   ","_X ","_XX","_1 ","_11"]
      
      varIndex = -1
      PetscCallA(exgvp(exoid,obj_type,num_vars,ierr))
      Do i=1,num_vars
         PetscCallA(exgvnm(exoid,obj_type,i,var_name,ierr))
         Do j=1,num_suffix
            ext_name = trim(name)//trim(suffix(j))
            If (ext_name .eq. var_name) Then
               varIndex = i
            End If
         End Do
      End Do
   
   End Subroutine MEF90EXOGetVarIndex_Private
   
#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecViewNodal_Private"
   Subroutine MEF90EXOVecViewNodal_Private(v,exoid,step,offset,ierr)
      Integer,Intent(IN)               :: exoid
      PetscInt,Intent(IN)              :: step, offset
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
         PetscCall(expnvs(exoid,step,offset,xs+1,xe-xs,varray,ierr))
         PetscCall(VecRestoreArrayReadF90(v,varray,ierr))
      Else 
         PetscCall(ISCreateStride(PETSC_COMM_WORLD,(xe-xs)/bs,xs,bs,compIS,ierr))
         Do c = 0,bs-1
            PetscCall(ISStrideSetStride(compIS,(xe-xs)/bs,xs+c,bs,ierr))
            PetscCall(VecGetSubVector(v, compIS, vComp,ierr))
            PetscCall(VecGetArrayReadF90(vComp, varray,ierr))
            PetscCall(expnvs(exoid,step,offset+c,xs/bs+1,(xe-xs)/bs,varray,ierr))
            PetscCall(VecRestoreArrayReadF90(vComp, varray,ierr))
            PetscCall(VecRestoreSubVector(v, compIS, vComp,ierr))
         End Do
         PetscCall(ISDestroy(compIS,ierr))
      End If
   End Subroutine MEF90EXOVecViewNodal_Private

#undef __FUNCT__
#define __FUNCT__ "MEF90EXOVecLoadNodal_Private"
   Subroutine MEF90EXOVecLoadNodal_Private(v,exoid,step,offset,ierr)
      Integer,Intent(IN)               :: exoid
      PetscInt,Intent(IN)              :: step, offset
      Type(tVec),Intent(IN)            :: v
      PetscErrorCode,Intent(INOUT)     :: ierr
   
      PetscInt                         :: xs,xe,bs,c
      PetscScalar,Dimension(:),Pointer :: varray
      Type(tVec)                       :: vComp
      Type(tIS)                        :: compIS
   
      PetscCall(VecGetOwnershipRange(v,xs,xe,ierr))
      PetscCall(VecGetBlockSize(v,bs,ierr))
      If (bs == 1) Then
         PetscCall(VecGetArrayReadF90(v,varray,ierr));
         PetscCall(exgnnv(exoid,step,offset,xs+1,xe-xs,varray,ierr))
         PetscCall(VecRestoreArrayReadF90(v,varray,ierr))
      Else 
         PetscCall(ISCreateStride(PETSC_COMM_WORLD,(xe-xs)/bs,xs,bs,compIS,ierr))
         Do c = 0,bs-1
            PetscCall(ISStrideSetStride(compIS,(xe-xs)/bs,xs+c,bs,ierr))
            PetscCall(VecGetSubVector(v, compIS, vComp,ierr))
            PetscCall(VecGetArrayReadF90(vComp, varray,ierr))
            PetscCall(exgnnv(exoid,step,offset+c,xs/bs+1,(xe-xs)/bs,varray,ierr))
            PetscCall(VecRestoreArrayReadF90(vComp, varray,ierr))
            PetscCallA(VecISCopy(v,compIS,SCATTER_FORWARD,vComp,ierr))
            PetscCall(VecRestoreSubVector(v, compIS, vComp,ierr))
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
   
      PetscInt                         :: xs,xe,bs,c,numCS,set,csLocalSize,csxs=0
      PetscScalar,Dimension(:),Pointer :: varray
      PetscInt,Dimension(:),Pointer    :: csID,csSize
      Type(tVec)                       :: vComp
      Type(tIS)                        :: compIS
      Character(len=MXSTLN)            :: elemType
      PetscMPIInt                      :: rank
   
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      numCS = exinqi(exoid,EX_INQ_ELEM_BLK)
      Allocate(csID(numCS))
      Allocate(csSize(numCS))
      PetscCall(exgebi(exoid,csID,ierr))
      Do set=1,numCS
         PetscCall(exgelb(exoid,csID(set),elemType,csSize(set),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))
      End Do
      PetscCallA(VecGetOwnershipRange(v,xs,xe,ierr))
      PetscCallA(VecGetBlockSize(v,bs,ierr))
      If (bs > 1) Then
         PetscCallA(ISCreateStride(PETSC_COMM_WORLD,(xe-xs)/bs,xs,bs,compIS,ierr))
      End If
      Do set=1,numCS
         !  range of indices for set setID[set]: csxs:csxs + csSize[set]-1
         !  local slice of zonal values:         xs/bs,xm/bs-1
         !  intersection:                        max(xs/bs,csxs),min(xm/bs-1,csxs + csSize[set]-1)
         csLocalSize = max(0, min(xe/bs, csxs+csSize(set)) - max(xs/bs, csxs))
         If (bs == 1) Then
            PetscCall(VecGetArrayReadF90(v,varray,ierr))
            PetscCall(expevs(exoid,step,offset,csID(set),max(xs-csxs,0)+1,csLocalSize,varray,ierr))
            PetscCall(VecRestoreArrayReadF90(v,varray,ierr))
         Else
            Do c = 0,bs-1
               PetscCall(ISStrideSetStride(compIS,(xe-xs)/bs,xs+c,bs,ierr))
               PetscCall(VecGetSubVector(v,compIS,vComp,ierr))
               PetscCall(VecGetArrayReadF90(vComp,varray,ierr))
               PetscCall(expevs(exoid,step,offset+c,csID(set),max(xs/bs-csxs,0)+1,csLocalSize,varray,ierr))
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
      Type(tVec),Intent(IN)            :: v
      PetscErrorCode,Intent(INOUT)     :: ierr
   
      PetscInt                         :: xs,xe,bs,c,numCS,set,csLocalSize,csxs=0
      PetscScalar,Dimension(:),Pointer :: varray
      PetscInt,Dimension(:),Pointer    :: csID,csSize
      Type(tVec)                       :: vComp
      Type(tIS)                        :: compIS
      Character(len=MXSTLN)            :: elemType
      PetscMPIInt                      :: rank
   
      PetscCallMPIA(MPI_Comm_rank(PETSC_COMM_WORLD,rank,ierr))
      numCS = exinqi(exoid,EX_INQ_ELEM_BLK)
      Allocate(csID(numCS))
      Allocate(csSize(numCS))
      PetscCall(exgebi(exoid,csID,ierr))
      Do set=1,numCS
         PetscCall(exgelb(exoid,csID(set),elemType,csSize(set),PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr))
      End Do
      PetscCallA(VecGetOwnershipRange(v,xs,xe,ierr))
      PetscCallA(VecGetBlockSize(v,bs,ierr))
      If (bs > 1) Then
         PetscCallA(ISCreateStride(PETSC_COMM_WORLD,(xe-xs)/bs,xs,bs,compIS,ierr))
      End If
      Do set=1,numCS
         !  range of indices for set setID[set]: csxs:csxs + csSize[set]-1
         !  local slice of zonal values:         xs/bs,xm/bs-1
         !  intersection:                        max(xs/bs,csxs),min(xm/bs-1,csxs + csSize[set]-1)
         csLocalSize = max(0, min(xe/bs, csxs+csSize(set)) - max(xs/bs, csxs))
         If (bs == 1) Then
            PetscCall(VecGetArrayReadF90(v,varray,ierr))
            PetscCall(exgnev(exoid,step,offset,csID(set),csSize(set),max(xs-csxs,0)+1,csLocalSize,varray,ierr))
            PetscCall(VecRestoreArrayReadF90(v,varray,ierr))
         Else
            Do c = 0,bs-1
               PetscCall(ISStrideSetStride(compIS,(xe-xs)/bs,xs+c,bs,ierr))
               PetscCall(VecGetSubVector(v,compIS,vComp,ierr))
               PetscCall(VecGetArrayReadF90(vComp,varray,ierr))
               PetscCall(exgnev(exoid,step,offset+c,csID(set),csSize(set),max(xs/bs-csxs,0)+1,csLocalSize,varray,ierr))
               PetscCall(VecRestoreArrayReadF90(vComp,varray,ierr))
               PetscCallA(VecISCopy(v,compIS,SCATTER_FORWARD,vComp,ierr))
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
End Module m_MEF90_EXO