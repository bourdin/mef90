Module m_MEF90_Utils
#include "petsc/finclude/petsc.h"
   Use m_MEF90_Parameters
   Use m_MEF90_MPI
   Use petsc
   Implicit None

   Interface MEF90FindIndexOrdered
      Module Procedure MEF90FindIndexOrderedPetscInt,MEF90FindIndexOrderedPetscReal
   End Interface MEF90FindIndexOrdered
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90FindIndexOrderedPetscReal"

   Subroutine MEF90FindIndexOrderedPetscReal(x,array,pos)
      PetscReal,Intent(IN)                :: x
      PetscReal,Dimension(:),Pointer      :: array
      Integer,intent(OUT)                 :: pos
      
      Integer                             :: i1,i2,i
      
      i1 = lbound(array,1)
      i2 = ubound(array,1)
      
      If (array(i1) > x) Then
         pos = 0
      Else If (array(i1) == x) Then
         pos = i1
      Else If (array(i2) < x) Then
         pos = i2
      Else
         Do
            i = (i1+i2)/2
            If (array(i) == x) Then
               pos = i
               EXIT
            Else If (array(i) < x) Then
               i1 = i
            Else 
               i2 = i
            End If
   
            If (i2 == i1+1) Then
               pos = i1
               EXIT
            End If  
         End Do
      End If
   End Subroutine MEF90FindIndexOrderedPetscReal

#undef __FUNCT__
#define __FUNCT__ "MEF90FindIndexOrderedPetscInt"

   Subroutine MEF90FindIndexOrderedPetscInt(x,array,pos)
      PetscInt,Intent(IN)                 :: x
      PetscInt,Dimension(:),Pointer       :: array
      Integer,intent(OUT)                 :: pos
      
      Integer                             :: i1,i2,i
      
      i1 = lbound(array,1)
      i2 = ubound(array,1)
      
      If (array(i1) > x) Then
         pos = 0
      Else If (array(i1) == x) Then
         pos = i1
      Else If (array(i2) < x) Then
         pos = i2
      Else
         Do
            i = (i1+i2)/2
            If (array(i) == x) Then
               pos = i
               EXIT
            Else If (array(i) < x) Then
               i1 = i
            Else 
               i2 = i
            End If
   
            If (i2 == i1+1) Then
               pos = i1
               EXIT
            End If  
         End Do
      End If
   End Subroutine MEF90FindIndexOrderedPetscInt

#undef __FUNCT__
#define __FUNCT__ "MEF90ISAllGatherMerge"
!!!
!!!  
!!!  MEF90ISAllGatherMerge: Merge all values of an IS, deleting duplicates
!!!  
!!!  (c) 2013-2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90ISAllGatherMerge(Comm,labels,ierr)
      MPI_Comm,intent(IN)              :: Comm
      Type(tIS),intent(INOUT)          :: labels
      PetscErrorCode,intent(OUT)       :: ierr
      Type(tIS)                        :: tmplabels
      
      Call ISAllGather(labels,tmplabels,ierr);CHKERRQ(ierr)
      Call ISSortRemoveDups(tmplabels,ierr);CHKERRQ(ierr)
      Call ISDestroy(labels,ierr);CHKERRQ(ierr)
      labels%v = tmplabels%v
      
      ierr = 0
   End Subroutine MEF90ISAllGatherMerge
   
!#undef __FUNCT__
!#define __FUNCT__ "MEF90ISCreateCelltoVertex"
!!!!
!!!!  
!!!!  MEF90ISCreateCelltoVertex: Create an IS indexing all vertices associated with a cell set
!!!!                              Assume that the cone of an elements is made only of vertices
!!!!                              i.e. non  interpolated mesh
!!!!  
!!!!  (c) 2013-2014 Blaise Bourdin bourdin@lsu.edu
!!!!
!
!   Subroutine MEF90ISCreateCelltoVertex(mesh,comm,cellSetIS,vertexSetIS,ierr)
!      Type(tDM),Intent(IN)                             :: mesh
!      MPI_Comm,Intent(IN)                             :: Comm
!      Type(tIS),Intent(INOUT)                          :: cellSetIS,vertexSetIS
!      PetscErrorCode,Intent(OUT)                      :: ierr
!      
!      PetscInt,Dimension(:),Pointer                   :: cellSetIdx,vertexSetIdx,cone
!      PetscInt                                        :: numVertices,cell,v,coneSize
!      
!      Call ISGetIndicesF90(cellSetIS,cellSetIdx,ierr);CHKERRQ(ierr)
!      numVertices = 0
!      Do cell = 1,size(cellSetIdx)
!         Call DMMeshGetConeSize(mesh,cellSetIdx(cell),coneSize,ierr);CHKERRQ(ierr)
!         numVertices = numVertices + coneSize
!      End Do ! cell
!      Allocate(vertexSetIdx(numVertices))
!      Do cell = 1,size(cellSetIdx)
!         Call DMMeshGetConeSize(mesh,cellSetIdx(cell),coneSize,ierr);CHKERRQ(ierr)
!         Call DMMeshGetConeF90(mesh,cellSetIdx(cell),Cone,ierr);CHKERRQ(ierr)
!         Do v = 1,coneSize
!            vertexSetIdx((cell-1)*coneSize+v) = Cone(v)
!         End Do ! v
!         Call DMMeshRestoreConeF90(mesh,cellSetIdx(cell),Cone,ierr);CHKERRQ(ierr)
!      End Do ! cell
!      Call PetscSortRemoveDupsInt(numVertices,vertexSetIdx,ierr);CHKERRQ(ierr)
!      Call ISCreateGeneral(Comm,numVertices,vertexSetIdx,PETSC_COPY_VALUES,vertexSetIS,ierr);CHKERRQ(ierr)
!      DeAllocate(vertexSetIdx)
!      Call ISRestoreIndicesF90(cellSetIS,cellSetIdx,ierr);CHKERRQ(ierr)
!   End Subroutine MEF90ISCreateCelltoVertex
   
#undef __FUNCT__
#define __FUNCT__ "MEF90AskInt"

   Subroutine MEF90AskInt(val,msg,ArgUnit,IsBatch)
      PetscInt                                  :: Val
      Character(len=*)                          :: msg 
      PetscInt                                  :: argunit
      PetscBool                                 :: IsBatch

      Character(len=MEF90_MXSTRLEN)             :: IOBuffer   
      PetscInt                                  :: ierr   
      
      If (IsBatch) Then
         If (MEF90_MyRank == 0) Then
            Read(ArgUnit,*) Val
         End If
         Call MPI_BCast(Val,1,MPIU_INTEGER,0,PETSC_COMM_WORLD,ierr)
      Else
         Write(IOBuffer,"(A,t60,':  ')") Trim(msg)
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         If (MEF90_MyRank == 0) Then
            Read(*,*) Val
            Write(ArgUnit,"(I4,t60,A)") val,Trim(msg)
         End If
         Call MPI_BCast(Val,1,MPIU_INTEGER,0,PETSC_COMM_WORLD,ierr)
      End If
   End Subroutine MEF90AskInt   
   
#undef __FUNCT__
#define __FUNCT__ "MEF90AskReal"

   Subroutine MEF90AskReal(val,msg,ArgUnit,IsBatch)
      PetscReal                                 :: Val
      Character(len=*)                          :: msg 
      PetscInt                                  :: argunit
      PetscBool                                 :: IsBatch

      Character(len=MEF90_MXSTRLEN)             :: IOBuffer      
      PetscInt                                  :: ierr
      If (IsBatch) Then
         If (MEF90_MyRank == 0) Then
            Read(ArgUnit,*) Val
         End If
         Call MPI_BCast(Val,1,MPIU_SCALAR,0,PETSC_COMM_WORLD,ierr)
      Else
         Write(IOBuffer,"(A,t60,':  ')") Trim(msg)
         Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
         If (MEF90_MyRank == 0) Then
            Read(*,*) Val
            Write(ArgUnit,"(ES12.5,t60,A)") val,Trim(msg)
         End If
         Call MPI_BCast(Val,1,MPIU_SCALAR,0,PETSC_COMM_WORLD,ierr)
      End If
   End Subroutine MEF90AskReal
 

!!! 
!!! This should not be needed anywhere anymore
!!!   
!#undef __FUNCT__
!#define __FUNCT__ "MEF90Uniq"
!
!   Subroutine MEF90Uniq(dComm, dMyVals, dVals)
!      MPI_Comm                         :: dComm
!      PetscInt, Dimension(:), Pointer  :: dMyVals, dVals
!      
!      Logical, Dimension(:), Pointer   :: ValCount1,ValCount2
!      PetscInt                         :: GlobMinVal, MyMinVal
!      PetscInt                         :: GlobMaxVal, MyMaxVal
!      PetscInt                         :: UniqCount
!      PetscMPIInt                      :: rank
!      PetscInt                         :: i, j, iErr
!      PetscInt                         :: MySize, MaxSize
!
!      Call MPI_Comm_Rank(PETSC_COMM_WORLD, rank, iErr)
!
!      MySize = Size(dMyVals)
!      Call MPI_AllReduce(MySize, MaxSize, 1, MPIU_INTEGER, MPI_MAX, dComm, iErr)
!      If (MaxSize>0) Then
!         MyMinVal = MinVal(dMyVals)
!         MyMaxVal = MaxVal(dMyVals)
!         Call MPI_AllReduce(MyMinVal, GlobMinVal, 1, MPIU_INTEGER, MPI_MIN, dComm, iErr)
!         Call MPI_AllReduce(MyMaxVal, GlobMaxVal, 1, MPIU_INTEGER, MPI_MAX, dComm, iErr)
!         Allocate(ValCount1(GlobMinVal:GlobMaxVal))
!         Allocate(ValCount1(GlobMinVal:GlobMaxVal))
!         ValCount1 = .FALSE.
!         Do i = 1, Size(dMyVals)
!            ValCount1(dMyVals(i)) = .TRUE.
!         End Do
!
!         Call MPI_AllReduce(ValCount1, ValCount2, GlobMaxVal-GlobMinVal+1, MPI_LOGICAL, MPI_LOR, dComm, iErr)
!         !!! This is suboptimal. I could gather only to CPU 0 and do everything else on CPU 0 before broadcasting
!         
!         UniqCount = Count(ValCount2)
!   
!         Allocate(dVals(UniqCount))
!         j = 1
!         Do i = GlobMinVal, GlobMaxVal
!            If (ValCount2(i)) Then
!               dVals(j) = i
!               j = j+1
!            End If
!         End Do
!         DeAllocate(ValCount1)
!         DeAllocate(ValCount2)
!      Else
!         Allocate(dVals(0))
!      End If
!   End Subroutine MEF90Uniq

#undef __FUNCT__
#define __FUNCT__ "MEF90GetFilePrefix"

   subroutine  MEF90GetFilePrefix(instring, outstring)
      character(len=*), intent(IN)  :: instring
      character(len=*), intent(OUT) :: outstring

      integer                       :: index

      index = scan(instring,'.',back=.TRUE.)
      write(*,*) index
      outstring = instring(:index-1)
   end subroutine MEF90GetFilePrefix

   
End Module m_MEF90_Utils
