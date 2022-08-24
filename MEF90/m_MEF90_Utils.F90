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
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   Subroutine MEF90ISAllGatherMerge(Comm,is,ierr)
      MPI_Comm,intent(IN)              :: Comm
      Type(tIS),intent(INOUT)          :: is
      PetscErrorCode,intent(INOUT)     :: ierr

      Type(tIS)                        :: tmpIS
      PetscInt,Dimension(:),pointer    :: indices
      PetscInt                         :: sz
      
      PetscCall(ISGetIndicesF90(is,indices,ierr))
      sz = size(indices)
      PetscCall(ISCreateGeneral(Comm,sz,indices,PETSC_COPY_VALUES,tmpIS,ierr))
      PetscCall(ISRestoreIndicesF90(is,indices,ierr))
      PetscCall(ISDestroy(is,ierr))
      PetscCall(ISAllGather(tmpIS,is,ierr))
      PetscCall(ISDestroy(tmpIS,ierr))
      PetscCall(ISSortRemoveDups(is,ierr))
      ierr = 0
   End Subroutine MEF90ISAllGatherMerge
   

#undef __FUNCT__
#define __FUNCT__ "MEF90FilePrefix"
   function MEF90FilePrefix(s)
      character(len=*),intent(IN)   :: s
      character(len=MEF90MXSTRLEN) :: MEF90FilePrefix

      character(len=MEF90MXSTRLEN) :: sChop
      integer                       :: i,l

      l = len(s)
      if (l > MEF90MXSTRLEN) then
         !write(*,*) 'Warning, choping input string'
         sChop = s(1:MEF90MXSTRLEN)
      else 
         sChop = s
      end if
      i = index(sChop,'.',.TRUE.)
      if (i ==  0) then
         MEF90FilePrefix = sChop
      else
         MEF90FilePrefix = sChop(1:i-1)
      end if
   end function MEF90FilePrefix


#undef __FUNCT__
#define __FUNCT__ "MEF90FileExtension"
   function MEF90FileExtension(s)
      character(len=*),intent(IN)   :: s
      character(len=MEF90MXSTRLEN) :: MEF90FileExtension

      character(len=MEF90MXSTRLEN) :: sChop
      integer                       :: i,l

      l = len(s)
      if (l > MEF90MXSTRLEN) then
         !write(*,*) 'Warning, choping input string'
         sChop = s(1:MEF90MXSTRLEN)
      else 
         sChop = s
      end if
      i = index(sChop,'.',.TRUE.)
      if (i ==  0) then
         MEF90FileExtension = ''
      else
         MEF90FileExtension = sChop(i+1:)
      end if
   end function MEF90FileExtension


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
!      PetscCallMPI(MPI_Comm_Rank(PETSC_COMM_WORLD, rank, iErr))
!
!      MySize = Size(dMyVals)
!      PetscCallMPI(MPI_AllReduce(MySize, MaxSize, 1, MPIU_INTEGER, MPI_MAX, dComm, iErr))
!      If (MaxSize>0) Then
!         MyMinVal = MinVal(dMyVals)
!         MyMaxVal = MaxVal(dMyVals)
!         PetscCallMPI(MPI_AllReduce(MyMinVal, GlobMinVal, 1, MPIU_INTEGER, MPI_MIN, dComm, iErr))
!         PetscCallMPI(MPI_AllReduce(MyMaxVal, GlobMaxVal, 1, MPIU_INTEGER, MPI_MAX, dComm, iErr))
!         Allocate(ValCount1(GlobMinVal:GlobMaxVal))
!         Allocate(ValCount1(GlobMinVal:GlobMaxVal))
!         ValCount1 = .FALSE.
!         Do i = 1, Size(dMyVals)
!            ValCount1(dMyVals(i)) = .TRUE.
!         End Do
!
!         PetscCallMPI(MPI_AllReduce(ValCount1, ValCount2, GlobMaxVal-GlobMinVal+1, MPI_LOGICAL, MPI_LOR, dComm, iErr))
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
