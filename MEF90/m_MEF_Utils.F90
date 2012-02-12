Module m_MEF_Utils
#include "finclude/petscdef.h"
   Use m_MEF_Parameters
   Use petsc
   Implicit None

Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90_ISAllGatherMerge"
   !!! Merge all values of an IS, deleting duplicates
   Subroutine MEF90_ISAllGatherMerge(Comm,labels)
      MPI_Comm                         :: Comm
      Type(IS)                         :: labels
      
      Type(IS)                         :: tmplabels
      PetscInt, Dimension(:), Pointer  :: tmplabels_ptr
      PetscInt                         :: numval
      PetscErrorCode                   :: ierr
      
      Call ISAllGather(labels,tmplabels,ierr);CHKERRQ(ierr)
      Call ISGetSize(tmplabels,numval,ierr);CHKERRQ(ierr)
      Call ISGetIndicesF90(tmplabels,tmplabels_ptr,ierr);CHKERRQ(ierr)
      Call PetscSortRemoveDupsInt(numval,tmplabels_ptr,ierr);CHKERRQ(ierr)
      call ISDestroy(labels,ierr);CHKERRQ(ierr)
      Call ISCreateGeneral(Comm,numval,tmplabels_ptr,PETSC_COPY_VALUES,labels,ierr);CHKERRQ(ierr)
      Call ISRestoreIndicesF90(tmplabels,tmplabels_ptr,ierr);CHKERRQ(ierr)
      Call ISDestroy(tmplabels,ierr);CHKERRQ(ierr)
   End Subroutine MEF90_ISAllGatherMerge
   
   
   Subroutine Uniq(dComm, dMyVals, dVals)
      MPI_Comm                         :: dComm
      PetscInt, Dimension(:), Pointer  :: dMyVals, dVals
      
      Logical, Dimension(:), Pointer   :: ValCount
      PetscInt                         :: GlobMinVal, MyMinVal
      PetscInt                         :: GlobMaxVal, MyMaxVal
      PetscInt                         :: UniqCount
      PetscMPIInt                      :: rank
      PetscInt                         :: i, j, iErr
      PetscInt                         :: MySize, MaxSize

      Call MPI_Comm_Rank(PETSC_COMM_WORLD, rank, iErr)

      MySize = Size(dMyVals)
      Call MPI_AllReduce(MySize, MaxSize, 1, MPIU_INTEGER, MPI_MAX, dComm, iErr)
      If (MaxSize>0) Then
         MyMinVal = MinVal(dMyVals)
         MyMaxVal = MaxVal(dMyVals)
         Call MPI_AllReduce(MyMinVal, GlobMinVal, 1, MPIU_INTEGER, MPI_MIN, dComm, iErr)
         Call MPI_AllReduce(MyMaxVal, GlobMaxVal, 1, MPIU_INTEGER, MPI_MAX, dComm, iErr)
!!!$         openmpi doesn't like these so we do a Reduce followed by a Bcast
!!!$         Call MPI_Reduce(MyMinVal, GlobMinVal, 1, MPIU_INTEGER, MPI_MIN, 0, dComm, iErr)
!!!$         Call MIP_Bcast(GlobMinVal, 1, MPIU_INTEGER, 0, dComm, iErr)   
!!!$         Call MPI_Reduce(MyMaxVal, GlobMaxVal, 1, MPIU_INTEGER, MPI_MAX, 0, dComm, iErr)
!!!$         Call MIP_Bcast(GlobMaxVal, 1, MPIU_INTEGER, 0, dComm, iErr)   
!!!$   
         Allocate(ValCount(GlobMinVal:GlobMaxVal))
         ValCount = .FALSE.
         Do i = 1, Size(dMyVals)
            ValCount(dMyVals(i)) = .TRUE.
         End Do
   
         Call MPI_AllReduce(MPI_IN_PLACE, ValCount, GlobMaxVal-GlobMinVal+1, MPI_LOGICAL, MPI_LOR, dComm, iErr)
         !!! This is suboptimal. I could gather only to CPU 0 and do everything else on CPU 0 before broadcasting
         
         UniqCount = Count(ValCount)
   
         Allocate(dVals(UniqCount))
         j = 1
         Do i = GlobMinVal, GlobMaxVal
            If (ValCount(i)) Then
               dVals(j) = i
               j = j+1
            End If
         End Do
         DeAllocate(ValCount)
      Else
         Allocate(dVals(0))
      End If
   End Subroutine Uniq
   
   Subroutine GaussJordan_Inverse(A, Status)
!!!
!!! TODO
!!! Use Status correctly
!!!
      !
      ! Gauss Jordan inversion
      ! Very closely based on the routine from Numerical Recipes
      ! 
      PetscReal, Dimension(:,:), Pointer          :: A
      PetscInt, Intent(OUT)                       :: Status
      
      Integer, Dimension(:), Pointer              :: ipiv,indxr,indxc
      Logical, Dimension(:), Pointer              :: lpiv 
      Logical, Dimension(:,:), Pointer            :: lMask
      PetscReal                                   :: pivinv 
      PetscReal, Dimension(:),Pointer             :: dumc 
!      PetscReal, Dimension(:,:), Pointer          :: DumC2
      Integer, Target                             :: irc(2) 
      Integer                                     :: i,l,n 
      Integer, Pointer                            :: irow,icol 
      
      
      Status = 0
      N = Size(A,1)
      If (N /= Size(A,2) ) Then
         Write(*,*) 'Gauss Jordan: A is not square...'
         Status = -1
         Return
      End If
      Allocate (ipiv(N))
      Allocate (indxr(N))
      Allocate (indxc(N))
      Allocate (lpiv(N))
      Allocate (dumc(N))
      Allocate (lmask(N,N))
      
      
      irow =>irc(1) 
      icol =>irc(2) 
      ipiv=0 
      
      Do i = 1,n 
         lpiv = (ipiv == 0)
         Do l = 1, N
            lmask(:,l) = (ipiv == 0)
         End Do
         Do l = 1, N
            lmask(l,:) = (lmask(l,:) .AND. (ipiv == 0))
         End Do
         irc  = maxloc(abs(a), lmask)
         ipiv(icol) = ipiv(icol)+1 
         If (ipiv(icol) > 1) Then
            Print*, 'Singular Matrix'
            Status = -2
            Return
         End If
         If (irow /= icol) then 
            DumC = A(irow,:)
            A(irow,:) = A(icol,:)
            A(icol,:) = DumC
         End If
         indxr(i) = irow 
         indxc(i) = icol 
         If (A(icol,icol) == 0.0_Kr) Then
            Print*, 'Singular Matrix'
            Status = -2
            Return
         End If
      
         pivinv = 1.0_Kr / A(icol,icol) 
         A(icol,icol) = 1.0_Kr 
         A(icol,:) = A(icol,:)*pivinv 
         dumc = A(:,icol)
         A(:,icol) = 0.0
         
         A(icol,icol)  = pivinv 
         Do l = 1, N
            A(1:icol-1,l) = A(1:icol-1,l) - dumc(1:icol-1) * A(icol,l) 
            A(icol+1:,l)  = A(icol+1:,l)  - dumc(icol+1:)  * A(icol,l) 
         End Do
      End Do
            
      Do l = n,1,-1 
         DumC = A(:, indxr(l)) 
         A(:,indxr(l)) = A(:,indxc(l))
         A(:,indxc(l)) = DumC
      End Do
      
      DeAllocate (ipiv)
      DeAllocate (indxr)
      DeAllocate (indxc)
      DeAllocate (lpiv)
      DeAllocate (dumc)
      DeAllocate (lmask)
   End Subroutine GaussJordan_Inverse

   Subroutine GaussJordan_Solve(A, b, Status)
      !
      ! Gauss Jordan inversion
      ! Very closely based on the routine from Numerical recipes
      ! 
      PetscReal, Dimension(:,:), Pointer          :: A
      PetscReal, Dimension(:), Pointer            :: b
      PetscBool, Intent(OUT)                      :: Status
      
      Integer, Dimension(:), Pointer              :: ipiv,indxr,indxc 
      Logical, Dimension(:), Pointer              :: lpiv 
      Logical, Dimension(:,:), Pointer            :: lMask
      PetscReal                                   :: pivinv , DumR
      PetscReal, Dimension(:),Pointer             :: dumc 
      PetscReal, Dimension(:,:), Pointer          :: DumC2
      Integer, Target                             :: irc(2) 
      Integer                                     :: i,l,n 
      Integer, Pointer                            :: irow,icol 
    
    
      Status = .TRUE.
      N = Size(A,1)
      If (N /= Size(A,2) ) Then
       Write(*,*) 'Gauss Jordan: A is not square...'
       Status = .FALSE.
       Return
      End If
      Allocate (ipiv(N))
      Allocate (indxr(N))
      Allocate (indxc(N))
      Allocate (lpiv(N))
      Allocate (dumc(N))
      Allocate (lmask(N,N))
      
      
      irow =>irc(1) 
      icol =>irc(2) 
      ipiv=0 
      
      Do i = 1,n 
         lpiv = (ipiv == 0)
         Do l = 1, N
            lmask(:,l) = (ipiv == 0)
         End Do
         Do l = 1, N
            lmask(l,:) = (lmask(l,:) .AND. (ipiv == 0))
         End Do
         irc  = maxloc(abs(a), lmask)
         ipiv(icol) = ipiv(icol)+1 
         If (ipiv(icol) > 1) Then
            Print*, 'Singular Matrix'
            Status = .FALSE.
            Return
         End If
         If (irow /= icol) then 
            DumC = A(irow,:)
            A(irow,:) = A(icol,:)
            A(icol,:) = DumC
            DumR = b(irow)
            b(irow) = b(icol)
            b(icol) = DumR
         End If
         indxr(i) = irow 
         indxc(i) = icol 
         If (A(icol,icol) == 0.0_Kr) Then
            Print*, 'Singular Matrix'
            Status = .FALSE.
            Return
         End If
         
         pivinv = 1.0_Kr / A(icol,icol) 
         A(icol,icol) = 1.0_Kr 
         A(icol,:) = A(icol,:)*pivinv 
         b(icol) = b(icol) * pivinv
         dumc = A(:,icol)
         A(:,icol) = 0.0
         
         A(icol,icol)  = pivinv 
         Do l = 1, N
            A(1:icol-1,l) = A(1:icol-1,l) - dumc(1:icol-1) * A(icol,l) 
            A(icol+1:,l)  = A(icol+1:,l)  - dumc(icol+1:)  * A(icol,l) 
         End Do
         b(1:icol-1) = b(1:icol-1) - dumC(1:icol-1) * b(icol)
         b(icol+1:) = b(icol+1:) - dumC(icol+1:) * b(icol)
      End Do
      
      Do l = n,1,-1 
         DumC = A(:, indxr(l)) 
         A(:,indxr(l)) = A(:,indxc(l))
         A(:,indxc(l)) = DumC
      End Do
      
      DeAllocate (ipiv)
      DeAllocate (indxr)
      DeAllocate (indxc)
      DeAllocate (lpiv)
      DeAllocate (dumc)
      DeAllocate (lmask)
   End Subroutine GaussJordan_Solve
End Module m_MEF_Utils
