Module m_Utils

   Use m_Constantes
   Implicit None
   Private
   
   Public :: GaussJordan_Inverse
   Public :: GaussJordan_Solve1

 Contains
   Subroutine GaussJordan_Inverse(A, Status)
      !
      ! Gauss Jordan inversion
      ! Very closely based on the routine from Numerical recipes
      ! 
      Implicit NONE 
      Real(Kind = Kr), Dimension(:,:), Pointer    :: A
      Logical, Intent(OUT)                        :: Status
      
      Integer, Dimension(:), Pointer              :: ipiv,indxr,indxc 
      Logical, Dimension(:), Pointer              :: lpiv 
      Logical, Dimension(:,:), Pointer            :: lMask
      Real(Kind = Kr)                             :: pivinv 
      Real(Kind = Kr), Dimension(:),Pointer       :: dumc 
      Real(Kind = Kr), Dimension(:,:), Pointer    :: DumC2
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

   Subroutine GaussJordan_Solve1(A, b, Status)
      !
      ! Gauss Jordan inversion
      ! Very closely based on the routine from Numerical recipes
      ! 
      Implicit NONE 
      Real(Kind = Kr), Dimension(:,:), Pointer    :: A
      Real(Kind = Kr), Dimension(:), Pointer      :: b
      Logical, Intent(OUT)                        :: Status
      
      Integer, Dimension(:), Pointer              :: ipiv,indxr,indxc 
      Logical, Dimension(:), Pointer              :: lpiv 
      Logical, Dimension(:,:), Pointer            :: lMask
      Real(Kind = Kr)                             :: pivinv , DumR
      Real(Kind = Kr), Dimension(:),Pointer       :: dumc 
      Real(Kind = Kr), Dimension(:,:), Pointer    :: DumC2
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
   End Subroutine GaussJordan_Solve1
End Module m_Utils
