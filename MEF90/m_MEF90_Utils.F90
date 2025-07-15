module m_MEF90_Utils
#include "petsc/finclude/petsc.h"
   use petscis
   use m_MEF90_Parameters
   implicit none(type, external)
   private
   public :: MEF90ISAllGatherMerge
   public :: MEF90FilePrefix
   public :: MEF90FileExtension
   public :: MEF90StrCount
   public :: MEF90StrTokenize

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90ISAllGatherMerge"
!!!
!!!
!!!  MEF90ISAllGatherMerge: Merge all values of an IS, deleting duplicates
!!!
!!!  (c) 2013-2014 Blaise Bourdin bourdin@lsu.edu
!!!      2022      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90ISAllGatherMerge(Comm, is, ierr)
      MPI_Comm, intent(IN)              :: Comm
      type(tIS), intent(INOUT)          :: is
      PetscErrorCode, intent(INOUT)     :: ierr

      type(tIS)                        :: tmpIS
      PetscInt, dimension(:), pointer    :: indices
      PetscInt                         :: sz

      PetscCall(ISGetIndices(is, indices, ierr))
      sz = size(indices)
      PetscCall(ISCreateGeneral(Comm, sz, indices, PETSC_COPY_VALUES, tmpIS, ierr))
      PetscCall(ISRestoreIndices(is, indices, ierr))
      PetscCall(ISDestroy(is, ierr))
      PetscCall(ISAllGather(tmpIS, is, ierr))
      PetscCall(ISDestroy(tmpIS, ierr))
      PetscCall(ISSortRemoveDups(is, ierr))
      ierr = 0
   end subroutine MEF90ISAllGatherMerge

#undef __FUNCT__
#define __FUNCT__ "MEF90FilePrefix"
   function MEF90FilePrefix(s)
      character(len=*), intent(IN)   :: s
      character(len=MEF90MXSTRLEN)  :: MEF90FilePrefix

      character(len=MEF90MXSTRLEN)  :: sChop
      integer                       :: i, l

      l = len(s)
      if (l > MEF90MXSTRLEN) then
         sChop = s(1:MEF90MXSTRLEN)
      else
         sChop = s
      end if
      i = index(sChop, '.', .true.)
      if (i == 0) then
         MEF90FilePrefix = sChop
      else
         MEF90FilePrefix = sChop(1:i - 1)
      end if
   end function MEF90FilePrefix

#undef __FUNCT__
#define __FUNCT__ "MEF90FileExtension"
   function MEF90FileExtension(s)
      character(len=*), intent(IN)   :: s
      character(len=MEF90MXSTRLEN)  :: MEF90FileExtension

      character(len=MEF90MXSTRLEN)  :: sChop
      integer                       :: i, l

      l = len(s)
      if (l > MEF90MXSTRLEN) then
         sChop = s(1:MEF90MXSTRLEN)
      else
         sChop = s
      end if
      i = index(sChop, '.', .true.)
      if (i == 0) then
         MEF90FileExtension = ''
      else
         MEF90FileExtension = sChop(i + 1:)
      end if
   end function MEF90FileExtension

#undef __FUNCT__
#define __FUNCT__ "MEF90StrCount"
!!!
!!!
!!!  MEF90StrCount: counts the number of occurrences of a substring in a string
!!!
!!!  (c) 2025      Blaise Bourdin bourdin@mcmaster.ca
!!!

   function MEF90StrCount(s, c)
      character(len=*), intent(IN)   :: s
      character(len=*), intent(IN)   :: c
      integer                       :: MEF90StrCount

      integer                       :: idx, pos, count
      count = 0
      pos = 1
      do
         idx = index(s(pos:), c)
         if (idx == 0) then
            exit
         else
            count = count + 1
            pos = pos + idx + len(c)
         end if
      end do
      MEF90StrCount = count
   end function MEF90StrCount

#undef __FUNCT__
#define __FUNCT__ "MEF90StrTokenize"
!!!
!!!
!!!  MEF90StrTokenize: Breaks a string into an array of substrings based on a delimiter
!!!
!!!  (c) 2025      Blaise Bourdin bourdin@mcmaster.ca
!!!

   subroutine MEF90StrTokenize(str, delim, tokens)
      character(len=*), intent(IN)                           :: str
      character(len=*), intent(IN)                           :: delim
      character(len=MEF90MXSTRLEN), dimension(:), allocatable :: tokens

      integer                                               :: i, n, pos, offset
      n = MEF90StrCount(str, delim)
      allocate (tokens(n + 1))
      if (n == 0) then
         tokens(1) = str
      else
         pos = 1
         do i = 1, n
            offset = index(str(pos:), delim)
            tokens(i) = str(pos:pos + offset - 2)
            pos = pos + offset - 1 + len(delim)
         end do
         tokens(n + 1) = str(pos:)
      end if
   end subroutine MEF90StrTokenize
end module m_MEF90_Utils
