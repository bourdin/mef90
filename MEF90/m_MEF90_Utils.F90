Module m_MEF90_Utils
#include "petsc/finclude/petsc.h"
   Use petscis
   Use m_MEF90_Parameters
   Implicit None
   private
   public :: MEF90ISAllGatherMerge
   public :: MEF90FilePrefix
   public :: MEF90FileExtension

Contains
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
      
      PetscCall(ISGetIndices(is,indices,ierr))
      sz = size(indices)
      PetscCall(ISCreateGeneral(Comm,sz,indices,PETSC_COPY_VALUES,tmpIS,ierr))
      PetscCall(ISRestoreIndices(is,indices,ierr))
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
      character(len=MEF90MXSTRLEN)  :: MEF90FilePrefix

      character(len=MEF90MXSTRLEN)  :: sChop
      integer                       :: i,l

      l = len(s)
      if (l > MEF90MXSTRLEN) then
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
      character(len=MEF90MXSTRLEN)  :: MEF90FileExtension

      character(len=MEF90MXSTRLEN)  :: sChop
      integer                       :: i,l

      l = len(s)
      if (l > MEF90MXSTRLEN) then
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

End Module m_MEF90_Utils
