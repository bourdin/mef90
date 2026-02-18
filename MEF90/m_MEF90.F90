module m_MEF90
#include "petsc/finclude/petsc.h"
#include "../mef90version.h"
   use m_MEF90_Ctx
   use m_MEF90_LinAlg
   use m_MEF90_Parameters
   use m_MEF90_baseClass
   use m_MEF90_Materials
   use m_MEF90_MPI
   use m_MEF90_EXO
   use m_MEF90_Utils
   use m_MEF90_DMPlex
   use m_MEF90_Elements
   use m_MEF90_MassMatrixInterface
   use m_MEF90_NormsInterface
   use m_MEF90_HookesLaw

   implicit none(type, external)

contains
#undef __FUNCT__
#define __FUNCT__ "MEF90Initialize"
!!!
!!!
!!!  MEF90Initialize:
!!!
!!!  (c) 2014-18 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90Initialize(comm, ierr)
      MPIU_Comm, intent(IN)                             :: comm
      PetscErrorCode, intent(INOUT)                     :: ierr

      character(len=MEF90MXSTRLEN)                      :: IOBuffer

      !Call PetscLogBegin(ierr);CHKERRQ(ierr)
      write (IOBuffer, *) "# MEF90: git changeset ", MEF90_GITVER, "\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# Copyright (c) 1998-2025 B. Bourdin <bourdin@mcmaster.ca> and co-authors\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# See CONTRIBUTORS.txt for a list of contributors\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# PETSC_ARCH=", PETSC_ARCH, "\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# PETSC_DIR=", PETSC_DIR, "\n\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))

      write (IOBuffer, *) "# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ""AS IS"" AND\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))

      write (IOBuffer, *) "# This software is released under the 2-clause BSD license (aka ""Simplified BSD"" \n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))
      write (IOBuffer, *) "# or ""FreeBSD"") license. See the LICENSE file in the root of the software distribution\n\n"
      PetscCall(PetscPrintf(comm, IOBuffer, ierr))

      !!! Individual modules runtime initialization should be called here
      PetscCall(MEF90MPIInitialize_Private(ierr))
      PetscCall(MEF90MaterialsInitialize_Private(ierr))
      PetscCall(MEF90ElementsInitialize_Private(ierr))
      PetscCall(MEF90CtxInitialize_Private(ierr))
      PetscCall(PetscLogDefaultBegin(ierr))
   end subroutine MEF90Initialize

#undef __FUNCT__
#define __FUNCT__ "MEF90Finalize"
!!!
!!!
!!!  MEF90Finalize:
!!!
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   subroutine MEF90Finalize(ierr)
      PetscErrorCode, intent(INOUT)                   :: ierr

      PetscCall(MEF90MPIFinalize_Private(ierr))
   end subroutine MEF90Finalize
end module m_MEF90

