Module m_MEF90
#include "petsc/finclude/petsc.h"
#include "../mef90version.h"
   Use petsc
   Use m_MEF90_Ctx
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_Materials
   Use m_MEF90_MPI
   Use m_MEF90_EXO  
   Use m_MEF90_Utils
   Use m_MEF90_DMPlex
   !Use m_MEF90_DiffusionInterface
   !Use m_MEF90_ElasticityInterface
   !Use m_MEF90_GradDamageInterface
   Use m_MEF90_Elements 
   Use m_MEF90_MassMatrixInterface
   !Use m_MEF90_Norm

   Implicit NONE
   Public :: MEF90Initialize
   Public :: MEF90Finalize
   
   
Contains
#undef __FUNCT__
#define __FUNCT__ "MEF90Initialize"
!!!
!!!  
!!!  MEF90Initialize:
!!!  
!!!  (c) 2014-18 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90Initialize(ierr)
      PetscErrorCode,Intent(INOUT)                      :: ierr

      Character(len=MEF90MXSTRLEN)                      :: IOBuffer
       
      !Call PetscLogBegin(ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# MEF90: git changeset ",MEF90_GITVER,"\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# Copyright (c) 1998-2022 B. Bourdin <bourdin@mcmaster.ca> and co-authors\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# See CONTRIBUTORS.txt for a list of contributors\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# PETSC_ARCH=", PETSC_ARCH ,"\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# PETSC_DIR=", PETSC_DIR ,"\n\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))

      Write(IOBuffer,*) "# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ""AS IS"" AND\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      
      Write(IOBuffer,*) "# This software is released under the 2-clause BSD license (aka ""Simplified BSD"" \n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))
      Write(IOBuffer,*) "# or ""FreeBSD"") license. See the LICENSE file in the root of the software distribution\n\n"
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr))

      !!! Individual modules runtime initialization should be called here
      PetscCall(MEF90MPIInitialize_Private(ierr))
      PetscCall(MEF90MaterialsInitialize_Private(ierr))
      PetscCall(MEF90ElementsInitialize_Private(ierr))
      PetscCall(MEF90CtxInitialize_Private(ierr))
   End Subroutine MEF90Initialize
   
   
#undef __FUNCT__
#define __FUNCT__ "MEF90Finalize"
!!!
!!!  
!!!  MEF90Finalize:
!!!  
!!!  (c) 2014 Blaise Bourdin bourdin@lsu.edu
!!!

   Subroutine MEF90Finalize(ierr)
      PetscErrorCode,Intent(INOUT)                   :: ierr
      
      PetscCall(MEF90MPIFinalize_Private(ierr))
   End Subroutine MEF90Finalize
End Module m_MEF90

