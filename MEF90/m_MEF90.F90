Module m_MEF90
#include "finclude/petscdef.h"
#include "../mef90version.h"
   Use petsc
   Use m_MEF90_Ctx
   Use m_MEF90_DiffusionInterface
   Use m_MEF90_ElasticityInterface
   Use m_MEF90_GradDamageInterface
   Use m_MEF90_Elements 
   Use m_MEF90_EXO  
   Use m_MEF90_LinAlg
   Use m_MEF90_Parameters
   Use m_MEF90_MassMatrixInterface
   Use m_MEF90_Materials
   Use m_MEF90_MPI
   Use m_MEF90_Norm
   Use m_MEF90_Utils

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
      PetscInt,Intent(OUT)                               :: ierr

      Character(len=MEF90_MXSTRLEN)                      :: IOBuffer
       
      Call PetscLogBegin(ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# MEF90: hg changeset ",MEF90_HGVER,"\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# Copyright (c) 1998-2018 B. Bourdin <bourdin@lsu.edu> and co-authors\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# See CONTRIBUTORS.txt for a list of contributors\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# PETSC_ARCH=", PETSC_ARCH ,"\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# PETSC_DIR=", PETSC_DIR ,"\n\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      
      Write(IOBuffer,*) "# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ""AS IS"" AND\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      
      Write(IOBuffer,*) "# This software is released under the 2-clause BSD license (aka ""Simplified BSD"" \n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      Write(IOBuffer,*) "# or ""FreeBSD"") license. See the LICENSE file in the root of the software distribution\n\n"
      Call PetscPrintf(PETSC_COMM_WORLD,IOBuffer,ierr);CHKERRQ(ierr)
      
      !!! Individual modules runtime initialization should be called here
      Call MEF90MPIInitialize_Private(ierr);CHKERRQ(ierr)
      Call MEF90MaterialsInitialize_Private(ierr);CHKERRQ(ierr)
      Call MEF90CtxInitialize_Private(ierr);CHKERRQ(ierr)
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
      PetscInt,Intent(OUT)                   :: ierr
      
      Call MEF90MPIFinalize_Private(ierr);CHKERRQ(ierr)
   End Subroutine MEF90Finalize
End Module m_MEF90

