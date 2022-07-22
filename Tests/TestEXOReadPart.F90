Program TestExoRead_Part
#include <petsc/finclude/petsc.h>
#include "exodusii.h90"
   Use m_mef90
   Use petsc
   Implicit NONE   

   Integer                             :: exoid,cpu_ws,io_ws,mod_sz,exoerr
   Integer                             :: num_dim,num_nodes,num_elem,num_elem_blk,num_node_sets,num_side_sets
   Integer                             :: time_step = 1,var_index = 3,i,nparts,istart,iend,len
   Real,dimension(:),Pointer           :: var_values
   character(len=MXLNLN)               :: title
   Real                                :: vers
   PetscErrorCode                      :: ierr
   Type(MEF90Ctx_Type),target          :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default





   PetscCallA(PetscInitialize(PETSC_NULL_CHARACTER,ierr))
   Call MEF90Initialize(ierr)

   MEF90GlobalOptions_default%verbose           = 0
   MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
   MEF90GlobalOptions_default%timeMin           = 0.0_Kr
   MEF90GlobalOptions_default%timeMax           = 1.0_Kr
   MEF90GlobalOptions_default%timeNumStep       = 11


   Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)

   cpu_ws = 0
   io_ws = 0
   exoid = ex_open (MEF90Ctx%geometryfile, EXREAD, cpu_ws, io_ws, vers, exoerr)
   write (*, '("after exopen, error = ",i3)') exoerr

   write (*, '("test.exo is an EXODUSII file; version ", f4.2)') vers
   write (*, '("  I/O word size",i2)') io_ws
   write (*, '("  CPU word size",i2)') cpu_ws
   mod_sz = exlgmd(exoid)
   write (*, '("  Model Size",i2)') mod_sz

   call ex_get_init (exoid, title, num_dim, num_nodes, num_elem,num_elem_blk, num_node_sets, num_side_sets, exoerr)
   write (*, '("after exgini, error = ", i3)' ) exoerr

   write (*, '("database parameters")')
   write (*, '("   title = ",a81)') title
   write (*, '("   num_dim = ", i3 )') num_dim
   write (*, '("   num_nodes = ", i3 )') num_nodes
   write (*, '("   num_elem = ", i3 )') num_elem
   write (*, '("   num_elem_blk = ", i3 )') num_elem_blk
   write (*, '("   num_node_sets = ", i3 )') num_node_sets
   write (*, '("   num_side_sets = ", i3)') num_side_sets

   Allocate(var_values(num_nodes))
   call exgnv (exoid, time_step, var_index, num_nodes, var_values,exoerr)
   write (*, '("after exgnv, error = ", i3)' ) exoerr
   write(*,*) var_values
   DeAllocate(var_values)

   nparts=4
   Do i = 1, nparts
      istart = (i-1)*num_nodes/nparts+1
      iend   = i*num_nodes/nparts
      len    = iend - istart +1
      write(*,'("   chunk ",i3," -- ",i3, ":")',advance = 'no') istart,iend
      Allocate(var_values(len))
      call exgnnv(exoid, time_step, var_index,istart,len, var_values,exoerr)
      write(*,*) var_values
      DeAllocate(var_values)
   End Do

   call ex_close (exoid, exoerr)
   write (*, '("after exclos, error = ", i3)' ) exoerr
   PetscCallA(PetscFinalize(ierr))
End Program TestExoRead_Part
