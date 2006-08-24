Program Elast2DA

  Use m_MEF90
  Use m_Rupt_Struct


  Use m_Elast2DA_Vars
  Use m_Elast2DA_Proc
  Implicit NONE

#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscmat.h90"
#include "include/finclude/petscksp.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscao.h"


  Integer               :: iE, iS, NbIter
  PetscScalar           :: ener
  PetscScalar           :: value

  PetscReal             :: One = 1.0_Kr
  PetscReal             :: Zero = 0.0_Kr
  PetscReal		:: delta_t = 0.05_Kr
  PetscReal             :: alpha = -0.2_Kr

  Integer                                             :: Nb_Gauss, Nb_DoF
    Integer                                             :: iSL1, iSG1
    Integer                                             :: iSL2, iSG2
    Integer                                             :: iELoc, iG
    Integer                                             :: iBlk
    

  Real(Kind = Kr), Dimension(:), Pointer       :: BC_Ptr
  Real(Kind = Kr), Dimension(:), Pointer       :: EXO_Ptr

  Vec 						:: Temp  	  
  Vec                                           :: Temp2
	
 
  Mat                   :: Loc_Mat

  Call PetscGetTime(TotalTS, iErr)
  Call PetscGetTime(InitTS, iErr)

  Call Init()

  Call PetscGetTime(InitTF, iErr)

  Call PetscGetTime(AssembTS, iErr)

  Call Assemb_Mat_Elast2DA(MR, Geom, Params, Elem_db, Node_db)

  Call PetscGetTime(AssembTF, iErr)

  Call Init_KSP()

  Write(CharBuffer,*) 'Number of timesteps:                        ',      &
       & Size(Params%Load), '\n'c
  Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)


     Call VecGetArrayF90(BC, BC_Ptr, iErr)
     Call Read_EXO_Result_Nodes(Geom, 4, 1, EXO_Ptr, 1)
     BC_Ptr = EXO_Ptr
     Call VecRestoreArrayF90(BC, BC_Ptr, iErr)
     DeAllocate(EXO_Ptr)

! *** initialize the solution vectors
Call VecDuplicate(BC, U, iErr)
Call VecCopy(BC, U, iErr)
Call VecDuplicate(U, U_old, iErr)
Call VecCopy(U, U_old, iErr)
Call VecDuplicate(U, U_t, iErr)
Call VecDuplicate(U, U_t_old, iErr)
Call VecZeroEntries(U_t, iErr)
Call VecZeroEntries(U_t_old, iErr)
Call VecDuplicate(U, Temp, iErr)
Call VecDuplicate(U, Temp2, iErr)

Call MatScale(MR, -1.0_Kr, iErr)

  Do TimeStep = 1, Size(Params%Load)

     Write(CharBuffer,*) '=== processing time step                 ',         &
          & TimeStep, '\n'c
     Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)

   ! *** Read in the boundary condition
     Call VecGetArrayF90(BC, BC_Ptr, iErr)
     Call Read_EXO_Result_Nodes(Geom, 4, TimeStep, EXO_Ptr, 1)
     BC_Ptr = EXO_Ptr


    Do iBlk = 1, Geom%Num_elem_blks
      Do iELoc = 1, Geom%Elem_blk(iBlk)%Num_Elems
          iE = Geom%Elem_blk(iBlk)%Elem_ID(iELoc)
          Do iSL1 = 1, Elem_db(iE)%Nb_DoF      
             iSG1 = Elem_db(iE)%ID_DoF(iSL1)             
             If ( Node_db(iSG1)%BC /= BC_Type_NONE ) Then
                Call VecSetValue(U_old, iSG1-1, BC_Ptr(iSG1), INSERT_VALUES, iErr)
             End if
          End Do 
       End Do 
    End Do 


     Call VecRestoreArrayF90(BC, BC_Ptr, iErr)
     DeAllocate(EXO_Ptr)
   
	Call VecWAXPY(U, delta_t, U_t_old, U_old, iErr)

	Call MatMult(MR, U_old, Temp, iErr)
        Call VecWAXPY(Temp2, alpha, U_t, Temp, iErr)
	Call VecWAXPY(U_t, delta_t, Temp2, U_t_old, iErr)

	
	Call VecCopy(U, U_old, iErr)
	Call VecCopy(U_t, U_t_old, iErr)
	Call Export()
  End Do
 

  Do TimeStep = Size(Params%Load) + 1, 3000

     Write(CharBuffer,*) '=== processing time step                 ',         &
          & TimeStep, '\n'c
     Call PetscPrintf(PETSC_COMM_WORLD, CharBuffer, iErr)
   
	Call VecWAXPY(U, delta_t, U_t_old, U_old, iErr)

	Call MatMult(MR, U_old, Temp, iErr)
        Call VecWAXPY(Temp2, alpha, U_t, Temp, iErr)
	Call VecWAXPY(U_t, delta_t, Temp2, U_t_old, iErr)

	
	Call VecCopy(U, U_old, iErr)
	Call VecCopy(U_t, U_t_old, iErr)
	Call Export()
  End Do
 
! *** clean up the vectors
Call VecDestroy(U ,iErr)
Call VecDestroy(U_t ,iErr)
Call VecDestroy(U_old ,iErr)
Call VecDestroy(U_t_old ,iErr)
Call VecDestroy(Temp ,iErr)

  Call PetscGetTime(TotalTF, iErr)
  Write(CharBuffer,*) 'Total time                                  ',      &
       & TotalTF - TotalTS, '\n'c
  Call PetscPrintf(PETSC_COMM_SELF, CharBuffer, iErr)

  Call Finalize()

  100 Format(A)
End Program Elast2DA
