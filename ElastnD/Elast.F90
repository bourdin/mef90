Program  Elast

#include "finclude/petscdef.h"
#include "finclude/petscvecdef.h"
#include "finclude/petscmatdef.h"
#include "finclude/petsckspdef.h"
#include "finclude/petscviewerdef.h"
#include "finclude/petscmeshdef.h"

   Use m_MEF90
#if defined PB_2D
   Use m_Elast2D
#elif defined PB_3D 
   Use m_Elast3D
#endif

   Use petsc
   Use petscvec
   Use petscmat
   Use petscksp
   Use petscmesh

   Implicit NONE   


   Type(AppCtx_Type)                            :: AppCtx
   PetscInt                                     :: i, iErr
   Character(len=MEF90_MXSTRLEN)                :: IOBuffer
   PetscReal                                    :: rDummy
   Character                                    :: cDummy
   PetscInt                                     :: vers
   
   !!!
   PetscReal, Dimension(:), Pointer             :: CoordElem
   Type(SectionReal)                            :: CoordSec
   PetscInt                                     :: iE

   Call ElastInit(AppCtx)
   
   If (AppCtx%AppParam%verbose > 1) Then
      Call EXOView(AppCtx%EXO, AppCtx%AppParam%LogViewer)
      Call EXOView(AppCtx%MyEXO, AppCtx%AppParam%MyLogViewer)
      Call MeshTopologyView(AppCtx%MeshTopology, AppCtx%AppParam%MyLogViewer)
   End If   

   If (.NOT. AppCtx%VarFracSchemeParam%U_UseTao) Then
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Assembling the matrix\n'c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If   
      Call MatAssembly(AppCtx)
      If (AppCtx%AppParam%verbose > 1) Then
         Call MatView(AppCtx%KU, AppCtx%AppParam%LogViewer, iErr); CHKERRQ(iErr)
      End If
   End If
   
   !!! Testing Gradient function
!!$   Call SectionRealSet(AppCtx%U, 1.0_Kr, iErr)
   Allocate(CoordElem(6))
   Call MeshGetSectionReal(AppCtx%MeshTopology%mesh, 'coordinates', CoordSec, iErr); CHKERRQ(ierr)
   Do iE = 1, AppCtx%MeshTopology%Num_Elems
      Call MeshRestrictClosure(AppCtx%MeshTopology%mesh, CoordSec, iE-1, 6, CoordElem, iErr); CHKERRQ(iErr)
!      CoordElem = (CoordElem+1.0)/2.0
      CoordElem(2)=0.0; CoordElem(4) = 0.0; CoordElem(6) = 0.0
      CoordElem(1) = 1.0-CoordElem(1)**2; CoordElem(3) = 1.0-CoordElem(3)**2; CoordElem(5) = 1.0-CoordElem(5)**2
      Call MeshUpdateClosure(AppCtx%MeshTopology%Mesh, AppCtx%U, iE-1, CoordElem, iErr); CHKERRQ(iErr) 
   End Do
   
   Call SectionRealToVec(AppCtx%U, AppCtx%ScatterVect, SCATTER_FORWARD, AppCtx%U_Vec, ierr); CHKERRQ(ierr)
 Call FormFunctionandGradient(AppCtx%taoU, AppCtx%U_Vec, rDummy, AppCtx%RHSU, AppCtx, iErr)
     
!  Call VecView(AppCtx%RHSU, PETSC_VIEWER_STDOUT_WORLD, iErr); CHKERRQ(iErr)
   Write(*,*) 'func', rdummy
   Call VecDot(AppCtx%RHSU, AppCtx%U_Vec, rdummy, ierr)
   write(*,*)'==========', rdummy
  
!  Call ElastFinalize(AppCtx)
!  STOP 


!!!$   Do i = 1, AppCtx%NumTimeSteps
i=1
      AppCtx%TimeStep = i
      AppCtx%MyEXO%exoid = EXOPEN(AppCtx%MyEXO%filename, EXREAD, exo_cpu_ws, exo_io_ws, vers, iErr)
      Call EXGTIM(AppCtx%MyEXO%exoid, i, AppCtx%Time, iErr)
      Call EXCLOS(AppCtx%MyEXO%exoid, iErr)
      AppCtx%MyEXO%exoid = 0
      Call Read_EXO_Result_Global(AppCtx%MyEXO, AppCtx%MyEXO%GlobVariable(VarFrac_GlobVar_Load)%Offset, AppCtx%TimeStep, AppCtx%Load)

      !!! Read U, F, and Temperature
!!!$      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_DisplacementX)%Offset, AppCtx%TimeStep, AppCtx%U) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_ForceX)%Offset, AppCtx%TimeStep, AppCtx%F) 
      Call Read_EXO_Result_Vertex(AppCtx%MyEXO, AppCtx%MeshTopology, AppCtx%MyEXO%VertVariable(VarFrac_VertVar_Temperature)%Offset, AppCtx%TimeStep, AppCtx%Theta) 
   
      Call Solve(AppCtx)
   
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Computing Elastic energy, strains and stresses\n'c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If
      Call ComputeEnergy(AppCtx)
      Write(IOBuffer, 108) AppCtx%TimeStep, AppCtx%Time, AppCtx%Load, AppCtx%ElasticEnergy, AppCtx%ExtForcesWork, AppCtx%ElasticEnergy - AppCtx%ExtForcesWork
108 Format('TS ',I4, ' Time:', ES10.3, ' Load:', ES10.3, ' Elast:', ES10.3, ' Work:', ES10.3, ' Total:', ES10.3, '\n'c)
      Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)

      Write(IOBuffer, 110) AppCtx%TimeStep, AppCtx%Time, AppCtx%Load, AppCtx%ElasticEnergy, AppCtx%ExtForcesWork, AppCtx%ElasticEnergy - AppCtx%ExtForcesWork
110 Format(I4, 5(ES13.5, '   '), '\n'c)
      Call PetscViewerASCIIPrintf(AppCtx%AppParam%EnergyViewer, IOBuffer, iErr); CHKERRQ(iErr)

      If ( (AppCtx%VarFracSchemeParam%SaveStress) .OR. ( AppCtx%VarFracSchemeParam%SaveStrain) ) Then
         Call ComputeStrainStress(AppCtx)
      End If   
      If (AppCtx%AppParam%verbose > 0) Then
         Write(IOBuffer, *) 'Saving results\n'c
         Call PetscPrintf(PETSC_COMM_WORLD, IOBuffer, iErr); CHKERRQ(iErr)
      End If      
      Call Save(AppCtx)
!!!$   End Do

   Call ElastFinalize(AppCtx)
End Program  Elast
