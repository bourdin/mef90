Program  TestDMPlexComputeGeometry
#include <petsc/finclude/petsc.h>
#include <petsc/finclude/petscdmlabel.h>
   Use m_MEF90
   Use petsc
   Implicit NONE   

   PetscErrorCode                      :: ierr
   Type(MEF90Ctx_Type),target          :: MEF90Ctx
   Type(MEF90CtxGlobalOptions_Type)    :: MEF90GlobalOptions_default
   Type(tDM),target                    :: dm,dmDist
   PetscBool                           :: interpolate = PETSC_TRUE
   Character(len=MEF90_MXSTRLEN)       :: IOBuffer
   PetscInt                            :: set,cell
   DMLabel                             :: CSLabel,FSLabel
   type(tIS)                           :: CSIS,FSIS,CellIS
   PetscInt,Dimension(:),pointer       :: setID,cellID
   PetscReal,Dimension(:),pointer      :: v0
   PetscReal,Dimension(:),pointer      :: B,Binv
   PetscReal,Dimension(:),pointer      :: detB
   PetscInt                            :: dim,nquad
   PetscQuadrature                     :: q
   PetscReal                           :: vol
   PetscReal,Dimension(:),pointer      :: centroid,normal

   MEF90GlobalOptions_default%verbose           = 1
   MEF90GlobalOptions_default%dryrun            = PETSC_FALSE
   MEF90GlobalOptions_default%timeInterpolation = MEF90TimeInterpolation_linear
   MEF90GlobalOptions_default%timeMin           = 0.0_Kr
   MEF90GlobalOptions_default%timeMax           = 1.0_Kr
   MEF90GlobalOptions_default%timeFrequency     = 0.0_Kr
   MEF90GlobalOptions_default%timeNumStep       = 11
   MEF90GlobalOptions_default%fileFormat        = MEF90FileFormat_EXOSingle

   Call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
   Call MEF90Initialize(ierr)
   Call MEF90CtxCreate(PETSC_COMM_WORLD,MEF90Ctx,MEF90GlobalOptions_default,ierr)

   Call DMPlexCreateFromFile(PETSC_COMM_WORLD,MEF90Ctx%inputmesh,interpolate,dm,ierr);CHKERRQ(ierr);
   Call DMPlexDistribute(dm,0,PETSC_NULL_SF,dmDist,ierr);CHKERRQ(ierr)

   !!!
   !!! I am not sure everybody would approve of this...
   !!!
   if (dmDist%v /= -1) then
      call DMDestroy(dm,ierr)
      dm%v = dmDist%v
   end if

   Call DMView(dm,PETSC_VIEWER_STDOUT_WORLD,ierr)
   call DMGetDimension(dm,dim,ierr)
   nquad=2

   call PetscQuadratureCreate(PETSC_COMM_SELF,q,ierr)
   call PetscDTGaussTensorQuadrature(dim,1,nquad,-1.0_Kr,1.0_Kr,q,ierr)
   call PetscQuadratureView(q,PETSC_VIEWER_STDOUT_SELF,ierr)


   allocate(v0(nquad*dim*dim))
   allocate(B(nquad*dim*dim**2))
   allocate(Binv(nquad*dim*dim**2))
   allocate(detB(nquad*dim))
   call DMGetLabel(dm,"Cell Sets",CSLabel,ierr)
   call DMLabelGetValueIS(CSLabel,CSIS,ierr)
   call ISGetIndicesF90(CSIS,setID,ierr)
   do set = 1, size(setID)
      write(IOBuffer,*) 'Cell Set', setID(set),'\n'
      call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr)
      call DMLabelGetStratumIS(CSLabel,setID(set),cellIS,ierr)
      call ISGetIndicesF90(cellIS,cellID,ierr)
      do cell = 1, size(cellID)
         call DMPlexComputeCellGeometryFEM(dm,cellID(cell),q,v0,B,Binv,detB,ierr)
         write(*,*) 'cell ',cellID(cell)
         write(*,*) '   v0      ',v0
         write(*,*) '   B       ',B
         write(*,*) '   Binv    ',Binv
         write(*,*) '   detB    ',detB
         !call DMPlexComputeCellGeometryFVM(dm,cellID(cell),vol,centroid,normal,ierr)
         !write(*,*) '   vol      ', vol
         !write(*,*) '   centroid ', centroid
         !write(*,*) '   normal   ', normal
      end do
      call PetscPrintf(PETSC_COMM_SELF,"\n",ierr)
      call ISRestoreIndicesF90(cellIS,cellID,ierr)
      call ISDestroy(cellIS,ierr)
   end do
   call ISRestoreIndicesF90(CSIS,setID,ierr)
   call ISDestroy(CSIS,ierr)
   deallocate(detB)
   deAllocate(Binv)
   deAllocate(B)
   deAllocate(v0)
   call PetscQuadratureDestroy(q,ierr)

   allocate(centroid(dim))
   allocate(normal(dim))
   call DMGetLabel(dm,"Face Sets",FSLabel,ierr)
   call DMLabelGetValueIS(FSLabel,FSIS,ierr)
   call ISGetIndicesF90(FSIS,setID,ierr)
   do set = 1, size(setID)
      write(IOBuffer,*) 'Face Set', setID(set),'\n'
      call PetscPrintf(PETSC_COMM_SELF,IOBuffer,ierr)
      call DMLabelGetStratumIS(FSLabel,setID(set),cellIS,ierr)
      call ISGetIndicesF90(cellIS,cellID,ierr)
      do cell = 1, size(cellID)
         !call DMPlexComputeCellGeometryAffineFEM(dm,cellID(cell),v0,B,Binv,detB,ierr)
         !write(*,*) 'cell ',cellID(cell)
         !write(*,*) '   v0   ',v0
         !write(*,*) '   B    ',B
         !write(*,*) '   Binv ',Binv
         !write(*,*) '   detB ',detB
         call DMPlexComputeCellGeometryFVM(dm,cellID(cell),vol,centroid,normal,ierr)
         write(*,*) '   vol         ', vol
         write(*,*) '   centroid    ', centroid
         write(*,*) '   normal      ', normal
      end do
      call ISRestoreIndicesF90(cellIS,cellID,ierr)
      call ISDestroy(cellIS,ierr)
   end do
   call ISRestoreIndicesF90(FSIS,setID,ierr)
   call ISDestroy(FSIS,ierr)

   Call MEF90CtxDestroy(MEF90Ctx,ierr)   
   Call MEF90Finalize(ierr)
   Call PetscFinalize()
End Program  TestDMPlexComputeGeometry
