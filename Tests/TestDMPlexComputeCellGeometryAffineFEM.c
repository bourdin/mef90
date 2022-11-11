static char help[] = "Trying to understand DMPlexComputeCellGeometryAffineFEM";
#include <petsc.h>
#include "petsc/private/petscfeimpl.h"

int main(int argc,char **argv)
{
  DM               dm;
  char             filename[2048];
  Vec              coordLoc;
  PetscReal        v0[3];
  PetscReal        J[9];
  PetscReal        invJ[9];
  PetscReal        detJ;
  PetscInt         dim;

  const PetscReal  X0[3] = {0.,0.,0.};
  const PetscReal  X1[3] = {1.,0.,0.};
  const PetscReal  X2[3] = {0.,1.,0.};
  const PetscReal  X3[3] = {0.,0.,1.};

  const PetscReal  Y0[3] = {-1., 1.,-1.};
  const PetscReal  Y1[3] = {-1.,-1.,-1.};
  const PetscReal  Y2[3] = { 1.,-1.,-1.};
  const PetscReal  Y3[3] = {-1.,-1., 1.};

  PetscReal        X[3],Y[3];
  

  PetscCall(PetscInitialize(&argc,&argv,NULL,help));
  PetscOptionsBegin(PETSC_COMM_WORLD,"","TestDMPlexComputeCellGeometryAffineFEM options","none");
  PetscCall(PetscOptionsString("-i","filename to read","",filename,filename,sizeof(filename),NULL));
  PetscOptionsEnd();

  PetscCall(PetscPrintf(PETSC_COMM_WORLD," filename %s\n",filename));
  PetscCall(DMPlexCreateFromFile(PETSC_COMM_WORLD, filename, NULL, PETSC_FALSE, &dm));

  PetscCall(DMGetCoordinatesLocal(dm,&coordLoc));
  PetscCall(VecView(coordLoc,PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(DMGetDimension(dm,&dim));
  PetscCall(DMPlexComputeCellGeometryAffineFEM(dm,0,v0,J,invJ,&detJ));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"v0\n"));
  PetscCall(PetscRealView(dim,v0,PETSC_VIEWER_STDOUT_WORLD));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"J\n"));
  for (int i = 0; i < dim; i++){
    PetscCall(PetscRealView(dim,&J[i*dim],PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"invJ\n"));
  for (int i = 0; i < dim; i++){
    PetscCall(PetscRealView(dim,&invJ[i*dim],PETSC_VIEWER_STDOUT_WORLD));
  }
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"detJ : %g\n",detJ));

  PetscCall(DMPlexCoordinatesToReference(dm,0,1,X0,X));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"X0\n"));
  PetscCall(PetscRealView(dim,X,PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(DMPlexCoordinatesToReference(dm,0,1,X1,X));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"X1\n"));
  PetscCall(PetscRealView(dim,X,PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(DMPlexCoordinatesToReference(dm,0,1,X2,X));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"X2\n"));
  PetscCall(PetscRealView(dim,X,PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(DMPlexCoordinatesToReference(dm,0,1,X3,X));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"X3\n"));
  PetscCall(PetscRealView(dim,X,PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(DMPlexReferenceToCoordinates(dm,0,1,Y0,Y));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Y0\n"));
  PetscCall(PetscRealView(dim,Y,PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(DMPlexReferenceToCoordinates(dm,0,1,Y1,Y));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Y1\n"));
  PetscCall(PetscRealView(dim,Y,PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(DMPlexReferenceToCoordinates(dm,0,1,Y2,Y));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Y2\n"));
  PetscCall(PetscRealView(dim,Y,PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(DMPlexReferenceToCoordinates(dm,0,1,Y3,Y));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,"Y3\n"));
  PetscCall(PetscRealView(dim,Y,PETSC_VIEWER_STDOUT_WORLD));


  // CoordinatesRefToReal(dim,dim,v0,X0,J,X1,X);
  // PetscCall(PetscPrintf(PETSC_COMM_WORLD,"X1\n"));
  // PetscCall(PetscRealView(dim,X,PETSC_VIEWER_STDOUT_WORLD));

  // CoordinatesRefToReal(dim,dim,v0,X0,J,X2,X);
  // PetscCall(PetscPrintf(PETSC_COMM_WORLD,"X2\n"));
  // PetscCall(PetscRealView(dim,X,PETSC_VIEWER_STDOUT_WORLD));

  // CoordinatesRefToReal(dim,dim,v0,X0,J,X3,X);
  // PetscCall(PetscPrintf(PETSC_COMM_WORLD,"X3\n"));
  // PetscCall(PetscRealView(dim,X,PETSC_VIEWER_STDOUT_WORLD));

  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}
