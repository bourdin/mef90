static char help[] = "Trying to understand DMPlex and sections - part I\n\n";

#include <petsc.h>

#undef __FUNCT__
#define __FUNCT__ "dmViewDAG"
int dmViewDAG(DM dm,PetscViewer viewer)
{
  PetscInt       pStart,pEnd,p1,p2;
  PetscInt       vStart,vEnd,v,c,nc;
  PetscErrorCode ierr;
  //const PetscInt *cone,*support;
  PetscInt       *transitiveClosure = NULL;
  PetscInt       numPoints;
  const PetscInt *cone;
  const PetscInt *support;
  PetscInt       dim,depth,height;
  PetscSection   coordSection;
  Vec            coord;
  PetscReal      *xyz=NULL;

  
  ierr = PetscViewerASCIIPrintf(viewer,"=== DAG ===\n");CHKERRQ(ierr);
  ierr = DMPlexGetChart(dm,&pStart,&pEnd);
  ierr = DMGetDimension(dm,&dim);
  ierr = PetscViewerASCIIPrintf(viewer,"\nDimension %d:\n",dim);CHKERRQ(ierr);
  ierr = DMGetCoordinateSection(dm,&coordSection);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm,&coord);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm,0,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\nvertex range in DAG %d-%d:\n",vStart,vEnd);CHKERRQ(ierr);
  for (v = vStart; v < vEnd; v++) {
    ierr = DMPlexVecGetClosure(dm,coordSection,coord,v,&nc,&xyz);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"   vertex %d:",v);CHKERRQ(ierr);
    for (c = 0; c < nc; c++) {
      ierr = PetscViewerASCIIPrintf(viewer," %g",xyz[c]);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
  }
  ierr = DMPlexVecRestoreClosure(dm,coordSection,coord,v,&nc,&xyz);CHKERRQ(ierr);
  for (p1 = pStart; p1 < pEnd; p1++) {
    ierr = PetscViewerASCIIPrintf(viewer,"\npoint %d:\n",p1);CHKERRQ(ierr);
    PetscCall(DMPlexGetPointDepth(dm,p1,&depth));
    PetscCall(DMPlexGetPointHeight(dm,p1,&height));
    PetscCall(PetscViewerASCIIPrintf(viewer,"   Depth: %d Height: %d\n",depth,height));
    ierr = DMPlexGetConeSize(dm,p1,&numPoints);CHKERRQ(ierr);
    ierr = DMPlexGetCone(dm,p1,&cone);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"   Cone: ");CHKERRQ(ierr);
    for (p2 = 0; p2 < numPoints; p2++) {
      ierr = PetscViewerASCIIPrintf(viewer,"%d,",cone[p2]);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
    
    ierr = DMPlexGetSupportSize(dm,p1,&numPoints);CHKERRQ(ierr);
    ierr = DMPlexGetSupport(dm,p1,&support);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"   Support: ");CHKERRQ(ierr);
    for (p2 = 0; p2 < numPoints; p2++) {
      ierr = PetscViewerASCIIPrintf(viewer,"%d,",support[p2]);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
    
    ierr = DMPlexGetTransitiveClosure(dm,p1,PETSC_TRUE,&numPoints,&transitiveClosure);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"   TransitiveClosure IN: ");CHKERRQ(ierr);
    for (p2 = 0; p2 < numPoints; p2++) {
      ierr = PetscViewerASCIIPrintf(viewer,"%d,",transitiveClosure[p2*2]);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
    ierr = DMPlexRestoreTransitiveClosure(dm,p1,PETSC_TRUE,&numPoints,&transitiveClosure);CHKERRQ(ierr);
    
    ierr = DMPlexGetTransitiveClosure(dm,p1,PETSC_FALSE,&numPoints,&transitiveClosure);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"   TransitiveClosure OUT: ");CHKERRQ(ierr);
    for (p2 = 0; p2 < numPoints; p2++) {
      ierr = PetscViewerASCIIPrintf(viewer,"%d,",transitiveClosure[p2*2]);CHKERRQ(ierr);
    }
    ierr = PetscViewerASCIIPrintf(viewer,"\n");CHKERRQ(ierr);
    ierr = DMPlexRestoreTransitiveClosure(dm,p1,PETSC_TRUE,&numPoints,&transitiveClosure);CHKERRQ(ierr);
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "dmViewSets"
int dmViewSets(DM dm,PetscViewer viewer)
{
  PetscInt       pStart,pEnd,cStart,cEnd,vStart,vEnd,eStart,eEnd;
  IS             setIS,pointsIS;
  const PetscInt *setIDs;
  PetscInt       numSets;
  PetscInt       set;
  PetscErrorCode ierr;

  ierr = PetscViewerASCIIPrintf(viewer,"=== SETS ===\n");CHKERRQ(ierr);
  ierr = DMPlexGetChart(dm,&pStart,&pEnd);
  ierr = DMPlexGetHeightStratum(dm,0,&cStart,&cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm,1,&eStart,&eEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm,0,&vStart,&vEnd);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"pStart: %d, pEnd: %d\n",pStart,pEnd);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"cStart: %d, cEnd: %d\n",cStart,cEnd);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"vStart: %d, vEnd: %d\n",vStart,vEnd);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"eStart: %d, eEnd: %d\n",eStart,eEnd);CHKERRQ(ierr);

  ierr = DMGetLabelSize(dm,"Cell Sets",&numSets);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\n%d Cell sets: \n",numSets);CHKERRQ(ierr);
  if (numSets > 0) {
    ierr = DMGetLabelIdIS(dm,"Cell Sets",&setIS);CHKERRQ(ierr);
    ierr = ISView(setIS,viewer);CHKERRQ(ierr);
    ierr = ISGetIndices(setIS,&setIDs);CHKERRQ(ierr);
    for (set = 0; set < numSets; set++){
      ierr = PetscViewerASCIIPrintf(viewer,"=== Set %d ===\n",setIDs[set]);CHKERRQ(ierr);
      ierr = DMGetStratumIS(dm,"Cell Sets",setIDs[set],&pointsIS);CHKERRQ(ierr);
      ierr = ISView(pointsIS,viewer);CHKERRQ(ierr);
      ierr = ISDestroy(&pointsIS);CHKERRQ(ierr);
    }
    ierr = ISRestoreIndices(setIS,&setIDs);CHKERRQ(ierr);
    ierr = ISDestroy(&setIS);CHKERRQ(ierr);
  }
  ierr = DMGetLabelSize(dm,"Face Sets",&numSets);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\n%d Face sets: \n",numSets);CHKERRQ(ierr);
  if (numSets > 0) {
    ierr = DMGetLabelIdIS(dm,"Face Sets",&setIS);CHKERRQ(ierr);
    ierr = ISView(setIS,viewer);CHKERRQ(ierr);
    ierr = ISGetIndices(setIS,&setIDs);CHKERRQ(ierr);
    for (set = 0; set < numSets; set++){
      ierr = PetscViewerASCIIPrintf(viewer,"=== Set %d ===\n",setIDs[set]);CHKERRQ(ierr);
      ierr = DMGetStratumIS(dm,"Face Sets",setIDs[set],&pointsIS);CHKERRQ(ierr);
      ierr = ISView(pointsIS,viewer);CHKERRQ(ierr);
      ierr = ISDestroy(&pointsIS);CHKERRQ(ierr);
    }
    ierr = ISRestoreIndices(setIS,&setIDs);CHKERRQ(ierr);
    ierr = ISDestroy(&setIS);CHKERRQ(ierr);
  }

  ierr = DMGetLabelSize(dm,"Vertex Sets",&numSets);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPrintf(viewer,"\n%d Vertex sets: \n",numSets);CHKERRQ(ierr);
  if (numSets > 0) {
    ierr = DMGetLabelIdIS(dm,"Vertex Sets",&setIS);CHKERRQ(ierr);
    ierr = ISView(setIS,viewer);CHKERRQ(ierr);
    ierr = ISGetIndices(setIS,&setIDs);CHKERRQ(ierr);
    for (set = 0; set < numSets; set++){
      ierr = PetscViewerASCIIPrintf(viewer,"=== Set %d ===\n",setIDs[set]);CHKERRQ(ierr);
      ierr = DMGetStratumIS(dm,"Vertex Sets",setIDs[set],&pointsIS);CHKERRQ(ierr);
      ierr = ISView(pointsIS,viewer);CHKERRQ(ierr);
      ierr = ISDestroy(&pointsIS);CHKERRQ(ierr);
    }
    ierr = ISRestoreIndices(setIS,&setIDs);CHKERRQ(ierr);
    ierr = ISDestroy(&setIS);CHKERRQ(ierr);
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  DM             dm,dmDist;
  PetscErrorCode ierr;
  char           filename[2048],outfilename[2048];
  PetscBool      interpolate=PETSC_FALSE;
  PetscMPIInt    rank;
  PetscViewer    viewer;
  PetscSection   sectionCoord; 
  Vec            coordLoc;
  PetscSF        migrationSF;

  ierr = PetscInitialize(&argc,&argv,NULL,help);CHKERRQ(ierr);
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscOptionsBegin(PETSC_COMM_WORLD,"","viewDAG options","none");
  PetscCall(PetscOptionsString("-i","filename to read","",filename,filename,sizeof(filename),NULL));
  PetscCall(PetscOptionsString("-o","filename to write","",outfilename,outfilename,sizeof(outfilename),NULL));
  PetscCall(PetscOptionsBool("-interpolate","Generate intermediate mesh elements","",interpolate,&interpolate,NULL));
  PetscOptionsEnd();

  //ierr = DMPlexCreateFromFile(PETSC_COMM_WORLD,filename,interpolate,&dm);CHKERRQ(ierr);
  PetscCall(DMPlexCreateFromFile(PETSC_COMM_WORLD, filename, NULL, interpolate, &dm));
  ierr = DMPlexDistribute(dm,0,&migrationSF,&dmDist);CHKERRQ(ierr);
  if (dmDist) {
    DMDestroy(&dm);
    dm   = dmDist;
  }
  ierr = DMSetFromOptions(dm);CHKERRQ(ierr);
  ierr = PetscSNPrintf(outfilename,sizeof(outfilename),outfilename,rank);CHKERRQ(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,outfilename,&viewer);CHKERRQ(ierr);
  ierr = dmViewSets(dm,viewer);CHKERRQ(ierr);
  ierr = dmViewDAG(dm,viewer);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  ierr = DMGetCoordinateSection(dm,&sectionCoord);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm,&coordLoc);CHKERRQ(ierr);
  ierr = VecView(coordLoc,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return 0;
}
