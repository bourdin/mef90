static char help[] = "Trying to understand DMPlex and sections - part I\n\n";

#include <petsc.h>

#undef __FUNCT__
#define __FUNCT__ "dmViewDAG"
int dmViewDAG(DM dm,PetscViewer viewer)
{
  PetscInt       pStart,pEnd,p1,p2;
  PetscInt       vStart,vEnd,v,c,nc;
  PetscInt       *transitiveClosure = NULL;
  PetscInt       numPoints;
  const PetscInt *cone;
  const PetscInt *support;
  PetscInt       dim,depth,height;
  PetscSection   coordSection;
  Vec            coord;
  PetscReal      *xyz=NULL;

  
  PetscCall(PetscViewerASCIIPrintf(viewer,"=== DAG ===\n"));
  PetscCall(DMPlexGetChart(dm,&pStart,&pEnd));
  PetscCall(DMGetDimension(dm,&dim));
  PetscCall(PetscViewerASCIIPrintf(viewer,"\nDimension %" PetscInt_FMT ":\n",dim));
  PetscCall(DMGetCoordinateSection(dm,&coordSection));
  PetscCall(DMGetCoordinatesLocal(dm,&coord));
  PetscCall(DMPlexGetDepthStratum(dm,0,&vStart,&vEnd));
  PetscCall(PetscViewerASCIIPrintf(viewer,"\nvertex range in DAG %" PetscInt_FMT "-%" PetscInt_FMT ":\n",vStart,vEnd));
  for (v = vStart; v < vEnd; v++) {
    PetscCall(DMPlexVecGetClosure(dm,coordSection,coord,v,&nc,&xyz));
    PetscCall(PetscViewerASCIIPrintf(viewer,"   vertex %" PetscInt_FMT ":",v));
    for (c = 0; c < nc; c++) {
      PetscCall(PetscViewerASCIIPrintf(viewer," %g",xyz[c]));
    }
    PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
  }
  PetscCall(DMPlexVecRestoreClosure(dm,coordSection,coord,v,&nc,&xyz));
  for (p1 = pStart; p1 < pEnd; p1++) {
    PetscCall(PetscViewerASCIIPrintf(viewer,"\npoint %" PetscInt_FMT ":\n",p1));
    PetscCall(DMPlexGetPointDepth(dm,p1,&depth));
    PetscCall(DMPlexGetPointHeight(dm,p1,&height));
    PetscCall(PetscViewerASCIIPrintf(viewer,"   Depth: %" PetscInt_FMT " Height: %" PetscInt_FMT "\n",depth,height));
    PetscCall(DMPlexGetConeSize(dm,p1,&numPoints));
    PetscCall(DMPlexGetCone(dm,p1,&cone));
    PetscCall(PetscViewerASCIIPrintf(viewer,"   Cone: "));
    for (p2 = 0; p2 < numPoints; p2++) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT ",",cone[p2]));
    }
    PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
    
    PetscCall(DMPlexGetSupportSize(dm,p1,&numPoints));
    PetscCall(DMPlexGetSupport(dm,p1,&support));
    PetscCall(PetscViewerASCIIPrintf(viewer,"   Support: "));
    for (p2 = 0; p2 < numPoints; p2++) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT ",",support[p2]));
    }
    PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
    
    PetscCall(DMPlexGetTransitiveClosure(dm,p1,PETSC_TRUE,&numPoints,&transitiveClosure));
    PetscCall(PetscViewerASCIIPrintf(viewer,"   TransitiveClosure IN: "));
    for (p2 = 0; p2 < numPoints; p2++) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT ",",transitiveClosure[p2*2]));
    }
    PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
    PetscCall(DMPlexRestoreTransitiveClosure(dm,p1,PETSC_TRUE,&numPoints,&transitiveClosure));
    
    PetscCall(DMPlexGetTransitiveClosure(dm,p1,PETSC_FALSE,&numPoints,&transitiveClosure));
    PetscCall(PetscViewerASCIIPrintf(viewer,"   TransitiveClosure OUT: "));
    for (p2 = 0; p2 < numPoints; p2++) {
      PetscCall(PetscViewerASCIIPrintf(viewer,"%" PetscInt_FMT ",",transitiveClosure[p2*2]));
    }
    PetscCall(PetscViewerASCIIPrintf(viewer,"\n"));
    PetscCall(DMPlexRestoreTransitiveClosure(dm,p1,PETSC_TRUE,&numPoints,&transitiveClosure));
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

  PetscCall(PetscViewerASCIIPrintf(viewer,"=== SETS ===\n"));
  PetscCall(DMPlexGetChart(dm,&pStart,&pEnd));
  PetscCall(DMPlexGetHeightStratum(dm,0,&cStart,&cEnd));
  PetscCall(DMPlexGetHeightStratum(dm,1,&eStart,&eEnd));
  PetscCall(DMPlexGetDepthStratum(dm,0,&vStart,&vEnd));
  PetscCall(PetscViewerASCIIPrintf(viewer,"pStart: %" PetscInt_FMT ", pEnd: %" PetscInt_FMT "\n",pStart,pEnd));
  PetscCall(PetscViewerASCIIPrintf(viewer,"cStart: %" PetscInt_FMT ", cEnd: %" PetscInt_FMT "\n",cStart,cEnd));
  PetscCall(PetscViewerASCIIPrintf(viewer,"vStart: %" PetscInt_FMT ", vEnd: %" PetscInt_FMT "\n",vStart,vEnd));
  PetscCall(PetscViewerASCIIPrintf(viewer,"eStart: %" PetscInt_FMT ", eEnd: %" PetscInt_FMT "\n",eStart,eEnd));

  PetscCall(DMGetLabelSize(dm,"Cell Sets",&numSets));
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n%" PetscInt_FMT " Cell sets: \n",numSets));
  if (numSets > 0) {
    PetscCall(DMGetLabelIdIS(dm,"Cell Sets",&setIS));
    PetscCall(ISView(setIS,viewer));
    PetscCall(ISGetIndices(setIS,&setIDs));
    for (set = 0; set < numSets; set++){
      PetscCall(PetscViewerASCIIPrintf(viewer,"=== Set %" PetscInt_FMT " ===\n",setIDs[set]));
      PetscCall(DMGetStratumIS(dm,"Cell Sets",setIDs[set],&pointsIS));
      PetscCall(ISView(pointsIS,viewer));
      PetscCall(ISDestroy(&pointsIS));
    }
    PetscCall(ISRestoreIndices(setIS,&setIDs));
    PetscCall(ISDestroy(&setIS));
  }
  PetscCall(DMGetLabelSize(dm,"Face Sets",&numSets));
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n%" PetscInt_FMT " Face sets: \n",numSets));
  if (numSets > 0) {
    PetscCall(DMGetLabelIdIS(dm,"Face Sets",&setIS));
    PetscCall(ISView(setIS,viewer));
    PetscCall(ISGetIndices(setIS,&setIDs));
    for (set = 0; set < numSets; set++){
      PetscCall(PetscViewerASCIIPrintf(viewer,"=== Set %" PetscInt_FMT " ===\n",setIDs[set]));
      PetscCall(DMGetStratumIS(dm,"Face Sets",setIDs[set],&pointsIS));
      PetscCall(ISView(pointsIS,viewer));
      PetscCall(ISDestroy(&pointsIS));
    }
    PetscCall(ISRestoreIndices(setIS,&setIDs));
    PetscCall(ISDestroy(&setIS));
  }

  PetscCall(DMGetLabelSize(dm,"Vertex Sets",&numSets));
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n%" PetscInt_FMT " Vertex sets: \n",numSets));
  if (numSets > 0) {
    PetscCall(DMGetLabelIdIS(dm,"Vertex Sets",&setIS));
    PetscCall(ISView(setIS,viewer));
    PetscCall(ISGetIndices(setIS,&setIDs));
    for (set = 0; set < numSets; set++){
      PetscCall(PetscViewerASCIIPrintf(viewer,"=== Set %" PetscInt_FMT " ===\n",setIDs[set]));
      PetscCall(DMGetStratumIS(dm,"Vertex Sets",setIDs[set],&pointsIS));
      PetscCall(ISView(pointsIS,viewer));
      PetscCall(ISDestroy(&pointsIS));
    }
    PetscCall(ISRestoreIndices(setIS,&setIDs));
    PetscCall(ISDestroy(&setIS));
  }
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  DM             dm,dmDist;
  char           filename[2048],outfilename[2048];
  PetscBool      interpolate=PETSC_FALSE,flg=PETSC_FALSE;
  PetscMPIInt    rank;
  PetscViewer    viewer;
  // PetscSection   sectionCoord; 
  // Vec            coordLoc;

  PetscCall(PetscInitialize(&argc,&argv,NULL,help));
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscOptionsBegin(PETSC_COMM_WORLD,"","viewDAG options","none");
  PetscCall(PetscOptionsString("-i","filename to read","",filename,filename,sizeof(filename),NULL));
  PetscCall(PetscOptionsString("-o","filename to write","",outfilename,outfilename,sizeof(outfilename),&flg));
  PetscCall(PetscOptionsBool("-interpolate","Generate intermediate mesh elements","",interpolate,&interpolate,NULL));
  PetscOptionsEnd();

  PetscCall(PetscPrintf(PETSC_COMM_WORLD," filename %s\n",filename));
  PetscCall(DMPlexCreateFromFile(PETSC_COMM_WORLD, filename, NULL, interpolate, &dm));
  PetscCall(DMPlexDistribute(dm,0,NULL,&dmDist));
  if (dmDist) {
    DMDestroy(&dm);
    dm   = dmDist;
  }
  //PetscCall(DMSetFromOptions(dm));
  viewer = PETSC_VIEWER_STDOUT_WORLD;
  if (flg) {
    PetscCall(PetscSNPrintf(outfilename,sizeof(outfilename),outfilename,rank));
    PetscCall(PetscViewerASCIIOpen(PETSC_COMM_SELF,outfilename,&viewer));
  }
  PetscCall(dmViewSets(dm,viewer));
  PetscCall(dmViewDAG(dm,viewer));

  // PetscCall(DMGetCoordinateSection(dm,&sectionCoord));
  // PetscCall(DMGetCoordinatesLocal(dm,&coordLoc));
  // PetscCall(VecView(coordLoc,PETSC_VIEWER_STDOUT_WORLD));

  if (flg) PetscCall(PetscViewerDestroy(&viewer));
  PetscCall(PetscFinalize());
  return 0;
}
