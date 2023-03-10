static char help[] = "Trying to understand DMPlex labels\n\n";

#include <petsc.h>

#undef __FUNCT__
#define __FUNCT__ "dmViewSets"
int dmViewSets(DM dm,PetscViewer viewer)
{
  PetscInt       pStart,pEnd,cStart,cEnd,vStart,vEnd,eStart,eEnd;
  IS             setIS,pointsIS;
  const PetscInt *setIDs;
  PetscInt       numSets;
  PetscInt       set,dim;

  PetscCall(PetscViewerASCIIPrintf(viewer,"=== SETS ===\n"));
  PetscCall(DMGetDimension(dm,&dim));
  PetscCall(DMPlexGetChart(dm,&pStart,&pEnd));
  PetscCall(DMPlexGetHeightStratum(dm,0,&cStart,&cEnd));
  PetscCall(DMPlexGetHeightStratum(dm,1,&eStart,&eEnd));
  PetscCall(DMPlexGetHeightStratum(dm,dim+1,&vStart,&vEnd));
  PetscCall(PetscViewerASCIIPrintf(viewer,"pStart: %" PetscInt_FMT ", pEnd: %" PetscInt_FMT "\n",pStart,pEnd));
  PetscCall(PetscViewerASCIIPrintf(viewer,"cStart: %" PetscInt_FMT ", cEnd: %" PetscInt_FMT "\n",cStart,cEnd));
  PetscCall(PetscViewerASCIIPrintf(viewer,"vStart: %" PetscInt_FMT ", vEnd: %" PetscInt_FMT "\n",vStart,vEnd));
  PetscCall(PetscViewerASCIIPrintf(viewer,"eStart: %" PetscInt_FMT ", eEnd: %" PetscInt_FMT "\n",eStart,eEnd));

  PetscCall(DMGetLabelSize(dm,"Cell Sets",&numSets));
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n%" PetscInt_FMT " Cell sets: \n",numSets));
  if (numSets > 0) {
    PetscCall(DMGetLabelIdIS(dm,"Cell Sets",&setIS));
    PetscCall(ISView(setIS,viewer));
    // PetscCall(ISGetIndices(setIS,&setIDs));
    // for (set = 0; set < numSets; set++){
    //   PetscCall(PetscViewerASCIIPrintf(viewer,"=== Cell Set %" PetscInt_FMT " ===\n",setIDs[set]));
    //   PetscCall(DMGetStratumIS(dm,"Cell Sets",setIDs[set],&pointsIS));
    //   PetscCall(ISView(pointsIS,viewer));
    //   PetscCall(ISDestroy(&pointsIS));
    // }
    // PetscCall(ISRestoreIndices(setIS,&setIDs));
    PetscCall(ISDestroy(&setIS));
  }
  PetscCall(DMGetLabelSize(dm,"Face Sets",&numSets));
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n%" PetscInt_FMT " Face sets: \n",numSets));
  if (numSets > 0) {
    PetscCall(DMGetLabelIdIS(dm,"Face Sets",&setIS));
    PetscCall(ISView(setIS,viewer));
    // PetscCall(ISGetIndices(setIS,&setIDs));
    // for (set = 0; set < numSets; set++){
    //   PetscCall(PetscViewerASCIIPrintf(viewer,"=== Face Set %" PetscInt_FMT " ===\n",setIDs[set]));
    //   PetscCall(DMGetStratumIS(dm,"Face Sets",setIDs[set],&pointsIS));
    //   PetscCall(ISView(pointsIS,viewer));
    //   PetscCall(ISDestroy(&pointsIS));
    // }
    // PetscCall(ISRestoreIndices(setIS,&setIDs));
    PetscCall(ISDestroy(&setIS));
  }

  PetscCall(DMGetLabelSize(dm,"Vertex Sets",&numSets));
  PetscCall(PetscViewerASCIIPrintf(viewer,"\n%" PetscInt_FMT " Vertex sets: \n",numSets));
  if (numSets > 0) {
    PetscCall(DMGetLabelIdIS(dm,"Vertex Sets",&setIS));
    PetscCall(ISView(setIS,viewer));
    // PetscCall(ISGetIndices(setIS,&setIDs));
    // for (set = 0; set < numSets; set++){
    //   PetscCall(PetscViewerASCIIPrintf(viewer,"=== Vertex Set %" PetscInt_FMT " ===\n",setIDs[set]));
    //   PetscCall(DMGetStratumIS(dm,"Vertex Sets",setIDs[set],&pointsIS));
    //   PetscCall(ISView(pointsIS,viewer));
    //   PetscCall(ISDestroy(&pointsIS));
    // }
    // PetscCall(ISRestoreIndices(setIS,&setIDs));
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

  PetscCall(PetscInitialize(&argc,&argv,NULL,help));
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
  PetscOptionsBegin(PETSC_COMM_WORLD,"","viewSets options","none");
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
  if (flg) PetscCall(PetscViewerDestroy(&viewer));

  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}
