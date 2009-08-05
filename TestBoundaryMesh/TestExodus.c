static char help[] = "Test distribution of properties using a mesh\n\n";

#include <stdio.h>
#include <stdlib.h>
#include "petsc.h"
#include <netcdf.h>
#include <exodusII.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main (int argc, char ** argv) {
  int   exoid;
  int   CPU_word_size = 0, IO_word_size = 0;
  //    const PetscMPIInt rank = mesh->commRank();
  float version;
  char  title[MAX_LINE_LENGTH+1], elem_type[MAX_STR_LENGTH+1];
  int   num_dim, num_nodes, num_elem, num_elem_blk, num_node_sets, num_side_sets;
  int   ierr;
  char filename[PETSC_MAX_PATH_LEN+1];

  ierr = PetscInitialize(&argc,&argv,(char *)0,help);  
  
  ierr = PetscOptionsGetString(PETSC_NULL, "-f", filename, PETSC_MAX_PATH_LEN, PETSC_NULL);CHKERRQ(ierr);

  // Read database parameters
  exoid = ex_open(filename, EX_READ, &CPU_word_size, &IO_word_size, &version);CHKERRQ(!exoid);
  ierr = ex_get_init(exoid, title, &num_dim, &num_nodes, &num_elem, &num_elem_blk, &num_node_sets, &num_side_sets);CHKERRQ(ierr);
  
  int   num_eb_props, num_ss_props, num_ns_props;
  int   **eb_prop_values, **ss_prop_values, **ns_prop_values;
  char  **eb_prop_names, **ss_prop_names, **ns_prop_names;
  float fdum;
  char  cdum;
  

  // Element Blocks
  ierr=ex_inquire(exoid, EX_INQ_EB_PROP, &num_eb_props, &fdum, &cdum);
  ierr=PetscPrintf(PETSC_COMM_WORLD, "num_eb_props: %i\n", num_eb_props);CHKERRQ(ierr);
  ierr=PetscMalloc2(num_eb_props,int*,&eb_prop_values,num_eb_props,char*,&eb_prop_names);CHKERRQ(ierr);
  for(int p=0; p<num_eb_props; ++p) {
    ierr=PetscMalloc((MAX_STR_LENGTH+1)* sizeof(char), &eb_prop_names[p]);CHKERRQ(ierr);
    ierr=PetscMalloc(num_elem_blk*sizeof(int), &eb_prop_values[p]);CHKERRQ(ierr);
  }
  ierr=ex_get_prop_names(exoid,EX_ELEM_BLOCK,eb_prop_names);
  for(int p=1; p<num_eb_props; ++p) {
    ierr=PetscPrintf(PETSC_COMM_WORLD, "eb_prop_name[%i]: %s\n", p, eb_prop_names[p]);CHKERRQ(ierr);
    ierr=ex_get_prop_array(exoid, EX_ELEM_BLOCK, eb_prop_names[p], eb_prop_values[p]);
    ierr=PetscPrintf(PETSC_COMM_WORLD, "eb_prop_values[%i]: ", p);CHKERRQ(ierr);
    for(int i=0; i<num_elem_blk; ++i){
      ierr=PetscPrintf(PETSC_COMM_WORLD, "%i, ", eb_prop_values[p][i]);CHKERRQ(ierr);
    }
    ierr=PetscPrintf(PETSC_COMM_WORLD, "\n");CHKERRQ(ierr);
  }
  ierr=PetscPrintf(PETSC_COMM_WORLD, "\n\n\n");CHKERRQ(ierr);

  // Side Sets
  ierr=ex_inquire(exoid, EX_INQ_SS_PROP, &num_ss_props, &fdum, &cdum);
  ierr=PetscPrintf(PETSC_COMM_WORLD, "num_ss_props: %i\n", num_ss_props);CHKERRQ(ierr);
  ierr=PetscMalloc2(num_ss_props,int*,&ss_prop_values,num_ss_props,char*,&ss_prop_names);CHKERRQ(ierr);
  for(int p=0; p<num_ss_props; ++p) {
    ierr=PetscMalloc((MAX_STR_LENGTH+1)* sizeof(char), &ss_prop_names[p]);CHKERRQ(ierr);
    ierr=PetscMalloc(num_elem_blk*sizeof(int), &ss_prop_values[p]);CHKERRQ(ierr);
  }
  ierr=ex_get_prop_names(exoid,EX_SIDE_SET,ss_prop_names);
  for(int p=1; p<num_ss_props; ++p) {
    ierr=PetscPrintf(PETSC_COMM_WORLD, "ss_prop_name[%i]: %s\n", p, ss_prop_names[p]);CHKERRQ(ierr);
    ierr=ex_get_prop_array(exoid, EX_SIDE_SET, ss_prop_names[p], ss_prop_values[p]);
    ierr=PetscPrintf(PETSC_COMM_WORLD, "ss_prop_values[%i]: ", p);CHKERRQ(ierr);
    for(int i=0; i<num_side_sets; ++i){
      ierr=PetscPrintf(PETSC_COMM_WORLD, "%i, ", ss_prop_values[p][i]);CHKERRQ(ierr);
    }
    ierr=PetscPrintf(PETSC_COMM_WORLD, "\n");CHKERRQ(ierr);
  }
  ierr=PetscPrintf(PETSC_COMM_WORLD, "\n\n\n");CHKERRQ(ierr);
  
  //Node Sets
  ierr=ex_inquire(exoid, EX_INQ_NS_PROP, &num_ns_props, &fdum, &cdum);
  ierr=PetscPrintf(PETSC_COMM_WORLD, "num_ns_props: %i\n", num_ns_props);CHKERRQ(ierr);
  ierr=PetscMalloc2(num_ns_props,int*,&ns_prop_values,num_ns_props,char*,&ns_prop_names);CHKERRQ(ierr);
  for(int p=0; p<num_ns_props; ++p) {
    ierr=PetscMalloc((MAX_STR_LENGTH+1)* sizeof(char), &ns_prop_names[p]);CHKERRQ(ierr);
    ierr=PetscMalloc(num_elem_blk*sizeof(int), &ns_prop_values[p]);CHKERRQ(ierr);
  }
  ierr=ex_get_prop_names(exoid,EX_NODE_SET,ns_prop_names);
  for(int p=1; p<num_ns_props; ++p) {
    ierr=PetscPrintf(PETSC_COMM_WORLD, "ns_prop_name[%i]: %s\n", p, ns_prop_names[p]);CHKERRQ(ierr);
    ierr=ex_get_prop_array(exoid, EX_NODE_SET, ns_prop_names[p], ns_prop_values[p]);
    ierr=PetscPrintf(PETSC_COMM_WORLD, "ns_prop_values[%i]: ", p);CHKERRQ(ierr);
    for(int i=0; i<num_node_sets; ++i){
      ierr=PetscPrintf(PETSC_COMM_WORLD, "%i, ", ns_prop_values[p][i]);CHKERRQ(ierr);
    }
    ierr=PetscPrintf(PETSC_COMM_WORLD, "\n");CHKERRQ(ierr);
  }


  PetscFinalize();
  return 0;
}