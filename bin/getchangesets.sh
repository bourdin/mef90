#!/bin/sh
# Q&D to retrieve the changesets of the sensitive packages. 
echo "PETSc-dev changeset:" ` cd $PETSC_DIR && hg log -l 1 | grep changeset | awk '{print $2}'`
echo "	Builder changeset:" ` cd $PETSC_DIR/config/BuildSystem/ && hg log -l 1 | grep changeset | awk '{print $2}'`
echo "MEF90 changeset:" ` cd $MEF90_DIR && hg log -l 1 | grep changeset | awk '{print $2}'`
echo "Tao-lite changeset:" ` cd $TAO_DIR && hg log -l 1 | grep changeset | awk '{print $2}'`
