# generate the mesh from the geo file
gmsh -2 Pacman.geo
# generate a .gen version of teh mesh as the python scripts cannot handle .msh files
rm Pacman.gen; gmsh2exo -geometry Pacman.msh -result Pacman.gen
# write the displacement field in the result file
vDefScalingBC.py -i Pacman.gen -o Pacman_out.gen --time_min 0. --time_max 2. --time_numstep 21 --cs 2
# run the computation
srun -n 4 vDef -geometry Pacman.msh -result Pacman_out.gen -options_file PacmanDamaged.yaml 
 
