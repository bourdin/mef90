# Using the mef90 docker container

This quick tutorial will demonstrate how to run some of the built-in examples in a docker container, and save the result to your local drive:
It assumes some familiarity with the basics of docker, for instance [this tutorial](https://docs.docker.com/get-started/).


## Starting the docker container

Start a container by typing the following command in a shell window (there may be a way to do this from the docker GUI).
```
mef90:~$ docker run -ti --rm -v ${HOME}/mef90Demo:/home/mef90  bourdin/mef90centos7mpicho
[root@05a4fc873968 /]# 
```
Note the ```-v ${HOME}/mef90Demo:/home/mef90``` argument binds the folder ```${HOME}/mef90Demo``` of my _local_ computer to the folder ```/home/mef90``` of the container. In other words, any file in ```${HOME}/mef90Demo``` of my _local_ computer can be accessed in the  ```/home/mef90``` folder of the container. Only the files created in this folder of teh container will remain once the container is stopped.


## Running the uniaxial tension example

* Move to the folder containing the example inside the container:
   ```
   cd $MEF90_DIR/doc/Examples/UniaxialTension 
   ```
* Generate a mesh of the geometry described in ```UniaxialTension2D.geo``` in ```.msh``` format using ```gmsh```:
    ```
    [root@05a4fc873968 UniaxialTension]# gmsh -2 -format  msh2 UniaxialTension2D.geo
    Info    : Running 'gmsh -2 -format msh2 UniaxialTension2D.geo' [Gmsh 4.5.5, 1 node, max. 1 thread]
    Info    : Started on Thu Feb  4 15:00:40 2021
    Info    : Reading 'UniaxialTension2D.geo'...
    Info    : Done reading 'UniaxialTension2D.geo'
    Info    : Meshing 1D...
    Info    : [  0 %] Meshing curve 1 (Line)
    Info    : [ 30 %] Meshing curve 2 (Line)
    Info    : [ 50 %] Meshing curve 3 (Line)
    Info    : [ 80 %] Meshing curve 4 (Line)
    Info    : Done meshing 1D (0.00264 s)
    Info    : Meshing 2D...
    Info    : Meshing surface 6 (Plane, Frontal)
    Info    : Done meshing 2D (0.775613 s)
    Info    : 15122 nodes 30246 elements
    Info    : Writing 'UniaxialTension2D.msh'...
    Info    : Done writing 'UniaxialTension2D.msh'
    Info    : Stopped on Thu Feb  4 15:00:41 2021
    ```
   From there, you should see a file named ```UniaxialTension2D.msh``` in the current folder.
* Convert this mesh into an exodusII-formatted mesh that can be used by ```mef90```:
    ```
    [root@b2092c5310ee UniaxialTension]# gmsh2exo.py UniaxialTension2D.msh UniaxialTension2D.gen

    You are using exodus.py v 1.19.1 (seacas-py3), a python wrapper of some of the exodus library.

    Copyright (c) 2013, 2014, 2015, 2016, 2017, 2018, 2019 National Technology &
    Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of
    Contract DE-NA0003525 with NTESS, the U.S. Government retains certain
    rights in this software.

    Opening exodus file: UniaxialTension2D.gen
    Closing exodus file: UniaxialTension2D.gen
    ```
* Run vDef:
    ```
    [root@b2092c5310ee UniaxialTension]#  mpirun -np 4 vDef -geometry UniaxialTension2D.gen -options_file_yaml UniaxialTensionAT2-Equality.yaml -result /home/mef90/UniaxialTension2D_out.gen 
    ```
    In this example, ```vDef``` will run on 4 cores (the ```-np 4``` option of ```mpirun```) and the result files will be placed in the ```/home/mef90/``` folder f the container, which corresponds to the ```${HOME}/mef90Demo``` folder of the host.

* Viewing results:
    <!-- * An energy plot can be generated in ```/home/mef90/``` using the command
    ```
    plotener.py /home/mef90/UniaxialTension2D_out.ener -o /home/mef90/UniaxialTension2D_out.pdf
    ``` -->
    The evolution of the phase field variable can be viewed using [visit](http://visit.llnl.gov) or [paraview](http://www.paraview.org) (or any other software capable of handling exodusII files) to render a pseudo-color plot of the ```Damage``` variable.