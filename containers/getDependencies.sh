#!/bin/bash
if [ -d tarball ]; then
	cd tarball
else
	mkdir tarball
	cd tarball
fi

echo 'Downloading netcdf  '; curl -O https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdf-c-4.7.4.tar.gz
echo 'Downloading chaco   '; curl -O http://ftp.mcs.anl.gov/pub/petsc/externalpackages/Chaco-2.2.tar.gz
echo 'Downloading exodusii'; curl -O https://github.com/gsjaardema/seacas/archive/0ac6baa46262ab17a7b3c0a07b8832ce8d088324.tar.gz
echo 'Downloading yaml    '; curl -O http://pyyaml.org/download/libyaml/yaml-0.1.4.tar.gz
echo 'Downloading hypre   '; curl -O http://ftp.mcs.anl.gov/pub/petsc/externalpackages/hypre-2.8.0b.tar.gz
echo 'Downloading sowing  '; curl -O http://ftp.mcs.anl.gov/pub/petsc/externalpackages/sowing-1.1.26-p1.tar.gz
echo 'Downloading ml      '; curl -O http://ftp.mcs.anl.gov/pub/petsc/externalpackages/ml-6.2-p3.tar.gz
echo 'Downloading parmetis'; curl -O http://ftp.mcs.anl.gov/pub/petsc/externalpackages/parmetis-4.0.2-p3.tar.gz
echo 'Downloading hdf5    '; curl -O http://ftp.mcs.anl.gov/pub/petsc/externalpackages/hdf5-1.8.8-p1.tar.gz
echo 'Downloading metis   '; curl -O http://ftp.mcs.anl.gov/pub/petsc/externalpackages/metis-5.0.2-p3.tar.gz
#echo 'Downloading gmsh    '; curl -O http://gmsh.info/bin/Linux/gmsh-4.7.1-Linux64.tgz

if [ -d petsc-3.3-mef90 ]; then
	cd petsc-3.3-mef90
	git pull
	cd ..
else
	git clone --single-branch --branch mef90-3.3 https://gitlab.com/blaisebourdin/petsc.git petsc-3.3-mef90
fi

if [ -d mef90 ]; then
	cd mef90
	git pull
	cd ..
else
	git clone https://github.com/bourdin/mef90.git
fi

if [ -d snlp ]; then
	cd snlp
	git pull
	cd ..
else
	git clone https://bitbucket.org/bourdin/snlp.git
fi

cd ..
