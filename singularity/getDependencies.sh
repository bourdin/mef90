#!/bin/bash
if [ -d tarball ]; then
	cd tarball
else
	mkdir tarball
	cd tarball
fi
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.3-p7.tar.gz
wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/Chaco-2.2.tar.gz
wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/exodusii-5.22b.tar.gz
wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/hypre-2.8.0b.tar.gz
wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/ml-6.2-win.tar.gz
wget ftp://ftp.unidata.ucar.edu/pub/netcdf/netcdf-4.3.3.1.tar.gz
wget http://pyyaml.org/download/libyaml/yaml-0.1.4.tar.gz
wget http://gmsh.info/bin/Linux/gmsh-3.0.6-Linux64.tgz
#wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/boost_minimal_1_42_0.tar.gz
#wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/sowing-1.1.16d.tar.gz

if [ -d mef90 ]; then
	cd mef90
	git pull
	cd ..
else
	git clone https://github.com/bourdin/mef90.git mef90
fi

if [ -d snlp ]; then
	cd snlp
	git pull
	cd ..
else
    git clone https://github.com/bourdin/SNLP.git snlp
fi

cd ..
