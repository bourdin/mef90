#!/bin/bash
if [ -d tarball ]; then
	cd tarball
else
	mkdir tarball
	cd tarball
fi
wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/Chaco-2.2.tar.gz
wget https://github.com/gsjaardema/seacas/archive/b37ad5be658bd447717faa714f77ea3e4aa9537b.tar.gz
wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/hypre-2.8.0b.tar.gz
wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/ml-6.2-win.tar.gz
wget https://www.unidata.ucar.edu/downloads/netcdf/ftp/netcdfc-4.7.4.tar.gz
wget http://pyyaml.org/download/libyaml/yaml-0.1.4.tar.gz
wget http://gmsh.info/bin/Linux/gmsh-4.5.6-Linux64.tgz
#wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/boost_minimal_1_42_0.tar.gz
#wget http://ftp.mcs.anl.gov/pub/petsc/externalpackages/sowing-1.1.16d.tar.gz

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
