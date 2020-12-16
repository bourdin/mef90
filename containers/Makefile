

.PHONY: centos7openmpi313base centos7openmpi313petsc33 mef90centos7openmpig mef90centos7openmpiO 
.PHONY: centos7mpich32base    centos7mpich32petsc33    mef90centos7mpichg   mef90centos7mpichO
.PHONY: build buildopenmpi buildmpich
.PHONY: publish publishopenmpi publishmpich

publish: publishopenmpi publishmpich

build: buildopenmpi buildmpich
buildopenmpi: mef90centos7openmpig mef90centos7openmpio
buildmpich: mef90centos7mpichg mef90centos7mpicho 


centos7openmpi313base: Dockerfile.centos7-petsc3.3-openmpi-3.1.3
	-@echo "Building centos openMPI base container"
	-@docker build -t $@ -f Dockerfile.centos7-petsc3.3-openmpi-3.1.3 --target $@ .
	-@docker tag $@:latest bourdin/$@:latest
centos7openmpi313petsc33: Dockerfile.centos7-petsc3.3-openmpi-3.1.3 centos7openmpi313base
	-@echo "Building openMPI PETSC-3.3 container"
	-@docker build -t $@ -f Dockerfile.centos7-petsc3.3-openmpi-3.1.3 --target $@ .
	-@docker tag $@:latest bourdin/$@:latest

mef90centos7openmpig: Dockerfile.mef90centos7openmpi-g centos7openmpi313petsc33
	-@echo "Building mef90 openMPI debugging container"
	-@docker build -t $@ -f Dockerfile.mef90centos7openmpi-g .
	-@docker tag $@:latest bourdin/$@:latest

mef90centos7openmpio: Dockerfile.mef90centos7openmpi-O centos7openmpi313petsc33
	-@echo "Building mef90 openMPI optimized container"
	-@docker build -t $@ -f Dockerfile.mef90centos7openmpi-O .
	-@docker tag $@:latest bourdin/$@:latest

publishopenmpi: mef90centos7openmpig mef90centos7openmpio
	-@echo "Publishing openMPI containers to dockerhub"
	-@docker push bourdin/centos7openmpi313base:latest
	-@docker push bourdin/centos7openmpi313petsc33:latest
	-@docker push bourdin/mef90centos7openmpig:latest
	-@docker push bourdin/mef90centos7openmpio:latest

centos7mpich32base: Dockerfile.centos7-petsc3.3-mpich-3.2
	-@echo "Building centos mpich base container"
	-@docker build -t $@ -f Dockerfile.centos7-petsc3.3-mpich-3.2 --target $@ .

centos7mpich32petsc33: Dockerfile.centos7-petsc3.3-mpich-3.2 centos7mpich32base
	-@echo "Building mpich PETSC-3.3 container"
	-@docker build -t $@ -f Dockerfile.centos7-petsc3.3-mpich-3.2 --target $@ .
	-@docker tag $@:latest bourdin/$@:latest

mef90centos7mpichg: Dockerfile.mef90centos7mpich-g centos7mpich32petsc33
	-@echo "Building mef90 mpich debugging container"
	-@docker build -t $@ -f Dockerfile.mef90centos7mpich-g .
	-@docker tag $@:latest bourdin/$@:latest

mef90centos7mpicho: Dockerfile.mef90centos7mpich-O centos7mpich32petsc33
	-@echo "Building mef90 mpich optimized container"
	-@docker build -t $@ -f Dockerfile.mef90centos7mpich-O .
	-@docker tag $@:latest bourdin/$@:latest

publishmpich: mef90centos7mpichg mef90centos7mpicho
	-@echo "Publishing mpich containers to dockerhub"
	-@docker push bourdin/centos7mpich32base:latest
	-@docker push bourdin/centos7mpich32petsc33:latest
	-@docker push bourdin/mef90centos7mpichg:latest
	-@docker push bourdin/mef90centos7mpicho:latest
