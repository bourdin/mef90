all: MEF90 m_HeatXfer HeatXfer m_DefMech ThermoElasticity vDef

${MEF90_DIR}/.hg/dirstate:
	@mkdir -p ${MEF90_DIR}/.hg; touch ${MEF90_DIR}/.hg/dirstate
	
mef90version.h: ${MEF90_DIR}/.hg/dirstate
	@bin/makeversion.sh ${MEF90_DIR}/mef90version.h

MEF90: mef90version.h chkpaths
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../MEF90/Makefile MEF90

m_HeatXfer: MEF90 chkpaths
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../m_HeatXfer/Makefile m_HeatXfer

HeatXfer: MEF90 m_HeatXfer chkpaths
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../HeatXfer/Makefile HeatXfer

m_DefMech: MEF90 chkpaths
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../m_DefMech/Makefile m_DefMech

m_Elasticity: MEF90 chkpaths
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../m_Elasticity/Makefile m_Elasticity

ThermoElasticity: MEF90 m_DefMech m_HeatXfer chkpaths
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../ThermoElasticity/Makefile ThermoElasticity

vDef: MEF90 m_DefMech m_HeatXfer chkpaths
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../vDef/Makefile vDef

test: MEF90 chkpaths
	-@make -s -C HeatXfer test
	-@make -s -C ThermoElasticity test
	-@make -s -C vDef test

runtests: MEF90 chkpaths
	-@make -C objs/${PETSC_ARCH} -f ../../Tests/Makefile runall

chkpaths: objs/${PETSC_ARCH} bin/${PETSC_ARCH} lib/${PETSC_ARCH}
objs/${PETSC_ARCH}:
	-@mkdir -p objs/${PETSC_ARCH}
bin/${PETSC_ARCH}:
	-@mkdir -p bin/${PETSC_ARCH}
lib/${PETSC_ARCH}:
	-@mkdir -p lib/${PETSC_ARCH}

tarball: clean
	$(eval hgver := $(shell hg parents | head -1 | cut -d : -f 2 | tr -d ' '))
	gtar --transform 's,^\.,mef90-${hgver},' --exclude-backups --exclude=.hg* --exclude=objs --exclude=lib --exclude=*pyc --exclude=bin/*/HeatXfer --exclude=bin/*/vDef --exclude=bin/*/ThermoElasticity --exclude=*so --exclude=*dylib -zcvhf ../mef90-${hgver}.tgz .

clean:
	-@rm -Rf objs/${PETSC_ARCH}
	-@rm -Rf bin/${PETSC_ARCH}
	-@make -C HeatXfer testclean
	-@make -C ThermoElasticity testclean
	-@make -C vDef testclean
	-@make -C Tests clean
