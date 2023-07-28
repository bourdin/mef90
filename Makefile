all: MEF90 m_HeatXfer HeatXfer m_DefMech ThermoElasticity vDef Utils

mef90version.h: chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h

MEF90: mef90version.h chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C ${PETSC_ARCH}/objs -f ../../MEF90/Makefile MEF90

m_HeatXfer: mef90version.h MEF90 chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C ${PETSC_ARCH}/objs -f ../../m_HeatXfer/Makefile m_HeatXfer

HeatXfer: mef90version.h MEF90 m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C ${PETSC_ARCH}/objs -f ../../HeatXfer/Makefile HeatXfer

m_DefMech: mef90version.h MEF90 chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C ${PETSC_ARCH}/objs -f ../../m_DefMech/Makefile m_DefMech

m_Elasticity: mef90version.h MEF90 chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C ${PETSC_ARCH}/objs -f ../../m_Elasticity/Makefile m_Elasticity

ThermoElasticity: mef90version.h MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C ${PETSC_ARCH}/objs -f ../../ThermoElasticity/Makefile ThermoElasticity

ThermoElastoPlasticity: mef90version.h MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C ${PETSC_ARCH}/objs -f ../../ThermoElastoPlasticity/Makefile ThermoElastoPlasticity

WorkControlled: mef90version.h MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C ${PETSC_ARCH}/objs -f ../../WorkControlled/Makefile WorkControlled

vDef: mef90version.h MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C ${PETSC_ARCH}/objs -f ../../vDef/Makefile vDef
#	-@make -C ${PETSC_ARCH}/objs -f ../../vDef/Makefile vDef vDefP vDefUpa vDefBT vDefHF

Utils: mef90version.h MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C ${PETSC_ARCH}/objs -f ../../Utils/Makefile all

test: MEF90 chkpaths
	-@make -s -C HeatXfer test
	-@make -s -C ThermoElasticity test
	-@make -s -C vDef test

runtests: MEF90 chkpaths
	-@make -C ${PETSC_ARCH}/objs -f ../../Tests/Makefile runall

chkpaths: ${PETSC_ARCH}/objs ${PETSC_ARCH}/bin
${PETSC_ARCH}/objs:
	-@mkdir -p ${PETSC_ARCH}/objs
${PETSC_ARCH}/bin:
	-@mkdir -p ${PETSC_ARCH}/bin

doc: doc/vDef.pdf doc/vDef.tex
	-@echo "Building documentation"
	-@cd doc; pdflatex -shell-escape vDef.tex; bibtext vDef.tex; pdflatex -shell-escape vDef.tex; pdflatex -shell-escape vDef.tex;

tarball: clean
	$(eval gitver := $(shell git describe --dirty --always --tags))
	gtar --transform 's,^\.,mef90-${gitver},' --exclude .nfs* --exclude-backups --exclude=.git* --exclude=objs --exclude=lib --exclude=*pyc --exclude=bin/*/HeatXfer --exclude=bin/*/vDef --exclude=bin/*/ThermoElasticity --exclude=*so --exclude=*dylib -zcvhf ../mef90-${gitver}.tgz .

clean:
	-@rm ${MEF90_DIR}/mef90version.h
	-@rm -Rf ${PETSC_ARCH}/objs
	-@rm -Rf ${PETSC_ARCH}/bin
	-@make -C HeatXfer testclean
	-@make -C ThermoElasticity testclean
	-@make -C vDef testclean
	-@make -C Tests clean
