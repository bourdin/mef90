all: MEF90 m_HeatXfer HeatXfer m_DefMech ThermoElasticity vDef WorkControlled

MEF90: chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../MEF90/Makefile MEF90

m_HeatXfer: MEF90 chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../m_HeatXfer/Makefile m_HeatXfer

HeatXfer: MEF90 m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../HeatXfer/Makefile HeatXfer

m_DefMech: MEF90 chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../m_DefMech/Makefile m_DefMech

m_Elasticity: MEF90 chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../m_Elasticity/Makefile m_Elasticity

ThermoElasticity: MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../ThermoElasticity/Makefile ThermoElasticity

ThermoElastoPlasticity: MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../ThermoElastoPlasticity/Makefile ThermoElastoPlasticity

WorkControlled: MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../WorkControlled/Makefile WorkControlled

vDef: MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../vDef/Makefile vDef vDefP vDefUpa vDefBT

YAMLValidator: MEF90 m_DefMech m_HeatXfer
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	-@make -C objs/${PETSC_ARCH} -f ../../Tests/Makefile YAMLValidator

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

doc: doc/vDef.pdf doc/vDef.tex
	-@echo "Building documentation"
	-@cd doc; pdflatex -shell-escape vDef.tex; bibtext vDef.tex; pdflatex -shell-escape vDef.tex; pdflatex -shell-escape vDef.tex;

tarball: clean
	$(eval gitver := $(shell git describe --dirty --always --tags))
	gtar --transform 's,^\.,mef90-${gitver},' --exclude .nfs* --exclude-backups --exclude=.git* --exclude=objs --exclude=lib --exclude=*pyc --exclude=bin/*/HeatXfer --exclude=bin/*/vDef --exclude=bin/*/ThermoElasticity --exclude=*so --exclude=*dylib -zcvhf ../mef90-${gitver}.tgz .

clean:
	-@rm ${MEF90_DIR}/mef90version.h
	-@rm -Rf objs/${PETSC_ARCH}
	-@rm -Rf bin/${PETSC_ARCH}
	-@make -C HeatXfer testclean
	-@make -C ThermoElasticity testclean
	-@make -C vDef testclean
	-@make -C Tests clean
