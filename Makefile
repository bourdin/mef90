DIRS: MEF90 m_HeatXfer HeatXfer m_DefMech ThermoElasticity vDef Utils

all: $(DIRS)

.PHONY: all $(DIRS)

$(DIRS):
	$(MAKE) -C $@ $(MFLAGS) all

mef90version.h: chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h

MEF90: mef90version.h chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	${MAKE} -C ${PETSC_ARCH}/objs -f ../../MEF90/Makefile MEF90 || exit 1

m_HeatXfer: mef90version.h MEF90 chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	${MAKE} -C ${PETSC_ARCH}/objs -f ../../m_HeatXfer/Makefile m_HeatXfer || exit 1

HeatXfer: mef90version.h MEF90 m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	${MAKE} -C ${PETSC_ARCH}/objs -f ../../HeatXfer/Makefile HeatXfer || exit 1

m_DefMech: mef90version.h MEF90 chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	${MAKE} -C ${PETSC_ARCH}/objs -f ../../m_DefMech/Makefile m_DefMech || exit 1

m_Elasticity: mef90version.h MEF90 chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	${MAKE} -C ${PETSC_ARCH}/objs -f ../../m_Elasticity/Makefile m_Elasticity || exit 1

ThermoElasticity: mef90version.h MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	${MAKE} -C ${PETSC_ARCH}/objs -f ../../ThermoElasticity/Makefile ThermoElasticity || exit 1

ThermoElastoPlasticity: mef90version.h MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	${MAKE} -C ${PETSC_ARCH}/objs -f ../../ThermoElastoPlasticity/Makefile ThermoElastoPlasticity || exit 1

WorkControlled: mef90version.h MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	${MAKE} -C ${PETSC_ARCH}/objs -f ../../WorkControlled/Makefile WorkControlled || exit 1

vDef: mef90version.h MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	${MAKE} -C ${PETSC_ARCH}/objs -f ../../vDef/Makefile vDef || exit 1
#	-@make -C ${PETSC_ARCH}/objs -f ../../vDef/Makefile vDef vDefP vDefUpa vDefBT vDefHF

Utils: mef90version.h MEF90 m_DefMech m_HeatXfer chkpaths
	-@bin/makeversion.sh ${MEF90_DIR}/mef90version.h
	-@echo "Building $@ with PETSC_ARCH=${PETSC_ARCH}"
	${MAKE} -C ${PETSC_ARCH}/objs -f ../../Utils/Makefile all

test: MEF90 chkpaths
	${MAKE} -s -C HeatXfer test
	${MAKE} -s -C ThermoElasticity test
	${MAKE} -s -C vDef test

runtests: MEF90 chkpaths
	${MAKE} -C ${PETSC_ARCH}/objs -f ../../Tests/Makefile runall

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
	${MAKE} -C HeatXfer testclean
	${MAKE} -C ThermoElasticity testclean
	${MAKE} -C vDef testclean
	${MAKE} -C Tests clean
