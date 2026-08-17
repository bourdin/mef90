DIRS: MEF90 m_HeatXfer HeatXfer m_DefMech ThermoElasticity vDef Utils

all: $(DIRS)

.PHONY: all $(DIRS) ford fordclean fordpublish

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

FORD_VENV = ${MEF90_DIR}/.venv-ford
FORD = ${FORD_VENV}/bin/ford

# PETSc source tree, read to link the PETSc calls in the source listings to
# their manual pages on petsc.org. PETSC_DIR is right when PETSc was built in
# place; point PETSC_SRC at the sources when PETSc is a prefix install. The
# PETSc links are simply left out if neither holds a src directory.
PETSC_SRC ?= ${PETSC_DIR}
FORD_PETSC = $(if $(strip ${PETSC_SRC}),--petsc-src ${PETSC_SRC},)

${FORD}:
	-@echo "Creating the FORD virtualenv in ${FORD_VENV}"
	python3 -m venv ${FORD_VENV} || exit 1
	${FORD_VENV}/bin/pip install --quiet ford --editable ${MEF90_DIR}/doc/doi-link || exit 1

ford: ${FORD}
	-@echo "Building the source documentation in ${MEF90_DIR}/doc/html"
	${FORD} ${MEF90_DIR}/ford.md || exit 1
	${MEF90_DIR}/bin/ford-xref ${MEF90_DIR}/doc/html ${FORD_PETSC} || exit 1

fordclean:
	-@rm -Rf ${MEF90_DIR}/doc/html

# Publishing to GitHub Pages. The generated pages live in their own repository,
# checked out as an ordinary clone next to this one, so there is nothing to
# understand beyond `git clone`: look inside MEF90_DOCS at any point and it is
# just the website, on a branch you can inspect, diff, and push by hand.
#
#   git clone ${MEF90_DOCS_URL} ${MEF90_DOCS}    (once)
#   make fordpublish
MEF90_DOCS_URL ?= ssh://github.com/bourdin/mef90-doc
MEF90_DOCS ?= ${MEF90_DIR}/../mef90-doc

# Files belonging to the repository rather than to the generated site, which
# --delete would otherwise carry off on every publish. CNAME is written by
# GitHub when a custom domain is set in the Pages settings, and losing it takes
# mef90.org down with it.
MEF90_DOCS_KEEP = --exclude .git --exclude CNAME --exclude README.md \
                  --exclude LICENSE

fordpublish: ford
	@test -d ${MEF90_DOCS}/.git || { \
	  echo "No clone at ${MEF90_DOCS}."; \
	  echo "Create the repository on GitHub, then:"; \
	  echo "    git clone ${MEF90_DOCS_URL} ${MEF90_DOCS}"; \
	  exit 1; }
	cd ${MEF90_DOCS} && git pull --ff-only
	rsync -a --delete ${MEF90_DOCS_KEEP} ${MEF90_DIR}/doc/html/ ${MEF90_DOCS}/
	touch ${MEF90_DOCS}/.nojekyll
	cd ${MEF90_DOCS} && git add -A && \
	  (git diff --cached --quiet || \
	   git commit -m "mef90 $$(cd ${MEF90_DIR} && git describe --dirty --always --tags)") && \
	  git push

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
