all: MEF90 m_HeatXfer HeatXfer m_DefMech ThermoElasticity mef90version.h

mef90version.h: ${MEF90_DIR}/.hg/dirstate
	@echo \#define MEF90_HGVER \"`hg tip | head -1 | cut -d : -f 2,3 | tr -d ' '`\" > ${MEF90_DIR}/mef90version.h

MEF90: chkpaths
	-@echo "Building $@"
	-@make -C objs/${PETSC_ARCH} -f ../../MEF90/Makefile MEF90

m_HeatXfer: MEF90 chkpaths
	-@echo "Building $@"
	-@make -C objs/${PETSC_ARCH} -f ../../m_HeatXfer/Makefile m_HeatXfer

HeatXfer: MEF90 m_HeatXfer chkpaths
	-@echo "Building $@"
	-@make -C objs/${PETSC_ARCH} -f ../../HeatXfer/Makefile HeatXfer

m_DefMech: MEF90 chkpaths
	-@echo "Building $@"
	-@make -C objs/${PETSC_ARCH} -f ../../m_DefMech/Makefile m_DefMech

m_Elasticity: MEF90 chkpaths
	-@echo "Building $@"
	-@make -C objs/${PETSC_ARCH} -f ../../m_Elasticity/Makefile m_Elasticity

ThermoElasticity: MEF90 m_DefMech chkpaths
	-@echo "Building $@"
	-@make -C objs/${PETSC_ARCH} -f ../../ThermoElasticity/Makefile ThermoElasticity

tests: MEF90 chkpaths
	-@echo "Building $@"
	-@make -C objs/${PETSC_ARCH} -f ../../Tests/Makefile all

runtests: MEF90 chkpaths
	-@make -C objs/${PETSC_ARCH} -f ../../Tests/Makefile runall

chkpaths: objs/${PETSC_ARCH} bin/${PETSC_ARCH}
objs/${PETSC_ARCH}:
	-@mkdir -p objs/${PETSC_ARCH}
bin/${PETSC_ARCH}:
	-@mkdir -p bin/${PETSC_ARCH}

clean:
	-@rm -Rf objs/${PETSC_ARCH}
	-@rm -Rf bin/${PETSC_ARCH}

