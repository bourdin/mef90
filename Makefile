all: SimplePoissonNG

SimplePoissonNG: SimplePoissonNG2D SimplePoissonNG3D

MEF90: chkpaths
	-@echo "Building $@"
	-@make -C objs/${PETSC_ARCH} -f ../../MEF90/Makefile MEF90

m_HeatXfer: MEF90 chkpaths
	-@echo "Building $@"
	-@make -C objs/${PETSC_ARCH} -f ../../m_HeatXfer/Makefile m_HeatXfer

SimplePoissonNG2D: MEF90 chkpaths
	-@echo "Building $@"
	-@make -C objs/${PETSC_ARCH} -f ../../SimplePoissonNG/Makefile SimplePoissonNG MEF90_DIM=2

SimplePoissonNG3D: MEF90 chkpaths
	-@echo "Building $@"
	-@make -C objs/${PETSC_ARCH} -f ../../SimplePoissonNG/Makefile SimplePoissonNG MEF90_DIM=3

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


