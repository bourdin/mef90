all: TestMaterials TestMEF90Ctx

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables
include ${MEF90_DIR}/MEF90/Makefile.include



${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}
	
TestMaterials:  ${PETSC_ARCH} TestMaterials.o ${PETSC_ARCH}/TestMaterials.o
	@${FLINKER} -o ${PETSC_ARCH}/TestMaterials ${PETSC_ARCH}/TestMaterials.o ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

TestMEF90Ctx:  ${PETSC_ARCH} TestMEF90Ctx.o ${PETSC_ARCH}/TestMEF90Ctx.o
	@${FLINKER} -o ${PETSC_ARCH}/TestMEF90Ctx ${PETSC_ARCH}/TestMEF90Ctx.o ${MEF90_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

runTestMaterials: TestMaterials
	@echo "Running " $@
	-@ARGS='';\
	echo ${MPIEXEC} -n 2 ${PETSC_ARCH}/$< $${ARGS};\
	${MPIEXEC} -n 2 ${PETSC_ARCH}/$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

runTestMEF90Ctx: TestMEF90Ctx
	@echo "Running " $@
	-@ARGS='-verbose 1';\
	echo ${MPIEXEC} -n 2 ${PETSC_ARCH}/$< $${ARGS};\
	${MPIEXEC} -n 2 ${PETSC_ARCH}/$< $${ARGS} > $@.out 2>&1;\
	if (${DIFF} -B results/$@.out $@.out) then true; \
	else echo ${PWD} ; echo "Possible problem with with $<, diffs above \n========================================="; fi;\
	${RM} -f $@.out

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}



