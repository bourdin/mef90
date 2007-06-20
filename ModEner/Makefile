PROJ=ModEner

$(PROJ): ${PROJ}.o chkopts
	-${FLINKER} -o $(PROJ) $(PROJ).o ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

include ${PETSC_DIR}/bmake/common/base

really_clean: clean
	- rm -f ${PROJ}-${VER} *~ core

tarball: really_clean
	tar zcvf ${PROJ}-${VER}-`date '+%Y-%m-%d'`-`hostname -s`.tgz *[^tgz] > /dev/null

