#vim :set ft=makefile

OBJS = ${SRCS:.F90=.o}
SOBJS = ${OBJS:.o=.so}
	

MEF90_FLAGS += ${FC_MODULE_OUTPUT_FLAG} . 
MEF90_INCS += ${PETSC_FC_INCLUDES}  

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FLAGS +=  -ffixed-line-length-none  -ffree-line-length-none -ffree-form 
endif

all: ${BIN} ${LIB}

${LIB}: lib${LIB}.a lib${LIB}.so

lib${LIB}.a: ${OBJS} 
	${AR} ${AR_FLAGS} $@ $?

lib${LIB}.so: ${SOBJS}
	${FLINKER} -shared -o $@ $?  

#-static library ? 
${BIN}: ${OBJS}
	${FLINKER}  -o $@ $?  

.F90.o:
	${FC} ${MEF90_FLAGS} ${MEF90_INCS} -c -o $@ $<

.F90.so:
	${FC} ${FFLAGS} ${MEF90_FLAGS} ${MEF90_INCS} -DPIC -c -o $@ $<


clean:
	-rm -f *.o *.mod *.so 

.PHONY: all
.SUFFIXES: .o .so .F90 


include ${PETSC_DIR}/conf/variables
