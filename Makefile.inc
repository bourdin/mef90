#vim :set ft=makefile

OBJS = ${SRCS:.F90=.o}
SOBJS = ${OBJS:.o=.so}

#OBJDIR = obj/${PETSC_ARCH}
OBJDIR = .

MEF90_FLAGS += ${FC_MODULE_OUTPUT_FLAG} ${OBJDIR}
MEF90_INCS += ${PETSC_FC_INCLUDES}  

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FLAGS +=  -ffixed-line-length-none  -ffree-line-length-none -ffree-form 
endif

all: ${OBJDIR} ${BIN} ${LIB}

${OBJDIR}:
	mkdir -p $@

${LIB}: lib${LIB}.a lib${LIB}.so

lib${LIB}.a: ${OBJS} 
	${AR} ${AR_FLAGS} ${OBJDIR}/$@ $?

lib${LIB}.so: ${SOBJS}
	${FLINKER} -shared -o ${OBJDIR}/$@ $?  

#-static library ? 
${BIN}: ${OBJS}
	${FLINKER}  -o $@ $? ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

.F90.o:
	${FC} ${MEF90_FLAGS} ${MEF90_INCS} -c -o ${OBJDIR}/$@ $<

.F90.so:
	${FC} ${FFLAGS} ${MEF90_FLAGS} ${MEF90_INCS} -DPIC -c -o ${OBJDIR}/$@ $<


clean:
	-cd ${OBJDIR} && rm -f *.o *.mod *.so 
	-rm -rf ${OBJDIR}

.PHONY: all
.SUFFIXES: .o .so .F90 


include ${PETSC_DIR}/conf/variables
