## -*- Mode: Makefile; -*-
## vim: set ft=automake :

OBJS = ${SRCS:.F90=.o}
SOBJS = ${OBJS:.o=.so}
#DIMOBJS = ${addprefix ${DIM}, ${OBJS}}

OBJDIR = ${MEF90_DIR}/obj/${PETSC_ARCH}/
LIBDIR = ${OBJDIR}
BINDIR = ${MEF90_DIR}/bin/${PETSC_ARCH}/
#OBJDIR = .
#DIROBJS = ${addprefix ${OBJDIR}, ${OBJS}}
#DIRSOBJS = ${addprefix ${OBJDIR}, ${SOBJS}}

MEF90_FLAGS += ${FC_MODULE_OUTPUT_FLAG} ${OBJDIR}
MEF90_INCS += ${PETSC_FC_INCLUDES}  

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FLAGS +=  -ffixed-line-length-none  -ffree-line-length-none -ffree-form 
endif

all: ${BINDIR} ${OBJDIR} ${BIN} ${LIB}

${OBJDIR}:
	mkdir -p $@

${BINDIR}:
	mkdir -p $@

${LIB}: lib${LIB}.a lib${LIB}.so

lib${LIB}.a: ${OBJS} 
	${AR} ${AR_FLAGS} ${OBJDIR}/$@ ${addprefix ${OBJDIR}, $?}

#TODO option -shared take automatic varaible from Petsc
#same for -DPic ? 
lib${LIB}.so: ${SOBJS}
	${FLINKER} -shared -o ${OBJDIR}/$@ ${addprefix ${OBJDIR}, $?}  

#-static library ? 
${BIN}: ${OBJS} 
	${FLINKER} -o ${BINDIR}/$@${DIM} ${addprefix ${OBJDIR},${addprefix ${DIM}, $?}} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

.F90.o:
	${FC} ${MEF90_FLAGS} ${MEF90_INCS} -c -o ${OBJDIR}/${DIM}$@ $<

.F90.so:
	${FC} ${FFLAGS} ${MEF90_FLAGS} ${MEF90_INCS} -DPIC -c -o ${OBJDIR}/$@ $<


clean:
	-cd ${OBJDIR} && rm -f *.o *.mod  
	-cd ${LIBDIR} && rm -f *.a *.so  
	-rm -rf ${OBJDIR}
	-rm -rf ${LIBDIR}
	-rm -rf ${BINDIR}

debug::
	@echo MEF90_FC_FLAGS is ${MEF90_FLAGS}
	@echo FC is ${FC}
	@echo FC_FLAGS is ${FC_FLAGS}
	@echo PETSC_MAKE_STOP_ON_ERROR is ${PETSC_MAKE_STOP_ON_ERROR}
	@echo FCPPFLAGS is ${FCPPFLAGS}
	@echo FC_MODULE_FLAG is ${FC_MODULE_FLAG}
	@echo MODULE_OUTPUT_FLAG is ${FC_MODULE_OUTPUT_FLAG}


.PHONY: all
.SUFFIXES: .o .so .F90 


include ${PETSC_DIR}/conf/variables
