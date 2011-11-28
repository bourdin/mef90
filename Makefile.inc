## -*- Mode: Makefile; -*-
## vim: set ft=automake :

OBJDIR = ${MEF90_DIR}/obj/${PETSC_ARCH}/
LIBDIR = ${OBJDIR}
BINDIR = ${MEF90_DIR}/bin/${PETSC_ARCH}/

OBJS =  ${addprefix ${OBJDIR},${addprefix ${DIM},${SRCS:.F90=.o}}}
SOBJS = ${OBJS:.o=.so}
TARGET = $(addprefix ${BINDIR},${BIN})
##################################################################
#                        Variables 
#Setting flags and include statements 
#################################################################

MEF90_FLAGS += ${FC_MODULE_OUTPUT_FLAG} ${OBJDIR}
MEF90_INCS += ${PETSC_FC_INCLUDES}  

#FFLAGS is used in FLINKER def in petsc/rules
#No link edition in libs 
ifdef BIN
	FFLAGS = ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
	FFLAGS += -L${LIBDIR} -lMEF90 
	FFLAGS += -lVar_Struct 
endif

#This should be improved !!!!! 
ifdef WITH_HEAT
	FFLAGS += ${OBJDIR}/${DIM}m_Poisson.o ${OBJDIR}/${DIM}m_TransientHeat.o ${OBJDIR}/${DIM}m_MEF_Integrate.o 
endif


#This is to accept long lines 
ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FLAGS +=  -ffixed-line-length-none  -ffree-line-length-none -ffree-form 
endif

ifneq (,$(findstring xlf,$(FC)))
	FC_D_FLAG=-WF,-D
else
  	FC_D_FLAG=-D
endif

#Dimension Flags 2D / 3D 
ifdef DIM 
	ifneq (,$(findstring gcc,$(PETSC_ARCH)))
		MEF90_FLAGS += ${FC_D_FLAG}PB_${DIM} 
	else
		MEF90_FLAGS += -DPB_${DIM}  
	endif
endif

#TAO Fllags 
ifdef TAO
	FFLAGS += ${TAO_LIB} ${TAO_FORTRAN_LIB} 
	ifneq (,$(findstring gcc,$(PETSC_ARCH)))
		MEF90_FLAGS +=  ${FC_D_FLAG}WITH_TAO 
	else
		MEF90_FLAGS += -DWITH_TAO 
	endif
endif 

################################################################## 
#                         Rules
##################################################################
all: ${BINDIR} ${OBJDIR} ${BIN} ${LIB}

${OBJDIR}:
	mkdir -p $@

${BINDIR}:
	mkdir -p $@

${LIB}: lib${LIB}.a lib${LIB}.so

lib${LIB}.a: ${OBJS}
	${AR} ${AR_FLAGS} ${OBJDIR}/$@  $?

#TODO option -shared take automatic variable from Petsc
#same for -DPic ? 
lib${LIB}.so: ${SOBJS}
	${FLINKER} -shared -o ${OBJDIR}/$@  $?  

#${PETSC_FORTRAN_LIB} ${PETSC_LIB} must be at the end
#FLINKER = ${FC_LINKER} ${FC_LINKER_FLAGS} ${FFLAGS} 
${BIN}: ${OBJS} 
	${FC} ${MEF90_FLAGS} ${MEF90_INCS} -c -o ${OBJDIR}/${DIM}$@.o $@.F90
	${FC_LINKER} ${FC_LINKER_FLAGS} -o ${BINDIR}/$@${DIM}  $? ${addprefix ${OBJDIR},${addprefix ${DIM}, $@.o}} ${FFLAGS} 

${OBJDIR}${DIM}%.o : %.F90
	${FC} ${MEF90_FLAGS} ${MEF90_INCS} -c -o $@ $<

${OBJDIR}${DIM}%.so : %.F90
	${FC} ${FFLAGS} ${MEF90_FLAGS} ${MEF90_INCS} -DPIC -c -o $@ $<


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

########################################################################
#                           Includes 
########################################################################

include ${PETSC_DIR}/conf/variables

ifdef TAO
	include ${TAO_DIR}/bmake/tao_common_variables
endif 
