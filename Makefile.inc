## -*- Mode: Makefile; -*-
## vim: set ft=automake :

OBJDIR = ${MEF90_DIR}/obj/${PETSC_ARCH}/
LIBDIR = ${OBJDIR}
BINDIR = ${MEF90_DIR}/bin/${PETSC_ARCH}/

OBJS =  ${addprefix ${OBJDIR},${addprefix ${DIM},${SRCS:.F90=.o}}}
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
	FFLAGS += -L${LIBDIR} -lmef90 
	FFLAGS += -lvarstruct
ifneq (,$(findstring Darwin,$(PETSC_ARCH)))
	FFLAGS += ${OBJDIR}/m_MEF_Parameters.o 
endif 
endif

#This should be improved !!!!! 
#ifdef WITH_HEAT
ifeq (${OPTS}, Heat)
	FFLAGS += ${OBJDIR}/${DIM}m_TransientHeat.o ${OBJDIR}/${DIM}m_MEF_Integrate.o 
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
#2011/12/09 why do I have to add next line to compile varFrac with Heat module ? 	
#Order of flags is important ???? 
	FFLAGS += -L${LIBDIR} -lmef90 
#2011/12/28 and now this line for the linkage of m_transientHeat with VarStructm_Heat_struct ? 
	FFLAGS += -lvarstruct
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

${LIB}: lib${LIB}.a 

lib${LIB}.a: ${OBJS}
	${AR} ${AR_FLAGS} ${OBJDIR}/$@  $?

#${PETSC_FORTRAN_LIB} ${PETSC_LIB} must be at the end
#FLINKER = ${FC_LINKER} ${FC_LINKER_FLAGS} ${FFLAGS} 
${BIN}: ${OBJS} 
	${FC} ${MEF90_FLAGS} ${MEF90_INCS} -c -o ${OBJDIR}/${DIM}$@.o $@.F90
	${FC_LINKER} ${FC_LINKER_FLAGS} -o ${BINDIR}/$@${DIM}${OPTS}  $? ${addprefix ${OBJDIR},${addprefix ${DIM}, $@.o}} ${FFLAGS} 

${OBJDIR}${DIM}%.o : %.F90
	${FC} ${MEF90_FLAGS} ${MEF90_INCS} -c -o $@ $<


clean:
	-cd ${OBJDIR} && rm -f *.o *.mod  
	-cd ${LIBDIR} && rm -f *.a   
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
	@echo DIM is ${DIM}
	@echo OPTS is ${OPTS}

.PHONY: all
.SUFFIXES: .o .so .F90 

########################################################################
#                           Includes 
########################################################################

include ${PETSC_DIR}/conf/variables

ifdef TAO
	include ${TAO_DIR}/bmake/tao_common_variables
endif 
