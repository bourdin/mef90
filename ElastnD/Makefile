PROJ      = Elast
MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 
M_VARFRACSTRUCTDIR  = ../m_VarFrac_Struct/${PETSC_ARCH}
M_VARFRACSTRUCTOBJS = ${M_VARFRACSTRUCTDIR}/m_VarFrac_Struct.o

all: ${PROJ}${DIM}

#include ${PETSC_DIR}/conf/rules
#include ${PETSC_DIR}/conf/variables
include ${TAO_DIR}/bmake/tao_common
include ../Makefile.inc

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/${DIM}$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

${PROJ}${DIM}: ${PETSC_ARCH} ${PROJ}.o   ${PETSC_ARCH/${DIM}${PROJ}.o m_${PROJ}.o ${PETSC_ARCH}/${DIM}m_${PROJ}.o
	${FLINKER} -o ${PETSC_ARCH}/${PROJ}${DIM} ${PETSC_ARCH}/${DIM}${PROJ}.o  ${PETSC_ARCH}/${DIM}m_${PROJ}.o  ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 

${PROJ}.o:   ${PETSC_ARCH} m_${PROJ}.o ${PETSC_ARCH}/${DIM}m_${PROJ}.o

m_${PROJ}.o: ${PETSC_ARCH}

clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}

debug:
	@echo FC is ${FC}
	@echo FC_FLAGS is ${FC_FLAGS}
	@echo MEF90_FC_FLAGS is ${MEF90_FC_FLAGS}
	@echo PETSC_MAKE_STOP_ON_ERROR is ${PETSC_MAKE_STOP_ON_ERROR}
	@echo FCPPFLAGS is ${FCPPFLAGS}
	@echo FC_MODULE_FLAG is ${FC_MODULE_FLAG}
	@echo FC_MODULE_OUTPUT_FLAG is ${FC_MODULE_OUTPUT_FLAG}
