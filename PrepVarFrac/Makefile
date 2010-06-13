MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 
M_VARFRACSTRUCTDIR  = ../m_VarFrac_Struct/${PETSC_ARCH}
M_VARFRACSTRUCTOBJS = ${M_VARFRACSTRUCTDIR}/m_VarFrac_Struct.o

all: PrepVarFrac PrepVarFracPolyCrystal

include ${PETSC_DIR}/conf/rules
include ${PETSC_DIR}/conf/variables

ifeq (${FC},gfortran)
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG} ${PETSC_ARCH} -I. -I../MEF90/${PETSC_ARCH} -I${M_VARFRACSTRUCTDIR} ${FC_FLAGS} --free-line-length-none
else
	MEF90_FC_FLAGS=${FC_MODULE_OUTPUT_FLAG} ${PETSC_ARCH} -I. -I../MEF90/${PETSC_ARCH} -I${M_VARFRACSTRUCTDIR} ${FC_FLAGS}
endif

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}

PrepVarFrac: ${PETSC_ARCH} PrepVarFrac.o ${PETSC_ARCH}/PrepVarFrac.o m_PrepVarFrac.o ${PETSC_ARCH}/m_PrepVarFrac.o
	${FLINKER} -o ${PETSC_ARCH}/PrepVarFrac ${PETSC_ARCH}/PrepVarFrac.o ${PETSC_ARCH}/m_PrepVarFrac.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}
	
PrepVarFracPolyCrystal: ${PETSC_ARCH} PrepVarFracPolyCrystal.o ${PETSC_ARCH}/PrepVarFracPolyCrystal.o m_PrepVarFrac.o ${PETSC_ARCH}/m_PrepVarFrac.o
	${FLINKER} -o ${PETSC_ARCH}/PrepVarFracPolyCrystal ${PETSC_ARCH}/PrepVarFracPolyCrystal.o ${PETSC_ARCH}/m_PrepVarFrac.o ${MEF90OBJS} ${M_VARFRACSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

PrepVarFrac.o: ${PETSC_ARCH} m_PrepVarFrac.o ${PETSC_ARCH}/m_PrepVarFrac.o
PrepVarFracPolyCrystal.o: ${PETSC_ARCH} m_PrepVarFrac.o ${PETSC_ARCH}/m_PrepVarFrac.o

m_PrepVarFrac.o: ${PETSC_ARCH}