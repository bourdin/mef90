PROJ      = AmbrosioTortorelli
MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 
M_RUPTSTRUCTDIR  = ../m_RuptStruct/${PETSC_ARCH}
M_RUPTSTRUCTOBJS = ${M_RUPTSTRUCTDIR}/m_RuptStruct.o

all: AmbrosioTortorelli${DIM}

include ${PETSC_DIR}/conf/base

FFLAGS = -DPB_${DIM} -module ${PETSC_ARCH} -I. -I${MEF90DIR} -I${M_RUPTSTRUCTDIR}
# -ftrace=full -ffree-line-length-huge

.F90.o: 
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${FFLAGS} ${FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/${DIM}$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}
	
AmbrosioTortorelli${DIM}: ${PETSC_ARCH} AmbrosioTortorelli.o   ${PETSC_ARCH}/${DIM}AmbrosioTortorelli.o \
                                      m_AmbrosioTortorelli_Types.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_Types.o \
                                      m_AmbrosioTortorelli_U.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_U.o \
                                      m_AmbrosioTortorelli_V.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_V.o \
                                      m_AmbrosioTortorelli.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli.o \
                                      m_AmbrosioTortorelli_Post.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_Post.o
	${FLINKER} -o ${PETSC_ARCH}/AmbrosioTortorelli${DIM} ${PETSC_ARCH}/${DIM}AmbrosioTortorelli.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_Types.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_U.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_V.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli.o  ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_Post.o ${MEF90OBJS} ${M_RUPTSTRUCTOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
	
AmbrosioTortorelli.o:   ${PETSC_ARCH} m_AmbrosioTortorelli.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli.o \
                                      m_AmbrosioTortorelli_Types.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_Types.o \
                                      m_AmbrosioTortorelli_U.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_U.o \
                                      m_AmbrosioTortorelli_V.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_V.o \
                                      m_AmbrosioTortorelli_Post.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_Post.o

m_AmbrosioTortorelli.o: ${PETSC_ARCH} m_AmbrosioTortorelli_Types.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_Types.o \
                                      m_AmbrosioTortorelli_U.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_U.o \
                                      m_AmbrosioTortorelli_V.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_V.o \
                                      m_AmbrosioTortorelli_Post.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_Post.o
m_AmbrosioTortorelli_U.o:  ${PETSC_ARCH} m_AmbrosioTortorelli_Types.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_Types.o
m_AmbrosioTortorelli_V.o:  ${PETSC_ARCH} m_AmbrosioTortorelli_Types.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_Types.o
m_AmbrosioTortorelli_Post.o:  ${PETSC_ARCH} m_AmbrosioTortorelli_Types.o ${PETSC_ARCH}/${DIM}m_AmbrosioTortorelli_Types.o
m_AmbrosioTortorelli_Types.o: ${PETSC_ARCH}


clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}
