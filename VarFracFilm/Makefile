PROJ      = VarFracFilm
MEF90DIR  = ../MEF90/${PETSC_ARCH}
MEF90OBJS = $(MEF90DIR)/m_MEF_Parameters.o $(MEF90DIR)/m_MEF_LinAlg.o    \
            $(MEF90DIR)/m_MEF_Types.o      $(MEF90DIR)/m_MEF_Elements.o  \
            $(MEF90DIR)/m_MEF90.o          $(MEF90DIR)/m_MEF_Utils.o         \
            $(MEF90DIR)/m_MEF_EXO.o        $(MEF90DIR)/m_MEF_MPI.o       \
            $(MEF90DIR)/m_MEF_Sieve.o 
M_VARFRACFILMDIR  = ../m_VarFracFilm_Struct/${PETSC_ARCH}
M_VARFRACFILMOBJS = ${M_VARFRACFILMDIR}/m_VarFracFilm_Struct.o

all: VarFracFilm${DIM}

include ${PETSC_DIR}/conf/base

FFLAGS = -DPB_${DIM} -module ${PETSC_ARCH} -I. -I${MEF90DIR} -I${M_VARFRACFILMDIR}
# -ftrace=full -ffree-line-length-huge

.F90.o: 
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${FFLAGS} ${FC_FLAGS} ${FCPPFLAGS} -o ${PETSC_ARCH}/${DIM}$@ $<

${PETSC_ARCH}:
	@mkdir ${PETSC_ARCH}
	
VarFracFilm${DIM}: ${PETSC_ARCH} VarFracFilm.o   ${PETSC_ARCH}/${DIM}VarFracFilm.o \
                                      m_VarFracFilm_Types.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Types.o \
                                      m_VarFracFilm_U.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_U.o \
                                      m_VarFracFilm_V.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_V.o \
                                      m_VarFracFilm_Phi.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Phi.o \
                                      m_VarFracFilm.o ${PETSC_ARCH}/${DIM}m_VarFracFilm.o \
                                      m_VarFracFilm_Post.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Post.o
	${FLINKER} -o ${PETSC_ARCH}/VarFracFilm${DIM} ${PETSC_ARCH}/${DIM}VarFracFilm.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Types.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_U.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_V.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Phi.o ${PETSC_ARCH}/${DIM}m_VarFracFilm.o  ${PETSC_ARCH}/${DIM}m_VarFracFilm_Post.o ${MEF90OBJS} ${M_VARFRACFILMOBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} 
	
VarFracFilm.o:   ${PETSC_ARCH} m_VarFracFilm.o ${PETSC_ARCH}/${DIM}m_VarFracFilm.o \
                                      m_VarFracFilm_Types.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Types.o \
                                      m_VarFracFilm_U.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_U.o \
                                      m_VarFracFilm_V.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_V.o \
                                      m_VarFracFilm_Phi.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Phi.o \
                                      m_VarFracFilm_Post.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Post.o

m_VarFracFilm.o: ${PETSC_ARCH} m_VarFracFilm_Types.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Types.o \
                                      m_VarFracFilm_U.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_U.o \
                                      m_VarFracFilm_V.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_V.o \
                                      m_VarFracFilm_Phi.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Phi.o \
                                      m_VarFracFilm_Post.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Post.o
m_VarFracFilm_U.o:  ${PETSC_ARCH} m_VarFracFilm_Types.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Types.o
m_VarFracFilm_V.o:  ${PETSC_ARCH} m_VarFracFilm_Types.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Types.o
m_VarFracFilm_Post.o:  ${PETSC_ARCH} m_VarFracFilm_Types.o ${PETSC_ARCH}/${DIM}m_VarFracFilm_Types.o
m_VarFracFilm_Types.o: ${PETSC_ARCH}


clean::
	@rm -Rf ${PETSC_ARCH}/*o ${PETSC_ARCH}/*mod
	@rmdir ${PETSC_ARCH}
