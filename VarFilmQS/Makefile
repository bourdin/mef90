PROJ      = VarFilmQS
OBJDIR_PROJ = ${OBJDIR}/${PROJ}

M_FILMSTRUCTDIR  = ${OBJDIR}/m_Film_Struct
M_FILMSTRUCTOBJS = ${M_FILMSTRUCTDIR}/m_Film_Struct.o

all: VarFilmQS

include ../Makefile.inc 
include ../Makefile.mk

include ${TAO_DIR}/bmake/tao_common

ifneq (,$(findstring gcc,$(PETSC_ARCH)))
	MEF90_FC_FLAGS=-DWITH_TAO ${FC_MODULE_OUTPUT_FLAG}${OBJDIR_PROJ} -I${OBJDIR_PROJ} -I${MEF90DIR} -I${M_FILMSTRUCTDIR} ${FC_FLAGS} -ffree-line-length-none -ffixed-line-length-none  -ffree-form 	
else
	MEF90_FC_FLAGS=-DWITH_TAO ${FC_MODULE_OUTPUT_FLAG} ${OBJDIR_PROJ} -I${OBJDIR_PROJ} -I${MEF90DIR} -I${M_FILMSTRUCTDIR} ${FC_FLAGS}
endif

.F90.o:
	${PETSC_MAKE_STOP_ON_ERROR}${FC} -c ${MEF90_FC_FLAGS} ${FCPPFLAGS} -o ${OBJDIR_PROJ}/$@ $<

dir:
	@mkdir -p ${BINDIR}
	@mkdir -p ${OBJDIR_PROJ}

VarFilmQS: dir VarFilmQS.o   ${OBJDIR_PROJ}/VarFilmQS.o \
                                      m_VarFilmQS_Types.o ${OBJDIR_PROJ}/m_VarFilmQS_Types.o \
                                      m_VarFilmQS_U.o ${OBJDIR_PROJ}/m_VarFilmQS_U.o \
                                      m_VarFilmQS_V.o ${OBJDIR_PROJ}/m_VarFilmQS_V.o \
                                      m_VarFilmQS_W.o ${OBJDIR_PROJ}/m_VarFilmQS_W.o \
                                      m_VarFilmQS.o ${OBJDIR_PROJ}/m_VarFilmQS.o \
                                      m_VarFilmQS_Post.o ${OBJDIR_PROJ}/m_VarFilmQS_Post.o
				${FLINKER} -o ${BINDIR}/VarFilmQS ${OBJDIR_PROJ}/VarFilmQS.o ${OBJDIR_PROJ}/m_VarFilmQS_Types.o ${OBJDIR_PROJ}/m_VarFilmQS_U.o ${OBJDIR_PROJ}/m_VarFilmQS_V.o ${OBJDIR_PROJ}/m_VarFilmQS_W.o ${OBJDIR_PROJ}/m_VarFilmQS.o  ${OBJDIR_PROJ}/m_VarFilmQS_Post.o ${MEF90OBJS} ${M_FILMSTRUCTOBJS} ${TAO_FORTRAN_LIB} ${TAO_LIB} ${PETSC_FORTRAN_LIB} ${PETSC_LIB}

VarFilmQS.o:  dir m_VarFilmQS.o ${OBJDIR_PROJ}/m_VarFilmQS.o \
                                      m_VarFilmQS_Types.o ${OBJDIR_PROJ}/m_VarFilmQS_Types.o \
                                      m_VarFilmQS_U.o ${OBJDIR_PROJ}/m_VarFilmQS_U.o \
                                      m_VarFilmQS_V.o ${OBJDIR_PROJ}/m_VarFilmQS_V.o \
                                      m_VarFilmQS_W.o ${OBJDIR_PROJ}/m_VarFilmQS_W.o \
                                      m_VarFilmQS_Post.o ${OBJDIR_PROJ}/m_VarFilmQS_Post.o


m_VarFilmQS.o: dir m_VarFilmQS_Types.o ${OBJDIR_PROJ}/m_VarFilmQS_Types.o \
                                      m_VarFilmQS_U.o ${OBJDIR_PROJ}/m_VarFilmQS_U.o \
                                      m_VarFilmQS_V.o ${OBJDIR_PROJ}/m_VarFilmQS_V.o \
                                      m_VarFilmQS_W.o ${OBJDIR_PROJ}/m_VarFilmQS_W.o \
                                      m_VarFilmQS_Post.o ${OBJDIR_PROJ}/m_VarFilmQS_Post.o
m_VarFilmQS_U.o:  dir m_VarFilmQS_Types.o ${OBJDIR_PROJ}/m_VarFilmQS_Types.o \
                                m_VarFilmQS_Post.o ${OBJDIR_PROJ}/m_VarFilmQS_Post.o \
                                m_VarFilmQS_W.o ${OBJDIR_PROJ}/m_VarFilmQS_W.o 
m_VarFilmQS_V.o:  dir m_VarFilmQS_Types.o ${OBJDIR_PROJ}/m_VarFilmQS_Types.o \
                                m_VarFilmQS_Post.o ${OBJDIR_PROJ}/m_VarFilmQS_Post.o
m_VarFilmQS_W.o:  dir m_VarFilmQS_Types.o ${OBJDIR_PROJ}/m_VarFilmQS_Types.o \
                                m_VarFilmQS_Post.o ${OBJDIR_PROJ}/m_VarFilmQS_Post.o
m_VarFilmQS_Post.o:  dir m_VarFilmQS_Types.o ${OBJDIR_PROJ}/m_VarFilmQS_Types.o
m_VarFilmQS_Types.o: dir

