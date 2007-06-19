PROJ=Prep_Poisson_Constant
MEF90DIR        = ../MEF90
MEF90OBJS       = $(MEF90DIR)/m_Constantes.o    $(MEF90DIR)/m_AlgebLin.o  \
                  $(MEF90DIR)/m_MEF_Types.o     $(MEF90DIR)/m_MEF_Gauss.o \
                  $(MEF90DIR)/m_MEF90.o         $(MEF90DIR)/m_Utils.o     \
                  $(MEF90DIR)/m_MEF_Ensight.o   $(MEF90DIR)/m_MEF_EXO.o   \
                  $(MEF90DIR)/m_MEF_SD.o        $(MEF90DIR)/m_MEF_PU.o    \
                  $(MEF90DIR)/m_MEF_Communication.o


M_POISSON_STRUCT_DIR=../m_Poisson_Struct
M_POISSON_STRUCT_OBJS=${M_POISSON_STRUCT_DIR}/m_Poisson.o

EXO_CFLAGS = -I$(HPCDIST_DIR)/include
EXO_LFLAGS = -L$(HPCDIST_DIR)/lib -lexoIIv2for -lexoIIv2c -lnetcdf -lmetis -lparmetis

all: ${PROJ}

include ${PETSC_DIR}/bmake/common/base

OBJS=    $(PROJ).o
FFLAGS = -I${MEF90DIR} -I${M_POISSON_STRUCT_DIR} ${EXO_CFLAGS}

${PROJ}: ${PROJ}.o
	-${FLINKER} -o ${PROJ} ${OBJS} ${MEF90OBJS} ${M_POISSON_STRUCT_OBJS} ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${EXO_LFLAGS}

# Dependences
$(PROJ).o: $(MEF90OBJS) 

really_clean: clean
	- rm -f ${PROJ}-${VER} *~ core

tarball: really_clean
	tar zcvf ${PROJ}-${VER}-`date '+%Y-%m-%d'`-`hostname -s`.tgz *[^tgz] > /dev/null

