PROJ=Prep_Rupt
VER=0.0.7
MEF90VER        = 0.12.0
MEF90DIR        = $(MEF90ROOT)/MEF90-$(MEF90VER)
MEF90OBJS       = $(MEF90DIR)/m_Constantes.o    $(MEF90DIR)/m_AlgebLin.o  \
                  $(MEF90DIR)/m_MEF_Types.o     $(MEF90DIR)/m_MEF_Tools.o \
                  $(MEF90DIR)/m_MEF_IO.o        $(MEF90DIR)/m_MEF_Gauss.o \
                  $(MEF90DIR)/m_MEF90.o         $(MEF90DIR)/m_Utils.o     \
                  $(MEF90DIR)/m_MEF_Ensight.o   $(MEF90DIR)/m_MEF_EXO.o   \
                  $(MEF90DIR)/m_MEF_SD.o        $(MEF90DIR)/m_MEF_PU.o    \
                  $(MEF90DIR)/m_MEF_Communication.o

M_Rupt_STRUCT_DIR = $(MEF90ROOT)/m_Rupt_Struct-0.0.6
M_Rupt_STRUCT_OBJS = $(M_Rupt_STRUCT_DIR)/m_Rupt_Struct.o

EXO_CFLAGS      = -I$(MATSCI_DIR)/SEAMS.SRC/inc
EXO_LFLAGS      = -L$(MATSCI_DIR)/SEAMS.SRC/lib -lexoIIv2for -lexoIIv2c   \
                  -lnetcdf -lexoIIv2for -lsuplib -lsupes

all: ${PROJ}-${VER}
include ${PETSC_DIR}/bmake/common/base

OBJS=    $(PROJ).o
FFLAGS = -I${MEF90DIR} -I${M_Rupt_STRUCT_DIR} ${EXO_CFLAGS}

$(PROJ)-${VER}: ${PROJ}.o
	${FLINKER} -o ${PROJ}-${VER} ${OBJS} ${MEF90OBJS}                     \
                    ${M_Rupt_STRUCT_OBJS}                                     \
                    ${PETSC_FORTRAN_LIB} ${PETSC_LIB} ${EXO_LFLAGS}


# Dependences
$(PROJ).o: $(MEF90OBJS) 

really_clean: clean
	- rm -f ${PROJ}-${VER} *~ core

tarball: really_clean
	tar zcvf ${PROJ}-${VER}-`date '+%Y-%m-%d'`-`hostname -s`.tgz *[^tgz] > /dev/null

