PROJ=MEF90
VER=0.12.0

MEF90DIR        = $(HOME)/eclipse/MEF90-$(VER)
MEF90OBJS	= $(MEF90DIR)/m_Constantes.o    $(MEF90DIR)/m_AlgebLin.o  \
                  $(MEF90DIR)/m_MEF_Types.o     $(MEF90DIR)/m_MEF_Tools.o \
                  $(MEF90DIR)/m_MEF_IO.o        $(MEF90DIR)/m_MEF_Gauss.o \
                  $(MEF90DIR)/m_MEF90.o         $(MEF90DIR)/m_Utils.o     \
                  $(MEF90DIR)/m_MEF_Ensight.o   $(MEF90DIR)/m_MEF_EXO.o   \
                  $(MEF90DIR)/m_MEF_SD.o        $(MEF90DIR)/m_MEF_PU.o    \
                  $(MEF90DIR)/m_MEF_Communication.o
EXO_CFLAGS      = -I$(MATSCI_DIR)/SEAMS.SRC/inc
EXO_LFLAGS      = -L$(MATSCI_DIR)/SEAMS.SRC/lib -lexoIIv2for -lexoIIv2c   \
                  -lnetcdf -lexoIIv2for -lsuplib -lsupes

all: m_MEF90.o

#include ../Makefile.${MEF90ARCH}
include ${PETSC_DIR}/bmake/common/base

FFLAGS = ${EXO_CFLAGS} -DUSE_PETSC_MAT -nochoicemod
# -CB -traceback -nothreads
#Dependencies:
m_Constantes:    
m_AlgebLin.o:    m_Constantes.o
m_Utils.o:       m_Constantes.o
m_MEF_Types.o:   m_AlgebLin.o
m_MEF_Gauss.o:   m_MEF_Types.o m_Utils.o
m_MEF_Tools.o:   m_MEF_Types.o
m_MEF_IO.o:      m_MEF_Types.o
m_MEF_Ensight.o: m_MEF_Types.o
m_MEF_EXO.o:     m_MEF_Types.o
m_MEF_Communication: m_MEF_Types.o m_MEF_EXO.o
m_MEF_PU.o:      m_MEF_Types.o m_MEF_EXO.o
m_MEF_SD.o:      m_MEF_Types.o m_MEF_EXO.o m_MEF_Communication.o
m_MEF90.o:       m_MEF_Types.o m_AlgebLin.o m_Constantes.o m_MEF_Ensight.o   \
                 m_MEF_PU.o m_MEF_Tools.o m_MEF_IO.o m_Utils.o m_MEF_Gauss.o \
                 m_MEF_EXO.o m_MEF_Communication.o m_MEF_SD.o

really_clean: clean
	- rm -f *~ core* *mod

