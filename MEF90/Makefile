PROJ=MEF90

MEF90OBJS	=     m_Constantes.o    m_AlgebLin.o  \
                  m_MEF_Types.o     m_MEF_Gauss.o \
                  m_MEF90.o         m_Utils.o     \
                  m_MEF_Ensight.o   m_MEF_EXO.o   \
                  m_MEF_SD.o        m_MEF_PU.o    \
                  m_MEF_Communication.o
EXO_CFLAGS      = -I$(HPCDIST_DIR)/include

all: m_MEF90.o

include ${PETSC_DIR}/bmake/common/base

FFLAGS = ${EXO_CFLAGS}
# -DUSE_PETSC_MAT -nochoicemod

#Dependencies:
m_Constantes:    
m_AlgebLin.o:    m_Constantes.o
m_Utils.o:       m_Constantes.o
m_MEF_Types.o:   m_AlgebLin.o
m_MEF_Gauss.o:   m_MEF_Types.o m_Utils.o
m_MEF_Ensight.o: m_MEF_Types.o
m_MEF_EXO.o:     m_MEF_Types.o
m_MEF_Communication: m_MEF_Types.o m_MEF_EXO.o
m_MEF_PU.o:      m_MEF_Types.o m_MEF_EXO.o
m_MEF_SD.o:      m_MEF_Types.o m_MEF_EXO.o m_MEF_Communication.o
m_MEF90.o:       m_MEF_Types.o m_AlgebLin.o m_Constantes.o m_MEF_Ensight.o   \
                 m_MEF_PU.o m_Utils.o m_MEF_Gauss.o \
                 m_MEF_EXO.o m_MEF_Communication.o m_MEF_SD.o

really_clean: clean
	- rm -f *~ core* *mod

