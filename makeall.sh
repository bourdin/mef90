#!/bin/sh
make -C MEF90
make -C m_VarFrac_Struct
make -C PrepVarFrac
make -C PrepVarFracNG
make -C VarFracQS DIM=2D
make -C VarFracQS DIM=3D
