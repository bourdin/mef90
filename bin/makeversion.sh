#!/bin/bash
#echo Checking version against $1
GIT=`which git`
if [ $GIT ] 
then
	GITVERSION=`git describe --dirty --always --tags`
else
	GITVERSION=unknown
fi
if [ ! -e $1 ] 
then
	echo \#define MEF90_GITVER \"$GITVERSION\" > $1
else
 	MEF90GITVERSION=`cut -d \" -f 2 $1`
 	if [ ! $MEF90GITVERSION == $GITVERSION ]
 	then
 		echo updating $1 $MEF90GITVERSION $GITVERSION
		echo \#define MEF90_GITVER \"$GITVERSION\" > $1
	fi
fi
