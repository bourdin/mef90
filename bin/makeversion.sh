#!/bin/bash
echo TESTING FOR $1
if [ ! -e $1 ] 
  then echo \#define MEF90_HGVER \"unknown\" >> $1
else
  HG=`which hg`
  if [ $HG ]
    then echo \#define MEF90_HGVER \"`hg parents | head -1 | cut -d : -f 2,3 | tr -d ' '`\" > $1
  fi
fi
