#!/bin/bash

# shell script invoked with the following arguments
# $(CXX) $(CXX_DEFINES) $(CXX_FLAGS) -o OBJECT_FILE -c SOURCE_FILE

# extract parameters
SOURCE_FILE="${@: -1:1}"
OBJECT_FILE="${@: -3:1}"

# invoke compiler
TIME_BEFORE=`date +%s%N | cut -b1-13`
"$@"
TIME_AFTER=`date +%s%N | cut -b1-13`
let TIME=($TIME_AFTER-$TIME_BEFORE)
seconds=`echo $TIME/1000 | bc -l | cut -b1-6`
echo " *** Built object `basename \"$OBJECT_FILE\"` from $SOURCE_FILE in $seconds s."

 
