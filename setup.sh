#!/bin/sh
if [ "${ROOTSYS}" = "" ]; then
  echo "!!!! Set ROOT environment; Abort." 
else
  echo "ROOT is ready :-)"
  if [ "${DONE}" = "" ]; then
    export BASEDIR=${PWD}
    export CPATH=${BASEDIR}/include
    export LD_LIBRARY_PATH=${BASEDIR}/lib:${LD_LIBRARY_PATH}
    export DYLD_LIBRARY_PATH=${LD_LIBRARY_PATH}
    export DONE=1
  else
    echo "The environment is already set."
  fi
fi
