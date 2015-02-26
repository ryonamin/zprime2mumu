#!/bin/sh
if [ "${ROOTSYS}" = "" ]; then
  echo "!!!! Set ROOT environment; Abort." 
else
  echo "ROOT is ready :-)"
  if [ "${BASEDIR}" != "${PWD}" ]; then
    export BASEDIR=${PWD}
    export CPATH=${BASEDIR}/include
    export LD_LIBRARY_PATH=${BASEDIR}/lib:${LD_LIBRARY_PATH}
    export DYLD_LIBRARY_PATH=${LD_LIBRARY_PATH}
    echo "The environment has been set."
  else
    echo "The environment is already set."
  fi
fi
