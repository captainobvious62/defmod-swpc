#!/bin/bash

export PETSC_ARCH=linux-gnu-debug
export PETSC_DIR=/home/dockimble/petsc
export PATH=${PATH}:/home/dockimble/defmod-swpc/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${PETSC_DIR}/include
