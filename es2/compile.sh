#!/bin/bash

rm -rf build
rm -rf install
mkdir build
mkdir install

CURRENT_PATH=$PWD

cd build

cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DCMAKE_INSTALL_PREFIX="$CURRENT_PATH/install" ..
make -j2
cd ..

