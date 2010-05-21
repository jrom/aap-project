#!/usr/bin/env bash

cd ../sources/
rm -rf fftw-2.1.3
git checkout fftw-2.1.3
cd 3D_Dock/progs/
make clean

