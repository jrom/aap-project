#!/usr/bin/env bash

OS=$(uname)

if [ "$OS" = "Darwin" ]; then
  OSX_FLAGS="--build=powerpc-apple-bsd --host=powerpc-apple-bsd"
else
  OSX_FLAGS=""
fi

cd ../sources/fftw-2.1.3/
mkdir installation
env CFLAGS="-O3 -march=core2" ./configure --prefix=$PWD/installation --enable-float $OSX_FLAGS
make
make install

