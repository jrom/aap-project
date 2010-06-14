#!/usr/bin/env bash

OS=$(uname)

if [ "$OS" = "Darwin" ]; then
  OSX_FLAGS="--build=powerpc-apple-bsd --host=powerpc-apple-bsd"
  C_FLAGS="-m32"
else
  OSX_FLAGS=""
  C_FLAGS=""
fi

cd ../sources/fftw-2.1.3/
mkdir installation
env CFLAGS="$C_FLAGS -O3 -march=core2" ./configure --prefix=$PWD/installation --enable-float $OSX_FLAGS
make
make install

