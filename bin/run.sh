#!/usr/bin/env bash

function check_coord {
  cat $1 | tail -n 4 | head -3 | cut -b 40- | sort
}

function comp_coord {
  if [ "$(check_coord $1.out)"="$(check_coord $OUTPUTS/$1.output)" ];
  then
    echo "ok"
  else
    echo "error"
  fi
}

FTDOCK=../sources/3D_Dock/progs/ftdock
INPUTS=../inputs/proteins
OUTPUTS=../outputs

$FTDOCK -static $INPUTS/2pka.parsed -mobile $INPUTS/5pti.parsed > test1.out
$FTDOCK -static $INPUTS/1hba.parsed -mobile $INPUTS/5pti.parsed > test2.out
$FTDOCK -static $INPUTS/4hhb.parsed -mobile $INPUTS/5pti.parsed > test3.out

comp_coord test1
comp_coord test2
comp_coord test3

rm test*.out
rm *.dat
