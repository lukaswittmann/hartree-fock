#!/bin/bash

# make
echo " >> make"
make -j4

# run
echo " >> run"
./scf molecules/h2.in