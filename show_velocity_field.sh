#!/bin/sh
./Fluid $1 > test.gpl
gnuplot test.gpl
