#!/bin/sh
./Fluid $1 $2> test.gpl
gnuplot test.gpl
