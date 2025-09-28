#!/bin/bash

file=$1
./easymesh $file
./fem_program cuadrado $file
./contour_plot $file.n $file.e &
./mesh_plot $file &
wait
