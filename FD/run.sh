#!/bin/bash

make clean 2> /dev/null # this will terminate the run if the make fails
set -e
make
# run the program and redirect the error output
./exe_FD_hMetal 2> a.err
# run the program and redirect screen and error output
# ./exe_FD_hMetal > a.out 2> a.err
