install lapack

install c/c++ core

sudo apt update

sudo apt install build-essential

install libraries

$ sudo apt install liblapack3
$ sudo apt install liblapack-dev
$ sudo apt install libopenblas-base
$ sudo apt install libopenblas-dev
$ sudo apt install liblapacke-dev

Install static linking libraries if needed

$ sudo apt install liblapack-dev

Add the necessary linking flags to gcc/g++

-lm -lblas -llapack -llapacke