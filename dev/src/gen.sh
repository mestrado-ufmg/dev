# Generate lib
gcc -O3 -c solver.c
gcc -shared -O3 -o ../bin/libsolver.so solver.o -lm -lblas -llapack -llapacke

# Remove intermediate file
rm solver.o