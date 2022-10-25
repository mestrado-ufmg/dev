# Generate lib
gcc -fPIC -Ofast -c solver.c
gcc -shared -Ofast -o libsolver.so solver.o -lm -lrt -lblas -llapack -llapacke

# Remove intermediate file
rm solver.o