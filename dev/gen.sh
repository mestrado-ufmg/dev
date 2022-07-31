gcc -std=c11 -Wall -O3 -Wextra -pedantic -c -fPIC ./src/solver_lib.c -o ./src/solver_lib.o
gcc -shared -O3 ./src/solver_lib.o -o solver_lib.so