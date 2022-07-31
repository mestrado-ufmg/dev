gcc -std=c11 -Wall -O3 -Wextra -pedantic -c -fPIC lib.c -o source.o
gcc -shared -O3 source.o -o source.so