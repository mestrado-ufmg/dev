gcc -std=c11 -Wall -Wextra -pedantic -c -fPIC lib.c -o lib.o
gcc -shared lib.o -o lib.so