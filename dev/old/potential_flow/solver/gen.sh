gcc -std=c11 -Wall -O3 -Wextra -pedantic -c -fPIC ./src/lib3.c -o ./src/lib3.o
gcc -shared -O3 ./src/lib3.o -o lib3.so