#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    int i;
    int n;
    char *p;

    long conv = strtol(argv[1], &p, 10);
    n = conv;

    int sum = 0;
    for (i = 0; i < n; i++) {
        sum += i;
    }
    printf("%d", sum);
    return 0;
}