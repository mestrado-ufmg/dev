#include <stdio.h>

// void func(int a, double* b, double* result) {
//     for (int i = 0; i < a; i++) {
//         result[i] = b[i] + 1;
//     }
// }

void func2(int a, double** b, double** result) {
    for (int i = 0; i < a; i++) {
        result[i][0] = b[i][0] + 1;
        if (i == 0) {
            printf("%f", b[i][0]);
        }
    }
}

// void cfun(const void * indatav, int rowcount, int colcount, void * outdatav) {
//     //void cfun(const double * indata, int rowcount, int colcount, double * outdata) {
//     const double * indata = (double *) indatav;
//     double * outdata = (double *) outdatav;
//     int i;
//     puts("Here we go!");
//     for (i = 0; i < rowcount * colcount; ++i) {
//         outdata[i] = indata[i] * 2;
//     }
//     puts("Done!");
// }