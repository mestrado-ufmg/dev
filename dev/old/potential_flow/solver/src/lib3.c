#include <stdio.h>
#include <stdlib.h>

void convert(int lines, int rows, double *in, double out[lines][rows]) {

    int i, j, k;

    j = 0;
    k = 0;

    for (i = 0; i < lines * rows; i++) {

        out[j][k] = in[i];

        k++;
        if (k == rows) {
            j++;
            k = 0;
        }
    }

}

void setValues(int lines, int rows, double in[lines][rows], double *out) {

    int i, j, k;

    j = 0;
    k = 0;

    for (i = 0; i < lines * rows; i++) {

        out[i] = in[j][k];

        k++;
        if (k == rows) {
            j++;
            k = 0;
        }
    }

}

void func1(int lines, int rows, double *in, double *out) {

    double inMatrix[lines][rows];
    double outMatrix[lines][rows];

    convert(lines, rows, in, inMatrix);

    int i, j;

    for (i = 0; i < lines; i++) {
        for (j = 0; j < rows; j++) {
            outMatrix[i][j] = inMatrix[i][j] + 1;
        }
    }

    setValues(lines, rows, outMatrix, out);

}

void func2(int lines, int rows, double **in, double **out) {

    int i, j;

    for (i = 0; i < lines; i++) {
        for (j = 0; j < rows; j++) {
            out[i][j] = in[i][j] + 1;
        }
    }

}

void func3(int lines, int rows, double *in, double *out) {

    int i, j;

    for (i = 0; i < lines; i++) {
        for (j = 0; j < rows; j++) {
            out[i * rows + j] = in[i * rows + j] + 1;
        }
    }

}

void solve_system(int n, double **a, double *b, double *sol) {

    int i, j, k;

    int max_index = n - 1;

    double mult;
    double sum;

    double **aux_a = (double**)malloc(n * sizeof(double*));
    for (i = 0; i < n; i++) {
        aux_a[i] = (double*)malloc(n * sizeof(double));
    }

    double *aux_b = (double*)malloc(n * sizeof(double*));

    // Insert first row on aux_a and aux_b
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            aux_a[i][j] = a[i][j];
        }
        aux_b[i] = b[i];
    }

    // Gauss elimination
    for (k = 0; k < n - 1; k++) {
        for (i = k + 1; i < n; i++) {
            mult = - aux_a[i][k] / aux_a[k][k];
            for (j = k; j < n; j++) {
                aux_a[i][j] = mult * aux_a[k][j] + aux_a[i][j];
            }
            aux_b[i] = mult * aux_b[k] + aux_b[i];
        }
    }

    // Solve system
    sol[max_index] = aux_b[max_index] / aux_a[max_index][max_index];
    for (i = n - 2; i >= 0; i--) {
        sum = aux_b[i];
        for (j = i + 1; j < n; j++) {
            sum -= aux_a[i][j] * sol[j];
        }
        sol[i] = sum / aux_a[i][i];
    }
}

void solve_system_1D(int n, double *a, double *b, double *sol) {

    int i, j, k;

    int max_index = n - 1;

    double mult;
    double sum;

    double *aux_a = (double*)malloc(n * n * sizeof(double*));
    double *aux_b = (double*)malloc(n * sizeof(double*));

    // Insert first row on aux_a and aux_b
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            aux_a[i * n + j] = a[i * n + j];
        }
        aux_b[i] = b[i];
    }

    // Gauss elimination
    for (k = 0; k < n - 1; k++) {
        for (i = k + 1; i < n; i++) {
            mult = - aux_a[i * n + k] / aux_a[k * n + k];
            for (j = k; j < n; j++) {
                aux_a[i * n + j] = mult * aux_a[k * n + j] + aux_a[i * n + j];
            }
            aux_b[i] = mult * aux_b[k] + aux_b[i];
        }
    }

    // Solve system
    sol[max_index] = aux_b[max_index] / aux_a[max_index * n + max_index];
    for (i = n - 2; i >= 0; i--) {
        sum = aux_b[i];
        for (j = i + 1; j < n; j++) {
            sum -= aux_a[i * n + j] * sol[j];
        }
        sol[i] = sum / aux_a[i * n + i];
    }
}