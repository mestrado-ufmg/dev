#ifndef DATA_POTENTIAL_FLOW_H
#define DATA_POTENTIAL_FLOW_H

#include <stdlib.h>

struct PotentialFlowData
{
    double *sigma;
    double *doublet;
    int n;
    int na;
    double *a;
    int *ia;
    int *ja;
    double *rhs;
    double *a_vel_x;
    double *a_vel_y;
    double *a_vel_z;
    double *rhs_vel_x;
    double *rhs_vel_y;
    double *rhs_vel_z;
    double *cp;
    double *vel_x, *vel_y, *vel_z;
    double *transpiration;
};

struct PotentialFlowData getPotentialFlowData(int nf, double *e3, double vel_x, double vel_y, double vel_z) {

    struct PotentialFlowData data;

    data.sigma = (double*)malloc(nf * sizeof(double));
    data.doublet = (double*)malloc(nf * sizeof(double));
    data.n = nf;
    data.na = nf * nf;
    data.a = (double*)malloc(nf * nf * sizeof(double));
    data.ia = (int*)malloc(nf * nf * sizeof(int));
    data.ja = (int*)malloc(nf * nf * sizeof(int));
    data.rhs = (double*)malloc(nf * sizeof(double));
    data.a_vel_x = (double*)malloc(nf * nf * sizeof(double));
    data.a_vel_y = (double*)malloc(nf * nf * sizeof(double));
    data.a_vel_z = (double*)malloc(nf * nf * sizeof(double));
    data.rhs_vel_x = (double*)malloc(nf * sizeof(double));
    data.rhs_vel_y = (double*)malloc(nf * sizeof(double));
    data.rhs_vel_z = (double*)malloc(nf * sizeof(double));
    data.cp = (double*)malloc(nf * sizeof(double));
    data.vel_x = (double*)malloc(nf * sizeof(double));
    data.vel_y = (double*)malloc(nf * sizeof(double));
    data.vel_z = (double*)malloc(nf * sizeof(double));
    data.transpiration = (double*)malloc(nf * sizeof(double));

    int i, j;

    for (i = 0; i < nf; i++)
    {
        data.sigma[i] = -(e3[3 * i] * vel_x + e3[3 * i + 1] * vel_y + e3[3 * i + 2] * vel_z);
    }

    return data;
}

void freePotentialFlowData(struct PotentialFlowData data) {
    printf("> 0\n");
    free(data.sigma);
    printf("> 1\n");
    free(data.doublet);
    printf("> 2\n");
    free(data.a);
    printf("> 3\n");
    free(data.ia);
    printf("> 4\n");
    free(data.ja);
    printf("> 5\n");
    free(data.rhs);
    printf("> 6\n");
    free(data.a_vel_x);
    printf("> 7\n");
    free(data.a_vel_y);
    printf("> 8\n");
    free(data.a_vel_z);
    printf("> 9\n");
    free(data.rhs_vel_x);
    printf("> 10\n");
    free(data.rhs_vel_y);
    printf("> 11\n");
    free(data.rhs_vel_z);
    printf("> 12\n");
    free(data.cp);
    printf("> 13\n");
    free(data.vel_x);
    printf("> 14\n");
    free(data.vel_y);
    printf("> 15\n");
    free(data.vel_z);
    printf("> 16\n");
    free(data.transpiration);
    printf("> 17\n");
}

#endif