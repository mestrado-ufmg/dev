#ifndef STRUCTS_H
#define STRUCTS_H

#include <stdlib.h>

struct SurfaceMesh
{
    int nv;
    int nf;
    double *vertices;
    int *faces;
    double *facesAreas;
    double *facesMaxDistance;
    double *facesCenter;
    double *controlPoints;
    double *p1;
    double *p2;
    double *p3;
    double *e1;
    double *e2;
    double *e3;
};

struct WakeMeshPart
{
    int nSpan;
    int nWake;
    int *grid;
    double *vertices;
    int *faces;
};

struct WakeMesh
{
    struct WakeMeshPart left;
    struct WakeMeshPart right;
    struct WakeMeshPart tail;
};

struct Mesh
{
    struct SurfaceMesh surface;
    struct WakeMesh wake;
};

struct Environment {
    double vel_x, vel_y, vel_z;
    double velNorm;
    double density;
    double viscosity;
    double soundSpeed;
};

struct Input {
    int type;
    struct Mesh mesh;
    struct Environment environment;
};

struct Output {
    double *cp_v;
    double *vel_x_v;
    double *vel_y_v;
    double *vel_z_v;
    double *transpiration_v;
    double *sigma_v;
    double *doublet_v;
};

struct Point {
    double x;
    double y;
    double z;
};

struct VerticesConnection {
    int n;
    double *coeffs;
    int *faces;
};

struct FacesConnection {
    int n;
    int *faces;
};

struct VerticesConnection* getVerticesConnectionData(int nv) {
    struct VerticesConnection *verticesConnection;
    verticesConnection = (struct VerticesConnection*)malloc(nv * sizeof(struct VerticesConnection));
    return verticesConnection;
}

#endif