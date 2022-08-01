#include "solver_lib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

///----------------------------------------------------------------///
/// Helpers
///----------------------------------------------------------------///
const double ZERO_ERROR = 1e-8;
const double PI = 3.14159265359;
const double FACTOR = 1 / (4 * PI);

struct Point {
	double x;
    double y;
    double z;
};

double division(double a,
                double b) {
    if ((-ZERO_ERROR < b) && (b < 0)) {
        return - a / ZERO_ERROR;
    } else if ((0 < b) && (b < ZERO_ERROR)) {
        return a / ZERO_ERROR;
    } else {
        return a / b;
    }
}

double norm(struct Point p) {
    return sqrt(pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2));
}

struct Point cross(struct Point p1,
                   struct Point p2) {
    struct Point p = {p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z, p1.x * p2.y - p1.y * p2.x};
    return p;
}

///----------------------------------------------------------------///
/// Potential flow
///----------------------------------------------------------------///
void sourceFunc(struct Point p,
                struct Point p1,
                struct Point p2,
                struct Point p3,
                struct Point e1,
                struct Point e2,
                struct Point e3,
                double area,
                double maxDistance,
                double* vel) {

    double u, v, w;
    double distance = norm(p);

    if (distance > maxDistance) {

        double pNorm3 = pow(norm(p), 3);

        u = FACTOR * area * p.x / pNorm3;
        v = FACTOR * area * p.y / pNorm3;
        w = FACTOR * area * p.z / pNorm3;

    } else {

        // Velocity parameters
        double r1, r2, r3;
        double l1, l2, l3;
        double h1, h2, h3;

        double d12, d23, d31;
        double m12, m23, m31;
        double ln12, ln23, ln31;

        // Calculate
        r1 = sqrt(pow(p.x - p1.x, 2) + pow(p.y - p1.y, 2) + pow(p.z, 2));
        r2 = sqrt(pow(p.x - p2.x, 2) + pow(p.y - p2.y, 2) + pow(p.z, 2));
        r3 = sqrt(pow(p.x - p3.x, 2) + pow(p.y - p3.y, 2) + pow(p.z, 2));

        l1 = pow(p.x - p1.x, 2) + pow(p.z, 2);
        l2 = pow(p.x - p2.x, 2) + pow(p.z, 2);
        l3 = pow(p.x - p3.x, 2) + pow(p.z, 2);

        h1 = (p.x - p1.x) * (p.y - p1.y);
        h2 = (p.x - p2.x) * (p.y - p2.y);
        h3 = (p.x - p3.x) * (p.y - p3.y);

        d12 = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
        m12 = division(p2.y - p1.y, p2.x - p1.x);

        d23 = sqrt(pow(p3.x - p2.x, 2) + pow(p3.y - p2.y, 2));
        m23 = division(p3.y - p2.y, p3.x - p2.x);

        d31 = sqrt(pow(p1.x - p3.x, 2) + pow(p1.y - p3.y, 2));
        m31 = division(p1.y - p3.y, p1.x - p3.x);

        ln12 = log(division(r1 + r2 - d12, r1 + r2 + d12));
        ln23 = log(division(r2 + r3 - d23, r2 + r3 + d23));
        ln31 = log(division(r3 + r1 - d31, r3 + r1 + d31));

        u = -FACTOR * ( division((p2.y - p1.y), d12) * ln12 + division((p3.y - p2.y), d23) * ln23 + division((p1.y - p3.y), d31) * ln31 );
        v = FACTOR * ( division((p2.x - p1.x), d12) * ln12 + division((p3.x - p2.x), d23) * ln23 + division((p1.x - p3.x), d31) * ln31 );
        w = -FACTOR * ( atan(division(m12 * l1 - h1, p.z * r1)) - atan(division(m12 * l2 - h2, p.z * r2)) + atan(division(m23 * l2 - h2, p.z * r2)) - atan(division(m23 * l3 - h3, p.z * r3)) + atan(division(m31 * l3 - h3, p.z * r3)) - atan(division(m31 * l1 - h1, p.z * r1)) );

    }

    vel[0] = u * e1.x + v * e2.x + w * e3.x;
    vel[1] = u * e1.y + v * e2.y + w * e3.y;
    vel[2] = u * e1.z + v * e2.z + w * e3.z;
}

void lineFunc(struct Point p,
              struct Point p1,
              struct Point p2,
              double* vel) {

    struct Point r1 = {p1.x - p.x, p1.y - p.y, p1.z - p.z};
    struct Point r2 = {p2.x - p.x, p2.y - p.y, p2.z - p.z};

    struct Point r1xr2 = cross(r1, r2);

    double r1Norm = norm(r1);
    double r2Norm = norm(r2);

    double r1xr2Norm2 = pow(norm(r1xr2), 2);

    double dot = (1 / r1xr2Norm2) * ((r1.x - r2.x) * (r1.x / r1Norm - r2.x / r2Norm) + (r1.y - r2.y) * (r1.y / r1Norm - r2.y / r2Norm) + (r1.z - r2.z) * (r1.z / r1Norm - r2.z / r2Norm));

    vel[0] = FACTOR * r1xr2.x * dot;
    vel[1] = FACTOR * r1xr2.y * dot;
    vel[2] = FACTOR * r1xr2.z * dot;
}

void doubletFunc(struct Point p,
                 struct Point p1,
                 struct Point p2,
                 struct Point p3,
                 struct Point e1,
                 struct Point e2,
                 struct Point e3,
                 double area,
                 double maxDistance,
                 double* vel) {

    double distance = norm(p);

    if (distance > maxDistance) {

        double pxLocal = p.x * e1.x + p.y * e1.y + p.z * e1.z;
        double pyLocal = p.x * e2.x + p.y * e2.y + p.z * e2.z;
        double pzLocal = p.x * e3.x + p.y * e3.y + p.z * e3.z;
        double den = pow(pxLocal * pxLocal + pyLocal * pyLocal + pzLocal * pzLocal, 2.5);

        double u = 0.75 * FACTOR * area * pzLocal * pxLocal / den;
        double v = 0.75 * FACTOR * area * pzLocal * pyLocal / den;
        double w = - FACTOR * area * (pxLocal * pxLocal + pyLocal * pyLocal - 2 * pzLocal * pzLocal) / den;

        vel[0] = u * e1.x + v * e2.x + w * e3.x;
        vel[1] = u * e1.y + v * e2.y + w * e3.y;
        vel[2] = u * e1.z + v * e2.z + w * e3.z;

    } else {

        double u, v, w;

        double *vel1 = (double*)malloc(3 * sizeof(double*));
        double *vel2 = (double*)malloc(3 * sizeof(double*));
        double *vel3 = (double*)malloc(3 * sizeof(double*));

        lineFunc(p, p1, p2, vel1);
        lineFunc(p, p2, p3, vel2);
        lineFunc(p, p3, p1, vel3);

        u = vel1[0] + vel2[0] + vel3[0];
        v = vel1[1] + vel2[1] + vel3[1];
        w = vel1[2] + vel2[2] + vel3[2];

        vel[0] = u * e1.x + v * e2.x + w * e3.x;
        vel[1] = u * e1.y + v * e2.y + w * e3.y;
        vel[2] = u * e1.z + v * e2.z + w * e3.z;

        free(vel1);
        free(vel2);
        free(vel3);

    }
}

///----------------------------------------------------------------///
/// Solves the potential and boundary layer correction
///----------------------------------------------------------------///
void solve(int n,
           double *facesAreas,
           double *facesMaxDistance,
           double *facesCenter,
           double *controlPoints,
           double *p1,
           double *p2,
           double *p3,
           double *e1,
           double *e2,
           double *e3,
           double *freestream,
           double *sigma,
           int nSpanLeftWing,
           int nWakeLeftWing,
           int *leftWingGrid,
           double *leftWingVertices,
           int *leftWingFaces,
           int nSpanRightWing,
           int nWakeRightWing,
           int *rightWingGrid,
           double *rightWingVertices,
           int *rightWingFaces,
           int nSpanTail,
           int nWakeTail,
           int *tailGrid,
           double *tailVertices,
           int *tailFaces,
           double *matrix,
           double *array,
           double *matrixVelx,
           double *matrixVely,
           double *matrixVelz,
           double *arrayVel) {
    
    ///---------------------------------------///
    /// Parameters
    
    // Loops
    int i, j, k, l;

    // Point
    int i3D1, i3D2, i3D3, j3D1, j3D2, j3D3;
    int j2D1, j2D2;
    int indexLine1, indexLine2;
    struct Point p;
    struct Point p1Line;
    struct Point p2Line;
    struct Point pLocal;
    struct Point p1Local; p1Local.z = 0.0;
    struct Point p2Local; p2Local.z = 0.0;
    struct Point p3Local; p3Local.z = 0.0;

    // Base vectors
    struct Point e3iPoint;

    struct Point e1jPoint;
    struct Point e2jPoint;
    struct Point e3jPoint;

    // Velocities
    double *sourceVel = (double*)malloc(3 * sizeof(double));
    double *doubletVel = (double*)malloc(3 * sizeof(double));
    double *lineVel = (double*)malloc(3 * sizeof(double));
    
    ///---------------------------------------///
    /// Linear system

    // Create
    for (i = 0; i < n; i++) {

        i3D1 = i * 3;
        i3D2 = i3D1 + 1;
        i3D3 = i3D1 + 2;

        e3iPoint.x = e3[i3D1]; e3iPoint.y = e3[i3D2]; e3iPoint.z = e3[i3D3];

        array[i] = 0;
        arrayVel[i * 3] = 0;
        arrayVel[i * 3 + 1] = 0;
        arrayVel[i * 3 + 2] = 0;
        
        // Surface
        // Effect of j on i
        for (j = 0; j < n; j++) {

            j3D1 = j * 3;
            j3D2 = j3D1 + 1;
            j3D3 = j3D1 + 2;

            j2D1 = j * 2;
            j2D2 = j2D1 + 1;

            // Points
            e1jPoint.x = e1[j3D1]; e1jPoint.y = e1[j3D2]; e1jPoint.z = e1[j3D3];
            e2jPoint.x = e2[j3D1]; e2jPoint.y = e2[j3D2]; e2jPoint.z = e2[j3D3];
            e3jPoint.x = e3[j3D1]; e3jPoint.y = e3[j3D2]; e3jPoint.z = e3[j3D3];

            p.x = controlPoints[i3D1] - facesCenter[j3D1];
            p.y = controlPoints[i3D2] - facesCenter[j3D2];
            p.z = controlPoints[i3D3] - facesCenter[j3D3];

            pLocal.x = p.x * e1jPoint.x + p.y * e1jPoint.y + p.z * e1jPoint.z;
            pLocal.y = p.x * e2jPoint.x + p.y * e2jPoint.y + p.z * e2jPoint.z;
            pLocal.z = p.x * e3jPoint.x + p.y * e3jPoint.y + p.z * e3jPoint.z;
            
            p1Local.x = p1[j2D1]; p1Local.y = p1[j2D2];
            p2Local.x = p2[j2D1]; p2Local.y = p2[j2D2];
            p3Local.x = p3[j2D1]; p3Local.y = p3[j2D2];

            sourceFunc(pLocal, p1Local, p2Local, p3Local, e1jPoint, e2jPoint, e3jPoint, facesAreas[j], facesMaxDistance[j], sourceVel);
            doubletFunc(pLocal, p1Local, p2Local, p3Local, e1jPoint, e2jPoint, e3jPoint, facesAreas[j], facesMaxDistance[j], doubletVel);

            matrix[i * n + j] = doubletVel[0] * e3iPoint.x + doubletVel[1] * e3iPoint.y + doubletVel[2] * e3iPoint.z;
            array[i] = array[i] - sigma[j] * (sourceVel[0] * e3iPoint.x + sourceVel[1] * e3iPoint.y + sourceVel[2] * e3iPoint.z);

            matrixVelx[i * n + j] = doubletVel[0];
            matrixVely[i * n + j] = doubletVel[1];
            matrixVelz[i * n + j] = doubletVel[2];

            arrayVel[i * 3] = arrayVel[i * 3] + sigma[j] * sourceVel[0];
            arrayVel[i * 3 + 1] = arrayVel[i * 3 + 1] + sigma[j] * sourceVel[1];
            arrayVel[i * 3 + 2] = arrayVel[i * 3 + 2] + sigma[j] * sourceVel[2];

        }

        array[i] = array[i] - (freestream[0] * e3iPoint.x + freestream[1] * e3iPoint.y + freestream[2] * e3iPoint.z);

        arrayVel[i * 3] = arrayVel[i * 3] + freestream[0];
        arrayVel[i * 3 + 1] = arrayVel[i * 3 + 1] + freestream[1];
        arrayVel[i * 3 + 2] = arrayVel[i * 3 + 2] + freestream[2];

        // Wake
        p.x = controlPoints[i3D1];
        p.y = controlPoints[i3D2];
        p.z = controlPoints[i3D3];

        // Left wing
        for (k = 0; k < nSpanLeftWing; k ++) {

            for (l = 0; l < nWakeLeftWing - 1; l ++) {
                
                indexLine1 = leftWingGrid[k * nWakeLeftWing + l] * 3;
                p1Line.x = leftWingVertices[indexLine1];
                p1Line.y = leftWingVertices[indexLine1 + 1];
                p1Line.z = leftWingVertices[indexLine1 + 2];

                indexLine2 = leftWingGrid[k * nWakeLeftWing + l + 1] * 3;
                p2Line.x = leftWingVertices[indexLine2];
                p2Line.y = leftWingVertices[indexLine2 + 1];
                p2Line.z = leftWingVertices[indexLine2 + 2];

                lineFunc(p, p1Line, p2Line, lineVel);

                if (k == 0) {
                    matrix[i * n + leftWingFaces[k * 2]] = matrix[i * n + leftWingFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[k * 2 + 1]] = matrix[i * n + leftWingFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                } else if (k == nSpanLeftWing - 1) {
                    matrix[i * n + leftWingFaces[(k - 1) * 2]] = matrix[i * n + leftWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                } else {
                    matrix[i * n + leftWingFaces[(k - 1) * 2]] = matrix[i * n + leftWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrix[i * n + leftWingFaces[k * 2]] = matrix[i * n + leftWingFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[k * 2 + 1]] = matrix[i * n + leftWingFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    indexLine1 = leftWingGrid[(k - 1) * nWakeLeftWing] * 3;
                    p1Line.x = leftWingVertices[indexLine1];
                    p1Line.y = leftWingVertices[indexLine1 + 1];
                    p1Line.z = leftWingVertices[indexLine1 + 2];

                    indexLine2 = leftWingGrid[k * nWakeLeftWing] * 3;
                    p2Line.x = leftWingVertices[indexLine2];
                    p2Line.y = leftWingVertices[indexLine2 + 1];
                    p2Line.z = leftWingVertices[indexLine2 + 2];

                    lineFunc(p, p1Line, p2Line, lineVel);

                    matrix[i * n + leftWingFaces[(k - 1) * 2]] = matrix[i * n + leftWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                }

            }

        }

        // Right wing
        for (k = 0; k < nSpanRightWing; k ++) {

            for (l = 0; l < nWakeRightWing - 1; l ++) {
                
                indexLine1 = rightWingGrid[k * nWakeRightWing + l] * 3;
                p1Line.x = rightWingVertices[indexLine1];
                p1Line.y = rightWingVertices[indexLine1 + 1];
                p1Line.z = rightWingVertices[indexLine1 + 2];

                indexLine2 = rightWingGrid[k * nWakeRightWing + l + 1] * 3;
                p2Line.x = rightWingVertices[indexLine2];
                p2Line.y = rightWingVertices[indexLine2 + 1];
                p2Line.z = rightWingVertices[indexLine2 + 2];

                lineFunc(p, p1Line, p2Line, lineVel);

                if (k == 0) {
                    matrix[i * n + rightWingFaces[k * 2]] = matrix[i * n + rightWingFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[k * 2 + 1]] = matrix[i * n + rightWingFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                } else if (k == nSpanRightWing - 1) {
                    matrix[i * n + rightWingFaces[(k - 1) * 2]] = matrix[i * n + rightWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                } else {
                    matrix[i * n + rightWingFaces[(k - 1) * 2]] = matrix[i * n + rightWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrix[i * n + rightWingFaces[k * 2]] = matrix[i * n + rightWingFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[k * 2 + 1]] = matrix[i * n + rightWingFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    indexLine1 = rightWingGrid[(k - 1) * nWakeRightWing] * 3;
                    p1Line.x = rightWingVertices[indexLine1];
                    p1Line.y = rightWingVertices[indexLine1 + 1];
                    p1Line.z = rightWingVertices[indexLine1 + 2];

                    indexLine2 = rightWingGrid[k * nWakeRightWing] * 3;
                    p2Line.x = rightWingVertices[indexLine2];
                    p2Line.y = rightWingVertices[indexLine2 + 1];
                    p2Line.z = rightWingVertices[indexLine2 + 2];

                    lineFunc(p, p1Line, p2Line, lineVel);

                    matrix[i * n + rightWingFaces[(k - 1) * 2]] = matrix[i * n + rightWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                }

            }

        }

        // Tail
        for (k = 0; k < nSpanTail; k ++) {

            for (l = 0; l < nWakeTail - 1; l ++) {
                
                indexLine1 = tailGrid[k * nWakeTail + l] * 3;
                p1Line.x = tailVertices[indexLine1];
                p1Line.y = tailVertices[indexLine1 + 1];
                p1Line.z = tailVertices[indexLine1 + 2];

                indexLine2 = tailGrid[k * nWakeTail + l + 1] * 3;
                p2Line.x = tailVertices[indexLine2];
                p2Line.y = tailVertices[indexLine2 + 1];
                p2Line.z = tailVertices[indexLine2 + 2];

                lineFunc(p, p1Line, p2Line, lineVel);

                if (k == 0) {
                    matrix[i * n + tailFaces[k * 2]] = matrix[i * n + tailFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[k * 2 + 1]] = matrix[i * n + tailFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                } else if (k == nSpanTail - 1) {
                    matrix[i * n + tailFaces[(k - 1) * 2]] = matrix[i * n + tailFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[(k - 1) * 2 + 1]] = matrix[i * n + tailFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                } else {
                    matrix[i * n + tailFaces[(k - 1) * 2]] = matrix[i * n + tailFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[(k - 1) * 2 + 1]] = matrix[i * n + tailFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrix[i * n + tailFaces[k * 2]] = matrix[i * n + tailFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[k * 2 + 1]] = matrix[i * n + tailFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    indexLine1 = tailGrid[(k - 1) * nWakeTail] * 3;
                    p1Line.x = tailVertices[indexLine1];
                    p1Line.y = tailVertices[indexLine1 + 1];
                    p1Line.z = tailVertices[indexLine1 + 2];

                    indexLine2 = tailGrid[k * nWakeTail] * 3;
                    p2Line.x = tailVertices[indexLine2];
                    p2Line.y = tailVertices[indexLine2 + 1];
                    p2Line.z = tailVertices[indexLine2 + 2];

                    lineFunc(p, p1Line, p2Line, lineVel);

                    matrix[i * n + tailFaces[(k - 1) * 2]] = matrix[i * n + tailFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[(k - 1) * 2 + 1]] = matrix[i * n + tailFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                }

            }

        }

    }
}
