#include "../helpers/structs.h"
#include "../helpers/customMath.h"
#include "data.h"
#include <math.h>

void sourceFunc(struct Point p, struct Point p1, struct Point p2, struct Point p3, struct Point e1, struct Point e2, struct Point e3, double area, double maxDistance, double *vel)
/* Indulced velocity by a source distribution over a triangular panel */
{

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

        u = -FACTOR * (division((p2.y - p1.y), d12) * ln12 + division((p3.y - p2.y), d23) * ln23 + division((p1.y - p3.y), d31) * ln31);
        v = FACTOR * (division((p2.x - p1.x), d12) * ln12 + division((p3.x - p2.x), d23) * ln23 + division((p1.x - p3.x), d31) * ln31);
        w = -FACTOR * (atan(division(m12 * l1 - h1, p.z * r1)) - atan(division(m12 * l2 - h2, p.z * r2)) + atan(division(m23 * l2 - h2, p.z * r2)) - atan(division(m23 * l3 - h3, p.z * r3)) + atan(division(m31 * l3 - h3, p.z * r3)) - atan(division(m31 * l1 - h1, p.z * r1)));
    
    }

    vel[0] = u * e1.x + v * e2.x + w * e3.x;
    vel[1] = u * e1.y + v * e2.y + w * e3.y;
    vel[2] = u * e1.z + v * e2.z + w * e3.z;

}

void lineFunc(struct Point p, struct Point p1, struct Point p2, double *vel)
/* Indulced velocity by a source distribution over a line */
{

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

void doubletFunc(struct Point p, struct Point p1, struct Point p2, struct Point p3, struct Point e1, struct Point e2, struct Point e3, double area, double maxDistance, double *vel)
/* Indulced velocity by a doublet distribution over a triangular panel */
{

    double distance = norm(p);

    if (distance > maxDistance) {

        double pxLocal = p.x * e1.x + p.y * e1.y + p.z * e1.z;
        double pyLocal = p.x * e2.x + p.y * e2.y + p.z * e2.z;
        double pzLocal = p.x * e3.x + p.y * e3.y + p.z * e3.z;
        double den = pow(pxLocal * pxLocal + pyLocal * pyLocal + pzLocal * pzLocal, 2.5);

        double u = 0.75 * FACTOR * area * pzLocal * pxLocal / den;
        double v = 0.75 * FACTOR * area * pzLocal * pyLocal / den;
        double w = -FACTOR * area * (pxLocal * pxLocal + pyLocal * pyLocal - 2 * pzLocal * pzLocal) / den;

        vel[0] = u * e1.x + v * e2.x + w * e3.x;
        vel[1] = u * e1.y + v * e2.y + w * e3.y;
        vel[2] = u * e1.z + v * e2.z + w * e3.z;

    } else {

        double u, v, w;

        double *vel1 = (double *)malloc(3 * sizeof(double *));
        double *vel2 = (double *)malloc(3 * sizeof(double *));
        double *vel3 = (double *)malloc(3 * sizeof(double *));

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

void addWakeCoefficients(struct Input input, double *lineVel, struct Point e3iPoint, int face, int nWake, int nSpanWake, double *wakeVertices, int *wakeGrid, int *wakeFaces, struct PotentialFlowData data)
{

    int k, l;
    int indexLine1, indexLine2;
    struct Point p1Line;
    struct Point p2Line;
    struct Point p;

    // Wake
    p.x = input.mesh.surface.controlPoints[3 * face];
    p.y = input.mesh.surface.controlPoints[3 * face + 1];
    p.z = input.mesh.surface.controlPoints[3 * face + 2];

    // Left wing
    for (k = 0; k < nSpanWake; k++) {

        for (l = 0; l < nWake - 1; l++) {

            indexLine1 = wakeGrid[k * nWake + l] * 3;
            p1Line.x = wakeVertices[indexLine1];
            p1Line.y = wakeVertices[indexLine1 + 1];
            p1Line.z = wakeVertices[indexLine1 + 2];

            indexLine2 = wakeGrid[k * nWake + l + 1] * 3;
            p2Line.x = wakeVertices[indexLine2];
            p2Line.y = wakeVertices[indexLine2 + 1];
            p2Line.z = wakeVertices[indexLine2 + 2];

            lineFunc(p, p1Line, p2Line, lineVel);

            if (k == 0) {

                data.a[face * input.mesh.surface.nf + wakeFaces[k * 2]] = data.a[face * input.mesh.surface.nf + wakeFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                data.a[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] = data.a[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[k * 2]] = data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[k * 2]] - lineVel[0];
                data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[k * 2]] = data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[k * 2]] - lineVel[1];
                data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[k * 2]] = data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[k * 2]] - lineVel[2];

                data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] = data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] + lineVel[0];
                data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] = data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] + lineVel[1];
                data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] = data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] + lineVel[2];
            } else if (k == nSpanWake - 1) {

                data.a[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] = data.a[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                data.a[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] = data.a[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] = data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] + lineVel[0];
                data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] = data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] + lineVel[1];
                data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] = data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] + lineVel[2];

                data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] = data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] - lineVel[0];
                data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] = data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] - lineVel[1];
                data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] = data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] - lineVel[2];
            } else {

                data.a[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] = data.a[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                data.a[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] = data.a[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] = data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] + lineVel[0];
                data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] = data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] + lineVel[1];
                data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] = data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] + lineVel[2];

                data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] = data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] - lineVel[0];
                data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] = data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] - lineVel[1];
                data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] = data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] - lineVel[2];

                data.a[face * input.mesh.surface.nf + wakeFaces[k * 2]] = data.a[face * input.mesh.surface.nf + wakeFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                data.a[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] = data.a[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[k * 2]] = data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[k * 2]] - lineVel[0];
                data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[k * 2]] = data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[k * 2]] - lineVel[1];
                data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[k * 2]] = data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[k * 2]] - lineVel[2];

                data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] = data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] + lineVel[0];
                data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] = data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] + lineVel[1];
                data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] = data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[k * 2 + 1]] + lineVel[2];

                indexLine1 = wakeGrid[(k - 1) * nWake] * 3;
                p1Line.x = wakeVertices[indexLine1];
                p1Line.y = wakeVertices[indexLine1 + 1];
                p1Line.z = wakeVertices[indexLine1 + 2];

                indexLine2 = wakeGrid[k * nWake] * 3;
                p2Line.x = wakeVertices[indexLine2];
                p2Line.y = wakeVertices[indexLine2 + 1];
                p2Line.z = wakeVertices[indexLine2 + 2];

                lineFunc(p, p1Line, p2Line, lineVel);

                data.a[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] = data.a[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                data.a[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] = data.a[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] = data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] + lineVel[0];
                data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] = data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] + lineVel[1];
                data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] = data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2]] + lineVel[2];

                data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] = data.a_vel_x[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] - lineVel[0];
                data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] = data.a_vel_y[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] - lineVel[1];
                data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] = data.a_vel_z[face * input.mesh.surface.nf + wakeFaces[(k - 1) * 2 + 1]] - lineVel[2];
            }
        }
    }

}

void getLinearSystemImp(struct Input input, struct PotentialFlowData data)
{
    /* Parameters */

    // Loops
    int i, j;

    // Point
    struct Point p;
    struct Point p1Line;
    struct Point p2Line;
    struct Point pLocal;
    struct Point p1Local;
    p1Local.z = 0.0;
    struct Point p2Local;
    p2Local.z = 0.0;
    struct Point p3Local;
    p3Local.z = 0.0;

    // Base vectors
    struct Point e3iPoint;

    struct Point e1jPoint;
    struct Point e2jPoint;
    struct Point e3jPoint;

    // Velocities
    double *sourceVel = (double *)malloc(3 * sizeof(double));
    double *doubletVel = (double *)malloc(3 * sizeof(double));
    double *lineVel = (double *)malloc(3 * sizeof(double));

    /* Create */
    for (i = 0; i < input.mesh.surface.nf; i++)
    {

        e3iPoint.x = input.mesh.surface.e3[i * 3];
        e3iPoint.y = input.mesh.surface.e3[i * 3 + 1];
        e3iPoint.z = input.mesh.surface.e3[i * 3 + 2];

        data.rhs[i] = 0.0;
        data.rhs_vel_x[i] = 0.0;
        data.rhs_vel_y[i] = 0.0;
        data.rhs_vel_z[i] = 0.0;

        /* Surface */
        for (j = 0; j < input.mesh.surface.nf; j++) // Effect of j on i
        {

            // Points
            e1jPoint.x = input.mesh.surface.e1[j * 3];
            e1jPoint.y = input.mesh.surface.e1[j * 3 + 1];
            e1jPoint.z = input.mesh.surface.e1[j * 3 + 2];

            e2jPoint.x = input.mesh.surface.e2[j * 3];
            e2jPoint.y = input.mesh.surface.e2[j * 3 + 1];
            e2jPoint.z = input.mesh.surface.e2[j * 3 + 2];

            e3jPoint.x = input.mesh.surface.e3[j * 3];
            e3jPoint.y = input.mesh.surface.e3[j * 3 + 1];
            e3jPoint.z = input.mesh.surface.e3[j * 3 + 2];

            p.x = input.mesh.surface.controlPoints[i * 3] - input.mesh.surface.facesCenter[j * 3];
            p.y = input.mesh.surface.controlPoints[i * 3 + 1] - input.mesh.surface.facesCenter[j * 3 + 1];
            p.z = input.mesh.surface.controlPoints[i * 3 + 2] - input.mesh.surface.facesCenter[j * 3 + 2];

            pLocal.x = p.x * e1jPoint.x + p.y * e1jPoint.y + p.z * e1jPoint.z;
            pLocal.y = p.x * e2jPoint.x + p.y * e2jPoint.y + p.z * e2jPoint.z;
            pLocal.z = p.x * e3jPoint.x + p.y * e3jPoint.y + p.z * e3jPoint.z;

            p1Local.x = input.mesh.surface.p1[j * 2];
            p1Local.y = input.mesh.surface.p1[j * 2 + 1];

            p2Local.x = input.mesh.surface.p2[j * 2];
            p2Local.y = input.mesh.surface.p2[j * 2 + 1];

            p3Local.x = input.mesh.surface.p3[j * 2];
            p3Local.y = input.mesh.surface.p3[j * 2 + 1];

            sourceFunc(pLocal, p1Local, p2Local, p3Local, e1jPoint, e2jPoint, e3jPoint, input.mesh.surface.facesAreas[j], input.mesh.surface.facesMaxDistance[j], sourceVel);
            doubletFunc(pLocal, p1Local, p2Local, p3Local, e1jPoint, e2jPoint, e3jPoint, input.mesh.surface.facesAreas[j], input.mesh.surface.facesMaxDistance[j], doubletVel);

            data.a[i * input.mesh.surface.nf + j] = doubletVel[0] * e3iPoint.x + doubletVel[1] * e3iPoint.y + doubletVel[2] * e3iPoint.z;
            data.rhs[i] = data.rhs[i] - data.sigma[j] * (sourceVel[0] * e3iPoint.x + sourceVel[1] * e3iPoint.y + sourceVel[2] * e3iPoint.z);

            data.a_vel_x[i * input.mesh.surface.nf + j] = doubletVel[0];
            data.a_vel_y[i * input.mesh.surface.nf + j] = doubletVel[1];
            data.a_vel_z[i * input.mesh.surface.nf + j] = doubletVel[2];
            
            data.rhs_vel_x[i] = data.rhs_vel_x[i] + data.sigma[j] * sourceVel[0];
            data.rhs_vel_y[i] = data.rhs_vel_y[i] + data.sigma[j] * sourceVel[1];
            data.rhs_vel_z[i] = data.rhs_vel_z[i] + data.sigma[j] * sourceVel[2];

            data.ia[i * input.mesh.surface.nf + j] = i;
            data.ja[i * input.mesh.surface.nf + j] = j;
        }

        data.rhs[i] = data.rhs[i] - (input.environment.vel_x * e3iPoint.x + input.environment.vel_y * e3iPoint.y + input.environment.vel_z * e3iPoint.z);

        data.rhs_vel_x[i] = data.rhs_vel_x[i] + input.environment.vel_x;
        data.rhs_vel_y[i] = data.rhs_vel_y[i] + input.environment.vel_y;
        data.rhs_vel_z[i] = data.rhs_vel_z[i] + input.environment.vel_z;

        /* Wake */
        // addWakeCoefficients(input, lineVel, e3iPoint, i, input.mesh.wake.tail.nWake, input.mesh.wake.tail.nSpan, input.mesh.wake.tail.vertices, input.mesh.wake.tail.grid, input.mesh.wake.tail.faces, data);

    }

    /* Free */
    free(sourceVel);
    free(doubletVel);
    free(lineVel);
}