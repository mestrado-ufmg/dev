#include "lib.h"
#include <math.h>
#include <stdio.h>

/// Constants
const double ZERO_ERROR = 1e-8;
const double PI = 3.14159265359;
const double FACTOR = 1 / (4 * PI);

/// Helpers
double division(double a, double b) {
    if ((-ZERO_ERROR < b) && (b < 0)) {
        return - a / ZERO_ERROR;
    } else if ((0 < b) && (b < ZERO_ERROR)) {
        return a / ZERO_ERROR;
    } else {
        return a / b;
    }
}

double norm(struct Point3D p) {
    return sqrt(pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2));
}

struct Point3D cross(struct Point3D p1, struct Point3D p2) {
    struct Point3D p = {p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z, p1.x * p2.y - p1.y * p2.x};
    return p;
}

/// Singularities
void source(struct Point3D p, struct Point2D p1, struct Point2D p2, struct Point2D p3, struct Point3D e1, struct Point3D e2, struct Point3D e3, double area, double maxDistance, double* vel) {

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

void line(struct Point3D p, struct Point3D p1, struct Point3D p2, double* vel) {

    struct Point3D r1 = {p1.x - p.x, p1.y - p.y, p1.z - p.z};
    struct Point3D r2 = {p2.x - p.x, p2.y - p.y, p2.z - p.z};

    struct Point3D r1xr2 = cross(r1, r2);

    double r1Norm = norm(r1);
    double r2Norm = norm(r2);

    double r1xr2Norm2 = pow(norm(r1xr2), 2);

    double dot = (1 / r1xr2Norm2) * ((r1.x - r2.x) * (r1.x / r1Norm - r2.x / r2Norm) + (r1.y - r2.y) * (r1.y / r1Norm - r2.y / r2Norm) + (r1.z - r2.z) * (r1.z / r1Norm - r2.z / r2Norm));

    vel[0] = FACTOR * r1xr2.x * dot;
    vel[1] = FACTOR * r1xr2.y * dot;
    vel[2] = FACTOR * r1xr2.z * dot;
}

void doublet(struct Point3D p, struct Point3D p1, struct Point3D p2, struct Point3D p3, struct Point3D e1, struct Point3D e2, struct Point3D e3, double area, double maxDistance, double* vel) {

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

        double vel1[3];
        double vel2[3];
        double vel3[3];

        line(p, p1, p2, vel1);
        line(p, p2, p3, vel2);
        line(p, p3, p1, vel3);

        vel[0] = vel1[0] + vel2[0] + vel3[0];
        vel[1] = vel1[1] + vel2[1] + vel3[1];
        vel[2] = vel1[2] + vel2[2] + vel3[2];

    }

}

/// Create linear system
void create(struct input_struct *input, struct output_struct *output) {

    // For loop
    int i, j; // , k;

    // Point parameters
    double px_inertial, py_inertial, pz_inertial;
    double px, py, pz;

    // Face points parameters
    double p1_inertial[input->nf][3], p2_inertial[input->nf][3], p3_inertial[input->nf][3];
    double p1[input->nf][2], p2[input->nf][2], p3[input->nf][2];
    struct Point3D point3D1, point3D2, point3D3;
    struct Point2D point2D1, point2D2, point2D3;
    struct Point3D e13D, e23D, e33D;
    struct Point3D point3D;
    struct Point3D pointInertial3D;

    // Veloity parameters
    double sourceVel[3];
    double doubletVel[3];
    // double lineVel[3];
    // double lineVelAux;

    // Face parameters
    for (i = 0; i < input->nf; i++) {
        p1_inertial[i][0] = input->vertices[input->faces[i][0]][0] - input->facesCenter[i][0];
        p1_inertial[i][1] = input->vertices[input->faces[i][0]][1] - input->facesCenter[i][1];
        p1_inertial[i][2] = input->vertices[input->faces[i][0]][2] - input->facesCenter[i][2];

        p2_inertial[i][0] = input->vertices[input->faces[i][1]][0] - input->facesCenter[i][0];
        p2_inertial[i][1] = input->vertices[input->faces[i][1]][1] - input->facesCenter[i][1];
        p2_inertial[i][2] = input->vertices[input->faces[i][1]][2] - input->facesCenter[i][2];

        p3_inertial[i][0] = input->vertices[input->faces[i][2]][0] - input->facesCenter[i][0];
        p3_inertial[i][1] = input->vertices[input->faces[i][2]][1] - input->facesCenter[i][1];
        p3_inertial[i][2] = input->vertices[input->faces[i][2]][2] - input->facesCenter[i][2];

        p1[i][0] = p1_inertial[i][0] * input->e1[i][0] + p1_inertial[i][1] * input->e1[i][1] + p1_inertial[i][2] * input->e1[i][2];
        p1[i][1] = p1_inertial[i][0] * input->e2[i][0] + p1_inertial[i][1] * input->e2[i][1] + p1_inertial[i][2] * input->e2[i][2];

        p2[i][0] = p2_inertial[i][0] * input->e1[i][0] + p2_inertial[i][1] * input->e1[i][1] + p2_inertial[i][2] * input->e1[i][2];
        p2[i][1] = p2_inertial[i][0] * input->e2[i][0] + p2_inertial[i][1] * input->e2[i][1] + p2_inertial[i][2] * input->e2[i][2];

        p3[i][0] = p3_inertial[i][0] * input->e1[i][0] + p3_inertial[i][1] * input->e1[i][1] + p3_inertial[i][2] * input->e1[i][2];
        p3[i][1] = p3_inertial[i][0] * input->e2[i][0] + p3_inertial[i][1] * input->e2[i][1] + p3_inertial[i][2] * input->e2[i][2];
    }

    // Surface
    for (i = 0; i < input->nf; i++) {

        // Face j acting on i
        for (j = 0; j < input->nf; j++) {

            // Point
            px_inertial = input->controlPoints[i][0] - input->facesCenter[j][0];
            py_inertial = input->controlPoints[i][1] - input->facesCenter[j][1];
            pz_inertial = input->controlPoints[i][2] - input->facesCenter[j][2];

            px = px_inertial * input->e1[j][0] + py_inertial * input->e1[j][1] + pz_inertial * input->e1[j][2];
            py = px_inertial * input->e2[j][0] + py_inertial * input->e2[j][1] + pz_inertial * input->e2[j][2];
            pz = px_inertial * input->e3[j][0] + py_inertial * input->e3[j][1] + pz_inertial * input->e3[j][2];

            // Velocity coefficients
            point3D.x = px; point3D.y = py; point3D.z = pz;
            pointInertial3D.x = px_inertial; pointInertial3D.y = py_inertial; pointInertial3D.z = pz_inertial;

            point3D1.x = p1_inertial[j][0]; point3D1.y = p1_inertial[j][1]; point3D1.z = p1_inertial[j][2];
            point3D2.x = p2_inertial[j][0]; point3D2.y = p2_inertial[j][1]; point3D2.z = p2_inertial[j][2];
            point3D3.x = p3_inertial[j][0]; point3D3.y = p3_inertial[j][1]; point3D3.z = p3_inertial[j][2];

            point2D1.x = p1[j][0]; point2D1.y = p1[j][1];
            point2D2.x = p2[j][0]; point2D2.y = p2[j][1];
            point2D3.x = p3[j][0]; point2D3.y = p3[j][1];

            e13D.x = input->e1[j][0]; e13D.y = input->e1[j][1]; e13D.z = input->e1[j][2];
            e23D.x = input->e2[j][0]; e23D.y = input->e2[j][1]; e23D.z = input->e2[j][2];
            e33D.x = input->e3[j][0]; e33D.y = input->e3[j][1]; e33D.z = input->e3[j][2];

            source(point3D, point2D1, point2D2, point2D3, e13D, e23D, e33D, input->facesAreas[j], input->facesMaxDistance[j], sourceVel);
            doublet(pointInertial3D, point3D1, point3D2, point3D3, e13D, e23D, e33D, input->facesAreas[j], input->facesMaxDistance[j], doubletVel);

            output->velXMatrix[i][j] = doubletVel[0];
            output->velYMatrix[i][j] = doubletVel[1];
            output->velZMatrix[i][j] = doubletVel[2];

            output->velArray[i][0] = output->velArray[i][0] + input->sigma[j] * sourceVel[0];
            output->velArray[i][1] = output->velArray[i][1] + input->sigma[j] * sourceVel[1];
            output->velArray[i][2] = output->velArray[i][2] + input->sigma[j] * sourceVel[2];

            output->matrix[i][j] = doubletVel[0] * input->e3[i][0] + doubletVel[1] * input->e3[i][1] + doubletVel[2] * input->e3[i][2];
            
            output->array[i] = output->array[i] - (input->sigma[j] * sourceVel[0] * input->e3[i][0] + input->sigma[j] * sourceVel[1] * input->e3[i][1] + input->sigma[j] * sourceVel[2] * input->e3[i][2]);
        }

        // // Left wing wake
        // px_inertial = input->controlPoints[i][0];
        // py_inertial = input->controlPoints[i][1];
        // pz_inertial = input->controlPoints[i][2];

        // pointInertial3D.x = px_inertial; pointInertial3D.y = py_inertial; pointInertial3D.z = pz_inertial;

        // // Loop over span
        // for (j = 0; j < input->nSpanLeftWing; j++) {

        //     // Loop over wake points
        //     for (k = 0; k < input->nWakeLeftWing - 1; j++) {

        //         point3D1.x = input->leftWingWakeVertices[input->leftWingWakeGrid[j][k]][0];
        //         point3D1.y = input->leftWingWakeVertices[input->leftWingWakeGrid[j][k]][1];
        //         point3D1.z = input->leftWingWakeVertices[input->leftWingWakeGrid[j][k]][2];

        //         point3D2.x = input->leftWingWakeVertices[input->leftWingWakeGrid[j][k + 1]][0];
        //         point3D2.y = input->leftWingWakeVertices[input->leftWingWakeGrid[j][k + 1]][1];
        //         point3D2.z = input->leftWingWakeVertices[input->leftWingWakeGrid[j][k + 1]][2];

        //         line(pointInertial3D, point3D1, point3D2, lineVel);

        //         // The filament always goes in the direction of the wake. Then is necessary to multiply by -1.
        //         lineVelAux = lineVel[0] * input->e3[i][0] + lineVel[1] * input->e3[i][1] + lineVel[2] * input->e3[i][2];
        //         output->matrix[i][input->nf + j] = -lineVelAux;

        //         if ((j != 0) && (j != input->nSpanLeftWing - 1)) {
        //             output->matrix[i][input->nf + j - 1] = lineVelAux;
        //         }
        //     }

        //     // Filament on the trailing edge
        //     if (j != input->nSpanLeftWing - 1) {
        //         point3D1.x = input->leftWingWakeVertices[input->leftWingWakeGrid[j][0]][0];
        //         point3D1.y = input->leftWingWakeVertices[input->leftWingWakeGrid[j][0]][1];
        //         point3D1.z = input->leftWingWakeVertices[input->leftWingWakeGrid[j][0]][2];

        //         point3D2.x = input->leftWingWakeVertices[input->leftWingWakeGrid[j + 1][0]][0];
        //         point3D2.y = input->leftWingWakeVertices[input->leftWingWakeGrid[j + 1][0]][1];
        //         point3D2.z = input->leftWingWakeVertices[input->leftWingWakeGrid[j + 1][0]][2];

        //         line(pointInertial3D, point3D1, point3D2, lineVel);

        //         lineVelAux = lineVel[0] * input->e3[i][0] + lineVel[1] * input->e3[i][1] + lineVel[2] * input->e3[i][2];
        //         output->matrix[i][input->nf + j] = -lineVelAux;
        //     }

        // }

        // Freestream
        output->array[i] = output->array[i] - (input->freestream[0] * input->e3[i][0] + input->freestream[1] * input->e3[i][1] + input->freestream[2] * input->e3[i][2]);

        output->velArray[i][0] = output->velArray[i][0] + input->freestream[0];
        output->velArray[i][1] = output->velArray[i][1] + input->freestream[1];
        output->velArray[i][2] = output->velArray[i][2] + input->freestream[2];

    }

}