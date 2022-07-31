#include "lib.h"
#include <math.h>
#include <stdio.h>

/// Constants
const double ZERO_ERROR = 1e-8;
const double PI = 3.14159265359;

/// Check division by 0
double division(double a, double b) {
    if ((-ZERO_ERROR < b) && (b < 0)) {
        return - a / ZERO_ERROR;
    } else if ((0 < b) && (b < ZERO_ERROR)) {
        return a / ZERO_ERROR;
    } else {
        return a / b;
    }
}

/// Create linear system
void create(struct input_struct *input, struct output_struct *output) {

    // For loop
    int i, j;

    // Point parameters
    double px_inertial, py_inertial, pz_inertial;
    double px, py, pz;

    // Face points parameters
    double p1x_inertial, p1y_inertial, p1z_inertial;
    double p2x_inertial, p2y_inertial, p2z_inertial;
    double p3x_inertial, p3y_inertial, p3z_inertial;
    double p1x, p1y;
    double p2x, p2y;
    double p3x, p3y;

    // Veloity parameters
    double r1, r2, r3;
    double l1, l2, l3;
    double h1, h2, h3;

    double d12, d23, d31;
    double m12, m23, m31;
    double ln12, ln23, ln31;
    double u, v, w;
    const double factor = 1 / (4 * PI);

    // Equations
    for (i = 0; i < input->n; i++) {

        // Face j acting on i
        for (j = 0; j < input->n; j++) {

            // Faces sides
            p1x_inertial = input->vertices[input->faces[j][0]][0] - input->facesCenter[j][0];
            p1y_inertial = input->vertices[input->faces[j][0]][1] - input->facesCenter[j][1];
            p1z_inertial = input->vertices[input->faces[j][0]][2] - input->facesCenter[j][2];

            p2x_inertial = input->vertices[input->faces[j][1]][0] - input->facesCenter[j][0];
            p2y_inertial = input->vertices[input->faces[j][1]][1] - input->facesCenter[j][1];
            p2z_inertial = input->vertices[input->faces[j][1]][2] - input->facesCenter[j][2];

            p3x_inertial = input->vertices[input->faces[j][2]][0] - input->facesCenter[j][0];
            p3y_inertial = input->vertices[input->faces[j][2]][1] - input->facesCenter[j][1];
            p3z_inertial = input->vertices[input->faces[j][2]][2] - input->facesCenter[j][2];

            p1x = p1x_inertial * input->e1[j][0] + p1y_inertial * input->e1[j][1] + p1z_inertial * input->e1[j][2];
            p1y = p1x_inertial * input->e2[j][0] + p1y_inertial * input->e2[j][1] + p1z_inertial * input->e2[j][2];
            // p1z = p1x_inertial * input->e3[j][0] + p1y_inertial * input->e3[j][1] + p1z_inertial * input->e3[j][2];

            p2x = p2x_inertial * input->e1[j][0] + p2y_inertial * input->e1[j][1] + p2z_inertial * input->e1[j][2];
            p2y = p2x_inertial * input->e2[j][0] + p2y_inertial * input->e2[j][1] + p2z_inertial * input->e2[j][2];
            // p2z = p2x_inertial * input->e3[j][0] + p2y_inertial * input->e3[j][1] + p2z_inertial * input->e3[j][2];

            p3x = p3x_inertial * input->e1[j][0] + p3y_inertial * input->e1[j][1] + p3z_inertial * input->e1[j][2];
            p3y = p3x_inertial * input->e2[j][0] + p3y_inertial * input->e2[j][1] + p3z_inertial * input->e2[j][2];
            // p3z = p3x_inertial * input->e3[j][0] + p3y_inertial * input->e3[j][1] + p3z_inertial * input->e3[j][2];

            // Velocity induced by face j on face i
            px_inertial = input->controlPoints[i][0] - input->facesCenter[j][0];
            py_inertial = input->controlPoints[i][1] - input->facesCenter[j][1];
            pz_inertial = input->controlPoints[i][2] - input->facesCenter[j][2];

            px = px_inertial * input->e1[j][0] + py_inertial * input->e1[j][1] + pz_inertial * input->e1[j][2];
            py = px_inertial * input->e2[j][0] + py_inertial * input->e2[j][1] + pz_inertial * input->e2[j][2];
            pz = px_inertial * input->e3[j][0] + py_inertial * input->e3[j][1] + pz_inertial * input->e3[j][2];

            // Velocity coefficients
            r1 = sqrt(pow(px - p1x, 2) + pow(py - p1y, 2) + pow(pz, 2));
            r2 = sqrt(pow(px - p2x, 2) + pow(py - p2y, 2) + pow(pz, 2));
            r3 = sqrt(pow(px - p3x, 2) + pow(py - p3y, 2) + pow(pz, 2));

            l1 = pow(px - p1x, 2) + pow(pz, 2);
            l2 = pow(px - p2x, 2) + pow(pz, 2);
            l3 = pow(px - p3x, 2) + pow(pz, 2);

            h1 = (px - p1x) * (py - p1y);
            h2 = (px - p2x) * (py - p2y);
            h3 = (px - p3x) * (py - p3y);

            d12 = sqrt(pow(p2x - p1x, 2) + pow(p2y - p1y, 2));
            m12 = division(p2y - p1y, p2x - p1x);

            d23 = sqrt(pow(p3x - p2x, 2) + pow(p3y - p2y, 2));
            m23 = division(p3y - p2y, p3x - p2x);

            d31 = sqrt(pow(p1x - p3x, 2) + pow(p1y - p3y, 2));
            m31 = division(p1y - p3y, p1x - p3x);

            ln12 = log(division(r1 + r2 - d12, r1 + r2 + d12));
            ln23 = log(division(r2 + r3 - d23, r2 + r3 + d23));
            ln31 = log(division(r3 + r1 - d31, r3 + r1 + d31));

            u = -factor * ( division((p2y - p1y), d12) * ln12 + division((p3y - p2y), d23) * ln23 + division((p1y - p3y), d31) * ln31 );
            v = factor * ( division((p2x - p1x), d12) * ln12 + division((p3x - p2x), d23) * ln23 + division((p1x - p3x), d31) * ln31 );
            w = -factor * ( atan(division(m12 * l1 - h1, pz * r1)) - atan(division(m12 * l2 - h2, pz * r2)) + atan(division(m23 * l2 - h2, pz * r2)) - atan(division(m23 * l3 - h3, pz * r3)) + atan(division(m31 * l3 - h3, pz * r3)) - atan(division(m31 * l1 - h1, pz * r1)) );

            output->matrix[i][j] = (u * input->e1[j][0] + v * input->e2[j][0] + w * input->e3[j][0]) * input->e3[i][0] + (u * input->e1[j][1] + v * input->e2[j][1] + w * input->e3[j][1]) * input->e3[i][1] + (u * input->e1[j][2] + v * input->e2[j][2] + w * input->e3[j][2]) * input->e3[i][2];
            
        }

        // Freestream
        output->array[i] = -(input->freestream[0] * input->e3[i][0] + input->freestream[1] * input->e3[i][1] + input->freestream[2] * input->e3[i][2]);

    }

}

/// Create linear system
void parameters(struct input_struct *input, struct output_vel_struct *output, double *sigma) {

    // For loop
    int i, j;

    // Point parameters
    double px_inertial, py_inertial, pz_inertial;
    double px, py, pz;

    // Face points parameters
    double p1x_inertial, p1y_inertial, p1z_inertial;
    double p2x_inertial, p2y_inertial, p2z_inertial;
    double p3x_inertial, p3y_inertial, p3z_inertial;
    double p1x, p1y;
    double p2x, p2y;
    double p3x, p3y;

    // Veloity parameters
    double r1, r2, r3;
    double l1, l2, l3;
    double h1, h2, h3;

    double d12, d23, d31;
    double m12, m23, m31;
    double ln12, ln23, ln31;
    double u, v, w;
    const double factor = 1 / (4 * PI);

    // Freestream
    double uNorm = sqrt(pow(input->freestream[0], 2) + pow(input->freestream[1], 2) + pow(input->freestream[2], 2));

    // Equations
    for (i = 0; i < input->n; i++) {

        // Face j acting on i
        for (j = 0; j < input->n; j++) {

            // Faces sides
            p1x_inertial = input->vertices[input->faces[j][0]][0] - input->facesCenter[j][0];
            p1y_inertial = input->vertices[input->faces[j][0]][1] - input->facesCenter[j][1];
            p1z_inertial = input->vertices[input->faces[j][0]][2] - input->facesCenter[j][2];

            p2x_inertial = input->vertices[input->faces[j][1]][0] - input->facesCenter[j][0];
            p2y_inertial = input->vertices[input->faces[j][1]][1] - input->facesCenter[j][1];
            p2z_inertial = input->vertices[input->faces[j][1]][2] - input->facesCenter[j][2];

            p3x_inertial = input->vertices[input->faces[j][2]][0] - input->facesCenter[j][0];
            p3y_inertial = input->vertices[input->faces[j][2]][1] - input->facesCenter[j][1];
            p3z_inertial = input->vertices[input->faces[j][2]][2] - input->facesCenter[j][2];

            p1x = p1x_inertial * input->e1[j][0] + p1y_inertial * input->e1[j][1] + p1z_inertial * input->e1[j][2];
            p1y = p1x_inertial * input->e2[j][0] + p1y_inertial * input->e2[j][1] + p1z_inertial * input->e2[j][2];
            // p1z = p1x_inertial * input->e3[j][0] + p1y_inertial * input->e3[j][1] + p1z_inertial * input->e3[j][2];

            p2x = p2x_inertial * input->e1[j][0] + p2y_inertial * input->e1[j][1] + p2z_inertial * input->e1[j][2];
            p2y = p2x_inertial * input->e2[j][0] + p2y_inertial * input->e2[j][1] + p2z_inertial * input->e2[j][2];
            // p2z = p2x_inertial * input->e3[j][0] + p2y_inertial * input->e3[j][1] + p2z_inertial * input->e3[j][2];

            p3x = p3x_inertial * input->e1[j][0] + p3y_inertial * input->e1[j][1] + p3z_inertial * input->e1[j][2];
            p3y = p3x_inertial * input->e2[j][0] + p3y_inertial * input->e2[j][1] + p3z_inertial * input->e2[j][2];
            // p3z = p3x_inertial * input->e3[j][0] + p3y_inertial * input->e3[j][1] + p3z_inertial * input->e3[j][2];

            // Velocity induced by face j on face i
            px_inertial = input->controlPoints[i][0] - input->facesCenter[j][0];
            py_inertial = input->controlPoints[i][1] - input->facesCenter[j][1];
            pz_inertial = input->controlPoints[i][2] - input->facesCenter[j][2];

            px = px_inertial * input->e1[j][0] + py_inertial * input->e1[j][1] + pz_inertial * input->e1[j][2];
            py = px_inertial * input->e2[j][0] + py_inertial * input->e2[j][1] + pz_inertial * input->e2[j][2];
            pz = px_inertial * input->e3[j][0] + py_inertial * input->e3[j][1] + pz_inertial * input->e3[j][2];

            // Velocity coefficients
            r1 = sqrt(pow(px - p1x, 2) + pow(py - p1y, 2) + pow(pz, 2));
            r2 = sqrt(pow(px - p2x, 2) + pow(py - p2y, 2) + pow(pz, 2));
            r3 = sqrt(pow(px - p3x, 2) + pow(py - p3y, 2) + pow(pz, 2));

            l1 = pow(px - p1x, 2) + pow(pz, 2);
            l2 = pow(px - p2x, 2) + pow(pz, 2);
            l3 = pow(px - p3x, 2) + pow(pz, 2);

            h1 = (px - p1x) * (py - p1y);
            h2 = (px - p2x) * (py - p2y);
            h3 = (px - p3x) * (py - p3y);

            d12 = sqrt(pow(p2x - p1x, 2) + pow(p2y - p1y, 2));
            m12 = division(p2y - p1y, p2x - p1x);

            d23 = sqrt(pow(p3x - p2x, 2) + pow(p3y - p2y, 2));
            m23 = division(p3y - p2y, p3x - p2x);

            d31 = sqrt(pow(p1x - p3x, 2) + pow(p1y - p3y, 2));
            m31 = division(p1y - p3y, p1x - p3x);

            ln12 = log(division(r1 + r2 - d12, r1 + r2 + d12));
            ln23 = log(division(r2 + r3 - d23, r2 + r3 + d23));
            ln31 = log(division(r3 + r1 - d31, r3 + r1 + d31));

            u = -sigma[j] * factor * ( division((p2y - p1y), d12) * ln12 + division((p3y - p2y), d23) * ln23 + division((p1y - p3y), d31) * ln31 );
            v = sigma[j] * factor * ( division((p2x - p1x), d12) * ln12 + division((p3x - p2x), d23) * ln23 + division((p1x - p3x), d31) * ln31 );
            w = -sigma[j] * factor * ( atan(division(m12 * l1 - h1, pz * r1)) - atan(division(m12 * l2 - h2, pz * r2)) + atan(division(m23 * l2 - h2, pz * r2)) - atan(division(m23 * l3 - h3, pz * r3)) + atan(division(m31 * l3 - h3, pz * r3)) - atan(division(m31 * l1 - h1, pz * r1)) );

            output->vel[i][0] = output->vel[i][0] + u;
            output->vel[i][1] = output->vel[i][1] + v;
            output->vel[i][2] = output->vel[i][2] + w;

            printf("%f, %f, %f\n", u, v, w);
            
        }

        output->vel[i][0] = output->vel[i][0] + input->freestream[0];
        output->vel[i][1] = output->vel[i][1] + input->freestream[1];
        output->vel[i][2] = output->vel[i][2] + input->freestream[2];

        output->velNorm[i] = sqrt(pow(output->vel[i][0], 2) + pow(output->vel[i][1], 2) + pow(output->vel[i][2], 2));
        
        output->pressure[i] = 0.5 * 1.225 * (uNorm - output->velNorm[i]);

    }

}