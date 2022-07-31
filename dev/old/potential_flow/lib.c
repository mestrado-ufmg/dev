#include "lib.h"
#include <math.h>
#include <stdio.h>

/// Constants
const double ZERO_ERROR = 1e-8;
const double PI = 3.14159265359;
const double CONVERGENCY = 1e-8;

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

/// Returns the absolute value
double abs_value(double a) {
    if (a < 0) {
        return -a;
    } else {
        return a;
    }
}

/// Solves a linear system
void linear_system_solver(int n, double matrix[n][n], double array[n], double *solution) {

    // For loop
    int i, j, k;

    // Relative error
    double error[n];
    int check = 0;

    // Max. interactions
    int max_int = 500;

    // Equations sum
    double sum;
    double new_val;

    // Loop
    for (k = 0; k < max_int; k++) {

        // Solve equations
        for (i = 0; i < n; i++) {

            sum = 0;

            for (j = 0; j < n; j++) {
                if (i != j) {
                    sum = sum + matrix[i][j] * solution[j];
                }
            }

            new_val = division(array[i] - sum, matrix[i][i]);
            error[i] = abs_value(division(new_val - solution[i], new_val));
            solution[i] = solution[i] * 0.8 + new_val * 0.2;

        }

        // Check convergency
        check = 1;
        for (j = 0; j < n; j++) {
            if (error[j] > CONVERGENCY) {
                check = 0;
                break;
            }
        }

        if (check == 1) {
            printf("Linear system converged\n");
            break;
        }

    }

}

void ls_wrapper(int n, double **a, double *b, double *x) {
    
    double matrix[n][n];
    double array[n];

    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            matrix[i][j] = a[i][j];
        }
        array[i] = b[i];
    }

    linear_system_solver(n, matrix, array, x);
}

/// Calcute the source distribution in order to generate a zero flow across each face.
void source(struct input_struct *input, struct output_struct *output, double velFieldMatrix[input->n][input->n][3]) {

    // For loop
    int i, j;

    // Linear system
    double matrix[input->n][input->n];
    double array[input->n];

    // Veloity parameters
    double d12, d23, d31;
    double m12, m23, m31;
    double r1, r2, r3;
    double e1, e2, e3;
    double h1, h2, h3;
    double px, py, pz;
    double ln12, ln23, ln31;
    double x, y, z;
    double u, v, w;
    const double factor = 1 / (4 * PI);

    // Equations
    for (i = 0; i < input->n; i++) {

        // Freestream
        array[i] = -(input->freestream[0] * input->e3[i][0] + input->freestream[1] * input->e3[i][1] + input->freestream[2] * input->e3[i][2]);

        // Face j acting on i
        for (j = 0; j < input->n; j++) {

            // Velocity induced by face j on face i
            px = input->controlPoints[i][0] - input->facesCenters[j][0];
            py = input->controlPoints[i][1] - input->facesCenters[j][1];
            pz = input->controlPoints[i][2] - input->facesCenters[j][2];

            x = px * input->e1[j][0] + py * input->e1[j][1] + pz * input->e1[j][2];
            y = px * input->e2[j][0] + py * input->e2[j][1] + pz * input->e2[j][2];
            z = px * input->e3[j][0] + py * input->e3[j][1] + pz * input->e3[j][2];

            d12 = sqrt(pow(input->p2[j][0] - input->p1[j][0], 2) + pow(input->p2[j][1] - input->p1[j][1], 2));
            d23 = sqrt(pow(input->p3[j][0] - input->p2[j][0], 2) + pow(input->p3[j][1] - input->p2[j][1], 2));
            d31 = sqrt(pow(input->p1[j][0] - input->p3[j][0], 2) + pow(input->p1[j][1] - input->p3[j][1], 2));

            m12 = division(input->p2[j][1] - input->p1[j][1], input->p2[j][0] - input->p1[j][0]);
            m23 = division(input->p3[j][1] - input->p2[j][1], input->p3[j][0] - input->p2[j][0]);
            m31 = division(input->p3[j][1] - input->p1[j][1], input->p3[j][0] - input->p1[j][0]);

            r1 = sqrt(pow(x - input->p1[j][0], 2) + pow(y - input->p1[j][1], 2) + pow(z, 2));
            r2 = sqrt(pow(x - input->p2[j][0], 2) + pow(y - input->p2[j][1], 2) + pow(z, 2));
            r3 = sqrt(pow(x - input->p3[j][0], 2) + pow(y - input->p3[j][1], 2) + pow(z, 2));

            e1 = pow(x - input->p1[j][0], 2) + pow(z, 2);
            e2 = pow(x - input->p2[j][0], 2) + pow(z, 2);
            e3 = pow(x - input->p3[j][0], 2) + pow(z, 2);

            h1 = (x - input->p1[j][0]) * (y - input->p1[j][1]);
            h2 = (x - input->p2[j][0]) * (y - input->p2[j][1]);
            h3 = (x - input->p3[j][0]) * (y - input->p3[j][1]);

            ln12 = log(division(r1 + r2 - d12, r1 + r2 + d12));
            ln23 = log(division(r2 + r3 - d23, r2 + r3 + d23));
            ln31 = log(division(r3 + r1 - d31, r3 + r1 + d31));

            u = -factor * ( division(input->p2[j][1] - input->p1[j][1], d12) * ln12 + division(input->p3[j][1] - input->p2[j][1], d23) * ln23 + division(input->p1[j][1] - input->p3[j][1], d31) * ln31 );
            v = -factor * ( division(input->p1[j][0] - input->p2[j][0], d12) * ln12 + division(input->p2[j][1] - input->p3[j][1], d23) * ln23 + division(input->p3[j][1] - input->p1[j][1], d31) * ln31 );
            w = -factor * ( atan(division(m12 * e1 - h1, z * r1)) - atan(division(m12 * e2 - h2, z * r2)) + atan(division(m23 * e2 - h2, z * r2)) - atan(division(m23 * e3 - h3, z * r3)) + atan(division(m31 * e3 - h3, z * r3)) - atan(division(m31 * e1 - h1, z * r1)) );

            matrix[i][j] = (u * input->e1[j][0] + v * input->e2[j][0] + w * input->e3[j][0]) * input->e3[i][0] + (u * input->e1[j][1] + v * input->e2[j][1] + w * input->e3[j][1]) * input->e3[i][1] + (u * input->e1[j][2] + v * input->e2[j][2] + w * input->e3[j][2]) * input->e3[i][2];
            
            velFieldMatrix[i][j][0] = u;
            velFieldMatrix[i][j][1] = v;
            velFieldMatrix[i][j][2] = w;
            
        }
    }

    // Solve linear system
    linear_system_solver(input->n, matrix, array, output->sigma);

}

/// Calculate the doublet distribution in order to stisfy the Kutta Condition.
void doublet() {}

/// Calculate the velocity field
void velocity(struct input_struct *input, struct output_struct *output, double sourceVelFieldMatrix[input->n][input->n][3]) {

    int i, j;

    double vx, vy, vz;

    for (i = 0; i < input->n; i++) {

        vx = 0;
        vy = 0;
        vz = 0;

        for (j = 0; j < input->n; j++) {
            vx += sourceVelFieldMatrix[i][j][0] * output->sigma[j];
            vy += sourceVelFieldMatrix[i][j][0] * output->sigma[j];
            vz += sourceVelFieldMatrix[i][j][0] * output->sigma[j];
        }

        output->velocityNorm[i] = sqrt(pow(vx, 2) + pow(vy, 2) + pow(vz, 2));
        output->velocityField[i][0] = vx + input->freestream[0];
        output->velocityField[i][1] = vy + input->freestream[1];
        output->velocityField[i][2] = vz + input->freestream[2];

    }

}

/// Calculate the pressure distribution over the bird's surface
void pressure(struct input_struct *input, struct output_struct *output) {

    int i;
    const double rho05 = 0.5 * input->density;
    const double constant = input->pressure + rho05 * (pow(input->freestream[0], 2) + pow(input->freestream[1], 2) + pow(input->freestream[2], 2));

    for (i = 0; i < input->n; i++) {
        output->pressure[i] = constant - rho05 * output->velocityNorm[i];
    }

}

/// Calculate the source and doublet distribution over the bird's surface
void solve(struct input_struct *input, struct output_struct *output) {

    double sourceVelFieldMatrix[input->n][input->n][3];

    source(input, output, sourceVelFieldMatrix);
    // doublet()

    velocity(input, output, sourceVelFieldMatrix);
    pressure(input, output);
}

// void solve(struct input_struct *input, struct output_struct *output) {

//     int i;

//     for (i = 0; i < input->n; i++) {
//         output->sigma[i] = -(input->freestream[0] * input->e3[i][0] + input->freestream[1] * input->e3[i][1] + input->freestream[2] * input->e3[i][2]);
//     }

// }

void create_source_linear_system(struct input_struct *input, double **matrix, double *array, double ***velFieldMatrix) {

    // For loop
    int i, j;

    // Veloity parameters
    double d12, d23, d31;
    double m12, m23, m31;
    double r1, r2, r3;
    double e1, e2, e3;
    double h1, h2, h3;
    double px, py, pz;
    double ln12, ln23, ln31;
    double x, y, z;
    double u, v, w;
    const double factor = 1 / (4 * PI);

    // Equations
    for (i = 0; i < input->n; i++) {

        // Freestream
        array[i] = -(input->freestream[0] * input->e3[i][0] + input->freestream[1] * input->e3[i][1] + input->freestream[2] * input->e3[i][2]);

        // Face j acting on i
        for (j = 0; j < input->n; j++) {

            // Velocity induced by face j on face i
            px = input->controlPoints[i][0] - input->facesCenters[j][0];
            py = input->controlPoints[i][1] - input->facesCenters[j][1];
            pz = input->controlPoints[i][2] - input->facesCenters[j][2];

            x = px * input->e1[j][0] + py * input->e1[j][1] + pz * input->e1[j][2];
            y = px * input->e2[j][0] + py * input->e2[j][1] + pz * input->e2[j][2];
            z = px * input->e3[j][0] + py * input->e3[j][1] + pz * input->e3[j][2];

            d12 = sqrt(pow(input->p2[j][0] - input->p1[j][0], 2) + pow(input->p2[j][1] - input->p1[j][1], 2));
            d23 = sqrt(pow(input->p3[j][0] - input->p2[j][0], 2) + pow(input->p3[j][1] - input->p2[j][1], 2));
            d31 = sqrt(pow(input->p1[j][0] - input->p3[j][0], 2) + pow(input->p1[j][1] - input->p3[j][1], 2));

            m12 = division(input->p2[j][1] - input->p1[j][1], input->p2[j][0] - input->p1[j][0]);
            m23 = division(input->p3[j][1] - input->p2[j][1], input->p3[j][0] - input->p2[j][0]);
            m31 = division(input->p3[j][1] - input->p1[j][1], input->p3[j][0] - input->p1[j][0]);

            r1 = sqrt(pow(x - input->p1[j][0], 2) + pow(y - input->p1[j][1], 2) + pow(z, 2));
            r2 = sqrt(pow(x - input->p2[j][0], 2) + pow(y - input->p2[j][1], 2) + pow(z, 2));
            r3 = sqrt(pow(x - input->p3[j][0], 2) + pow(y - input->p3[j][1], 2) + pow(z, 2));

            e1 = pow(x - input->p1[j][0], 2) + pow(z, 2);
            e2 = pow(x - input->p2[j][0], 2) + pow(z, 2);
            e3 = pow(x - input->p3[j][0], 2) + pow(z, 2);

            h1 = (x - input->p1[j][0]) * (y - input->p1[j][1]);
            h2 = (x - input->p2[j][0]) * (y - input->p2[j][1]);
            h3 = (x - input->p3[j][0]) * (y - input->p3[j][1]);

            ln12 = log(division(r1 + r2 - d12, r1 + r2 + d12));
            ln23 = log(division(r2 + r3 - d23, r2 + r3 + d23));
            ln31 = log(division(r3 + r1 - d31, r3 + r1 + d31));

            u = -factor * ( division(input->p2[j][1] - input->p1[j][1], d12) * ln12 + division(input->p3[j][1] - input->p2[j][1], d23) * ln23 + division(input->p1[j][1] - input->p3[j][1], d31) * ln31 );
            v = -factor * ( division(input->p1[j][0] - input->p2[j][0], d12) * ln12 + division(input->p2[j][1] - input->p3[j][1], d23) * ln23 + division(input->p3[j][1] - input->p1[j][1], d31) * ln31 );
            w = -factor * ( atan(division(m12 * e1 - h1, z * r1)) - atan(division(m12 * e2 - h2, z * r2)) + atan(division(m23 * e2 - h2, z * r2)) - atan(division(m23 * e3 - h3, z * r3)) + atan(division(m31 * e3 - h3, z * r3)) - atan(division(m31 * e1 - h1, z * r1)) );

            matrix[i][j] = (u * input->e1[j][0] + v * input->e2[j][0] + w * input->e3[j][0]) * input->e3[i][0] + (u * input->e1[j][1] + v * input->e2[j][1] + w * input->e3[j][1]) * input->e3[i][1] + (u * input->e1[j][2] + v * input->e2[j][2] + w * input->e3[j][2]) * input->e3[i][2];
            
            velFieldMatrix[i][j][0] = u;
            velFieldMatrix[i][j][1] = v;
            velFieldMatrix[i][j][2] = w;
            
        }
    }

}