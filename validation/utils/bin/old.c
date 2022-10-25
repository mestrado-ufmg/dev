#include <math.h>
#include <stdio.h>
#include <lapacke.h>

/*
#####################################################
    CONSTANTS
#####################################################
*/
const double ZERO_ERROR = 1e-12;
const double PI = 3.14159265359;
const double FACTOR = 1 / (4 * PI);

/*
#####################################################
    STRUCTS
#####################################################
*/
struct Point {
    double x;
    double y;
    double z;
};

struct VerticeConnection {
    int n;
    double *coeffs;
    int *faces;
};

struct FacesConnection {
    int n;
    int *faces;
};

/*
#####################################################
    MATH FUNCTIONS
#####################################################
*/
double division(double a, double b) {
    if (b < 0) {
        return a / (b - ZERO_ERROR);
    } else {
        return a / (b + ZERO_ERROR);
    }
}

double norm(struct Point p) {
    return sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
}

struct Point cross(struct Point p1, struct Point p2) {
    struct Point p = {p1.y * p2.z - p1.z * p2.y, p1.z * p2.x - p1.x * p2.z, p1.x * p2.y - p1.y * p2.x};
    return p;
}

double dot(struct Point p1, struct Point p2) {
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

double angleBetweenVectors(struct Point p1, struct Point p2) {

    double norm1 = norm(p1);
    double norm2 = norm(p2);
    double dot12 = dot(p1, p2);

    return acos(dot12 / (norm1 * norm2));
}

/*
#####################################################
    HELPER FUNCTIONS
#####################################################
*/
void calculateVerticesConnection(int nv, int nf, double *vertices, int *faces, struct VerticeConnection *connection) {

    /* Parameters */
    int i, j, k, n, faceLine, verticeLine1, verticeLine2, verticeLine3;
    int *facesIds;
    double *angles;
    struct Point point1, point2;
    double angle;
    double sum;

    /* Initialize */
    facesIds = (int *)malloc(50 * sizeof(int));
    angles = (double *)malloc(50 * sizeof(double));

    /* Loop over vertices */
    for (i = 0; i < nv; i++) {

        // Reset the number of faces and angles
        n = 0;
        sum = 0.0;

        /* Loop over faces */
        for (j = 0; j < nf; j++) {

            faceLine = j * 3;

            /* Check if the face contain the vertice */
            if ((faces[faceLine] == i) || (faces[faceLine + 1] == i) || (faces[faceLine + 2] == i)) {

                /* Calculate the angle */
                if (faces[faceLine] == i)
                {
                    verticeLine1 = 3 * faces[faceLine];
                    verticeLine2 = 3 * faces[faceLine + 1];
                    verticeLine3 = 3 * faces[faceLine + 2];
                }
                else if (faces[faceLine + 1] == i)
                {
                    verticeLine3 = 3 * faces[faceLine];
                    verticeLine1 = 3 * faces[faceLine + 1];
                    verticeLine2 = 3 * faces[faceLine + 2];
                }
                else
                {
                    verticeLine2 = 3 * faces[faceLine];
                    verticeLine3 = 3 * faces[faceLine + 1];
                    verticeLine1 = 3 * faces[faceLine + 2];
                }

                point1.x = vertices[verticeLine2] - vertices[verticeLine1];
                point1.y = vertices[verticeLine2 + 1] - vertices[verticeLine1 + 1];
                point1.z = vertices[verticeLine2 + 2] - vertices[verticeLine1 + 2];

                point2.x = vertices[verticeLine3] - vertices[verticeLine1];
                point2.y = vertices[verticeLine3 + 1] - vertices[verticeLine1 + 1];
                point2.z = vertices[verticeLine3 + 2] - vertices[verticeLine1 + 2];

                angle = angleBetweenVectors(point1, point2);

                /* Increase the number of faces and save the angle */
                angles[n] = angle;
                facesIds[n] = j;
                sum = sum + angle;

                // Increase the number of faces
                n++;
            }
        }

        /* Save the number of faces, faces id's and angles in connection */
        connection[i].n = n;
        connection[i].coeffs = (double *)malloc(n * sizeof(double));
        connection[i].faces = (int *)malloc(n * sizeof(int));

        for (k = 0; k < n; k++)
        {
            connection[i].coeffs[k] = angles[k] / sum;
            connection[i].faces[k] = facesIds[k];
        }
    }

    free(facesIds);
    free(angles);
}

void calculateFacesConnection(int nv, int nf, int *faces, struct VerticeConnection *verticesConnection, struct FacesConnection *facesConnection) {

    /* Parameters */
    int i, j, k;
    int check;
    int index;
    int *connectedFaces;
    int n;

    /* Initialize */
    connectedFaces = (int *)malloc(50 * sizeof(int));

    for (i = 0; i < nf; i++)
    {

        n = 0;

        // First index
        index = faces[3 * i];

        for (j = 0; j < verticesConnection[index].n; j++)
        {
            if (verticesConnection[index].faces[j] != i)
            {
                connectedFaces[n] = verticesConnection[index].faces[j];
                n = n + 1;
            }
        }

        // Second index
        index = faces[3 * i + 1];

        for (j = 0; j < verticesConnection[index].n; j++)
        {

            check = 0;

            for (k = 0; k < n; k++)
            {
                if (connectedFaces[k] == verticesConnection[index].faces[j])
                {
                    check = 1;
                    break;
                }
            }

            if ((verticesConnection[index].faces[j] != i) && (check == 0))
            {
                connectedFaces[n] = verticesConnection[index].faces[j];
                n = n + 1;
            }
        }

        // Second index
        index = faces[3 * i + 2];

        for (j = 0; j < verticesConnection[index].n; j++)
        {

            check = 0;

            for (k = 0; k < n; k++)
            {
                if (connectedFaces[k] == verticesConnection[index].faces[j])
                {
                    check = 1;
                    break;
                }
            }

            if ((verticesConnection[index].faces[j] != i) && (check == 0))
            {
                connectedFaces[n] = verticesConnection[index].faces[j];
                n = n + 1;
            }
        }

        facesConnection[i].n = n;
        facesConnection[i].faces = (int *)malloc(n * sizeof(int));

        for (j = 0; j < n; j++)
        {
            facesConnection[i].faces[j] = connectedFaces[j];
        }
    }

    free(connectedFaces);
}

/*
#####################################################
    LINEAR SYSTEM
#####################################################
*/
void solveLinearSystem(lapack_int n, double *A, double *b) {

    /* Locals */
    lapack_int info;

    /* Local arrays */
    lapack_int *ipiv;

    /* Initialization */
    ipiv = (lapack_int *)malloc(n * sizeof(lapack_int));

    /* Solve the equations A*X = B */
    info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, 1, A, n, ipiv, b, 1);

    /* Free */
    free(ipiv);
}

/*
#####################################################
    POTENTIAL FLOW
#####################################################
*/
void sourceFunc(struct Point p, struct Point p1, struct Point p2, struct Point p3, struct Point e1, struct Point e2, struct Point e3, double area, double maxDistance, double *vel) {

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

void lineFunc(struct Point p, struct Point p1, struct Point p2, double *vel) {

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

void doubletFunc(struct Point p, struct Point p1, struct Point p2, struct Point p3, struct Point e1, struct Point e2, struct Point e3, double area, double maxDistance, double *vel) {

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

void createLinearSystem(int n, double *facesAreas, double *facesMaxDistance, double *facesCenter, double *controlPoints, double *p1, double *p2, double *p3, double *e1, double *e2, double *e3, double *freestream, double *sigma, int nSpanWake, int nWake, int *wakeGrid, double *wakeVertices, int *wakeFaces, double *matrix, double *array, double *matrixVelx, double *matrixVely, double *matrixVelz, double *arrayVel) {

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

    ///---------------------------------------///
    /// Linear system

    // Create
    for (i = 0; i < n; i++) {

        i3D1 = i * 3;
        i3D2 = i3D1 + 1;
        i3D3 = i3D1 + 2;

        e3iPoint.x = e3[i3D1];
        e3iPoint.y = e3[i3D2];
        e3iPoint.z = e3[i3D3];

        array[i] = 0;
        arrayVel[i * 3] = 0;
        arrayVel[i * 3 + 1] = 0;
        arrayVel[i * 3 + 2] = 0;

        // Surface
        // Effect of j on i
        for (j = 0; j < n; j++)
        {

            j3D1 = j * 3;
            j3D2 = j3D1 + 1;
            j3D3 = j3D1 + 2;

            j2D1 = j * 2;
            j2D2 = j2D1 + 1;

            // Points
            e1jPoint.x = e1[j3D1];
            e1jPoint.y = e1[j3D2];
            e1jPoint.z = e1[j3D3];
            e2jPoint.x = e2[j3D1];
            e2jPoint.y = e2[j3D2];
            e2jPoint.z = e2[j3D3];
            e3jPoint.x = e3[j3D1];
            e3jPoint.y = e3[j3D2];
            e3jPoint.z = e3[j3D3];

            p.x = controlPoints[i3D1] - facesCenter[j3D1];
            p.y = controlPoints[i3D2] - facesCenter[j3D2];
            p.z = controlPoints[i3D3] - facesCenter[j3D3];

            pLocal.x = p.x * e1jPoint.x + p.y * e1jPoint.y + p.z * e1jPoint.z;
            pLocal.y = p.x * e2jPoint.x + p.y * e2jPoint.y + p.z * e2jPoint.z;
            pLocal.z = p.x * e3jPoint.x + p.y * e3jPoint.y + p.z * e3jPoint.z;

            p1Local.x = p1[j2D1];
            p1Local.y = p1[j2D2];
            p2Local.x = p2[j2D1];
            p2Local.y = p2[j2D2];
            p3Local.x = p3[j2D1];
            p3Local.y = p3[j2D2];

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

                    matrix[i * n + wakeFaces[k * 2]] = matrix[i * n + wakeFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + wakeFaces[k * 2 + 1]] = matrix[i * n + wakeFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + wakeFaces[k * 2]] = matrixVelx[i * n + wakeFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + wakeFaces[k * 2]] = matrixVely[i * n + wakeFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + wakeFaces[k * 2]] = matrixVelz[i * n + wakeFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + wakeFaces[k * 2 + 1]] = matrixVelx[i * n + wakeFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + wakeFaces[k * 2 + 1]] = matrixVely[i * n + wakeFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + wakeFaces[k * 2 + 1]] = matrixVelz[i * n + wakeFaces[k * 2 + 1]] + lineVel[2];
                } else if (k == nSpanWake - 1) {

                    matrix[i * n + wakeFaces[(k - 1) * 2]] = matrix[i * n + wakeFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + wakeFaces[(k - 1) * 2 + 1]] = matrix[i * n + wakeFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + wakeFaces[(k - 1) * 2]] = matrixVelx[i * n + wakeFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + wakeFaces[(k - 1) * 2]] = matrixVely[i * n + wakeFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + wakeFaces[(k - 1) * 2]] = matrixVelz[i * n + wakeFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + wakeFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + wakeFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + wakeFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + wakeFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + wakeFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + wakeFaces[(k - 1) * 2 + 1]] - lineVel[2];
                } else {

                    matrix[i * n + wakeFaces[(k - 1) * 2]] = matrix[i * n + wakeFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + wakeFaces[(k - 1) * 2 + 1]] = matrix[i * n + wakeFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + wakeFaces[(k - 1) * 2]] = matrixVelx[i * n + wakeFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + wakeFaces[(k - 1) * 2]] = matrixVely[i * n + wakeFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + wakeFaces[(k - 1) * 2]] = matrixVelz[i * n + wakeFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + wakeFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + wakeFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + wakeFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + wakeFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + wakeFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + wakeFaces[(k - 1) * 2 + 1]] - lineVel[2];

                    matrix[i * n + wakeFaces[k * 2]] = matrix[i * n + wakeFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + wakeFaces[k * 2 + 1]] = matrix[i * n + wakeFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + wakeFaces[k * 2]] = matrixVelx[i * n + wakeFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + wakeFaces[k * 2]] = matrixVely[i * n + wakeFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + wakeFaces[k * 2]] = matrixVelz[i * n + wakeFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + wakeFaces[k * 2 + 1]] = matrixVelx[i * n + wakeFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + wakeFaces[k * 2 + 1]] = matrixVely[i * n + wakeFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + wakeFaces[k * 2 + 1]] = matrixVelz[i * n + wakeFaces[k * 2 + 1]] + lineVel[2];

                    indexLine1 = wakeGrid[(k - 1) * nWake] * 3;
                    p1Line.x = wakeVertices[indexLine1];
                    p1Line.y = wakeVertices[indexLine1 + 1];
                    p1Line.z = wakeVertices[indexLine1 + 2];

                    indexLine2 = wakeGrid[k * nWake] * 3;
                    p2Line.x = wakeVertices[indexLine2];
                    p2Line.y = wakeVertices[indexLine2 + 1];
                    p2Line.z = wakeVertices[indexLine2 + 2];

                    lineFunc(p, p1Line, p2Line, lineVel);

                    matrix[i * n + wakeFaces[(k - 1) * 2]] = matrix[i * n + wakeFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + wakeFaces[(k - 1) * 2 + 1]] = matrix[i * n + wakeFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + wakeFaces[(k - 1) * 2]] = matrixVelx[i * n + wakeFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + wakeFaces[(k - 1) * 2]] = matrixVely[i * n + wakeFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + wakeFaces[(k - 1) * 2]] = matrixVelz[i * n + wakeFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + wakeFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + wakeFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + wakeFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + wakeFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + wakeFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + wakeFaces[(k - 1) * 2 + 1]] - lineVel[2];
                }
            }
        }

    }

    free(sourceVel);
    free(doubletVel);
    free(lineVel);
}

void calculateDoubletDistribution(int n, double *A, double *b, double *transpiration, double *sol) {

    /* Loop parameter */
    int i;

    /* Copy the input */
    double *aux = (double *)malloc(n * n * sizeof(double));

    for (i = 0; i < n * n; i++)
        aux[i] = A[i];
    for (i = 0; i < n; i++)
        sol[i] = b[i] + transpiration[i];

    /* Solve the linear system */
    solveLinearSystem(n, aux, sol);

    free(aux);
}

void calculateSurfaceParameters(int n, double *matrixVelx, double *matrixVely, double *matrixVelz, double *arrayVel, double *doublet, double freestream, double *velx, double *vely, double *velz, double *velNorm, double *cp, double *mach, double sound_speed) {

    /* Loop parameter */
    int i, j, line1, line2, point;

    /* Calculate parameters */
    for (i = 0; i < n; i++) {

        line1 = i * 3;

        // Set velocity equal 0
        velx[i] = arrayVel[line1];
        vely[i] = arrayVel[line1 + 1];
        velz[i] = arrayVel[line1 + 2];

        // Current line
        line2 = i * n;

        for (j = 0; j < n; j++) {

            point = line2 + j;

            velx[i] = velx[i] + matrixVelx[point] * doublet[j];
            vely[i] = vely[i] + matrixVely[point] * doublet[j];
            velz[i] = velz[i] + matrixVelz[point] * doublet[j];
        }

        // Velocity norm
        velNorm[i] = sqrt(velx[i] * velx[i] + vely[i] * vely[i] + velz[i] * velz[i]);

        // Pressure coefficient
        cp[i] = 1 - pow(velNorm[i] / freestream, 2);

        // Mach
        mach[i] = velNorm[i] / sound_speed;
    }
}

/*
#####################################################
    SOLVER
#####################################################
*/
void solve(int type,
           int nv,
           int nf,
           double *vertices,
           int *faces,
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
           int nSpanWake,
           int nWake,
           int *wakeGrid,
           double *wakeVertices,
           int *wakeFaces,
           double *doublet,
           double *velx, double *vely, double *velz, double *velNorm,
           double *cp,
           double *mach,
           double *delta,
           double *A, double *B,
           double *Psi,
           double *Ctau1, double *Ctau2,
           double *tau_x, double *tau_y, double *tau_z,
           double density, double viscosity, double sound_speed,
           double *transpiration,
           double *sigma_v,
           double *doublet_v,
           double *cp_v,
           double *velx_v, double *vely_v, double *velz_v, double *velNorm_v,
           double *transpiration_v,
           double *delta_v,
           double *A_v, double *B_v,
           double *Psi_v,
           double *Ctau1_v, double *Ctau2_v,
           double *tau_x_v, double *tau_y_v, double *tau_z_v) {

    printf("  > Potential flow\n");

    /* Potential flow parameters */
    double *matrix, *array;
    double *matrixVelx, *matrixVely, *matrixVelz, *arrayVel;
    double freestreamNorm;

    /* Initialize */
    matrix = (double *)malloc(nf * nf * sizeof(double));
    array = (double *)malloc(nf * sizeof(double));
    matrixVelx = (double *)malloc(nf * nf * sizeof(double));
    matrixVely = (double *)malloc(nf * nf * sizeof(double));
    matrixVelz = (double *)malloc(nf * nf * sizeof(double));
    arrayVel = (double *)malloc(nf * 3 * sizeof(double));
    freestreamNorm = sqrt(freestream[0] * freestream[0] + freestream[1] * freestream[1] + freestream[2] * freestream[2]);

    /* Potential flow system */
    printf("    * Creating linear system\n");
    createLinearSystem(nf, facesAreas, facesMaxDistance, facesCenter, controlPoints, p1, p2, p3, e1, e2, e3, freestream, sigma, nSpanWake, nWake, wakeGrid, wakeVertices, wakeFaces, matrix, array, matrixVelx, matrixVely, matrixVelz, arrayVel);

    /* Solve linear system with zero transpiration */
    printf("    * Solving linear system\n");
    calculateDoubletDistribution(nf, matrix, array, transpiration, doublet);

    /* Calculate potential surface parameters */
    calculateSurfaceParameters(nf, matrixVelx, matrixVely, matrixVelz, arrayVel, doublet, freestreamNorm, velx, vely, velz, velNorm, cp, mach, sound_speed);

    /* Vertices values */
    int i, j;
    struct VerticeConnection *vertices_connection;

    /* Initialize */
    vertices_connection = (struct VerticeConnection*)malloc(nv * sizeof(struct VerticeConnection));
    calculateVerticesConnection(nv, nf, vertices, faces, vertices_connection);

    /* Transpiration */
    for (i = 0; i < nf; i++) transpiration[i] = e3[3 * i] * velx[i] + e3[3 * i + 1] * vely[i] + e3[3 * i + 2] * velz[i];
    
    /* Vertices values */
    for (i = 0; i < nv; i++) {

        sigma_v[i] = 0.0;
        doublet_v[i] = 0.0;
        cp_v[i] = 0.0;
        velx_v[i] = 0.0;
        vely_v[i] = 0.0;
        velz_v[i] = 0.0;
        velNorm_v[i] = 0.0;
        transpiration_v[i] = 0.0;
        delta_v[i] = 0.0;
        A_v[i] = 0.0;
        B_v[i] = 0.0;
        Psi_v[i] = 0.0;
        Ctau1_v[i] = 0.0;
        Ctau2_v[i] = 0.0;
        tau_x_v[i] = 0.0;
        tau_y_v[i] = 0.0;
        tau_z_v[i] = 0.0;

        for (j = 0; j < vertices_connection[i].n; j++) {
            sigma_v[i] = sigma_v[i] + sigma[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            doublet_v[i] = doublet_v[i] + doublet[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            cp_v[i] = cp_v[i] + cp[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            velx_v[i] = velx_v[i] + velx[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            vely_v[i] = vely_v[i] + vely[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            velz_v[i] = velz_v[i] + velz[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            velNorm_v[i] = velNorm_v[i] + velNorm[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            transpiration_v[i] = transpiration_v[i] + transpiration[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
        }
    
    }

    free(matrix);
    free(array);
    free(matrixVelx);
    free(matrixVely);
    free(matrixVelz);
    free(arrayVel);
    free(vertices_connection);
}