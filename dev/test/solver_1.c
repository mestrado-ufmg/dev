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
const int LAYERS = 300;

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
}

/*
#####################################################
    POTENTIAL FLOW
#####################################################
*/
void sourceFunc(struct Point p, struct Point p1, struct Point p2, struct Point p3, struct Point e1, struct Point e2, struct Point e3, double area, double maxDistance, double *vel) {

    double u, v, w;
    double distance = norm(p);

    if (distance > maxDistance)
    {

        double pNorm3 = pow(norm(p), 3);

        u = FACTOR * area * p.x / pNorm3;
        v = FACTOR * area * p.y / pNorm3;
        w = FACTOR * area * p.z / pNorm3;
    }
    else
    {

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

    if (distance > maxDistance)
    {

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
    }
    else
    {

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

void createLinearSystem(int n, double *facesAreas, double *facesMaxDistance, double *facesCenter, double *controlPoints, double *p1, double *p2, double *p3, double *e1, double *e2, double *e3, double *freestream, double *sigma, int nSpanLeftWing, int nWakeLeftWing, int *leftWingGrid, double *leftWingVertices, int *leftWingFaces, int nSpanRightWing, int nWakeRightWing, int *rightWingGrid, double *rightWingVertices, int *rightWingFaces, int nSpanTail, int nWakeTail, int *tailGrid, double *tailVertices, int *tailFaces, double *matrix, double *array, double *matrixVelx, double *matrixVely, double *matrixVelz, double *arrayVel) {

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
    for (i = 0; i < n; i++)
    {

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
        for (k = 0; k < nSpanLeftWing; k++)
        {

            for (l = 0; l < nWakeLeftWing - 1; l++)
            {

                indexLine1 = leftWingGrid[k * nWakeLeftWing + l] * 3;
                p1Line.x = leftWingVertices[indexLine1];
                p1Line.y = leftWingVertices[indexLine1 + 1];
                p1Line.z = leftWingVertices[indexLine1 + 2];

                indexLine2 = leftWingGrid[k * nWakeLeftWing + l + 1] * 3;
                p2Line.x = leftWingVertices[indexLine2];
                p2Line.y = leftWingVertices[indexLine2 + 1];
                p2Line.z = leftWingVertices[indexLine2 + 2];

                lineFunc(p, p1Line, p2Line, lineVel);

                if (k == 0)
                {

                    matrix[i * n + leftWingFaces[k * 2]] = matrix[i * n + leftWingFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[k * 2 + 1]] = matrix[i * n + leftWingFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + leftWingFaces[k * 2]] = matrixVelx[i * n + leftWingFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + leftWingFaces[k * 2]] = matrixVely[i * n + leftWingFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + leftWingFaces[k * 2]] = matrixVelz[i * n + leftWingFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + leftWingFaces[k * 2 + 1]] = matrixVelx[i * n + leftWingFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + leftWingFaces[k * 2 + 1]] = matrixVely[i * n + leftWingFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + leftWingFaces[k * 2 + 1]] = matrixVelz[i * n + leftWingFaces[k * 2 + 1]] + lineVel[2];
                }
                else if (k == nSpanLeftWing - 1)
                {

                    matrix[i * n + leftWingFaces[(k - 1) * 2]] = matrix[i * n + leftWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + leftWingFaces[(k - 1) * 2]] = matrixVelx[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + leftWingFaces[(k - 1) * 2]] = matrixVely[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + leftWingFaces[(k - 1) * 2]] = matrixVelz[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[2];
                }
                else
                {

                    matrix[i * n + leftWingFaces[(k - 1) * 2]] = matrix[i * n + leftWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + leftWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + leftWingFaces[(k - 1) * 2]] = matrixVelx[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + leftWingFaces[(k - 1) * 2]] = matrixVely[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + leftWingFaces[(k - 1) * 2]] = matrixVelz[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[2];

                    matrix[i * n + leftWingFaces[k * 2]] = matrix[i * n + leftWingFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + leftWingFaces[k * 2 + 1]] = matrix[i * n + leftWingFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + leftWingFaces[k * 2]] = matrixVelx[i * n + leftWingFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + leftWingFaces[k * 2]] = matrixVely[i * n + leftWingFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + leftWingFaces[k * 2]] = matrixVelz[i * n + leftWingFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + leftWingFaces[k * 2 + 1]] = matrixVelx[i * n + leftWingFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + leftWingFaces[k * 2 + 1]] = matrixVely[i * n + leftWingFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + leftWingFaces[k * 2 + 1]] = matrixVelz[i * n + leftWingFaces[k * 2 + 1]] + lineVel[2];

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

                    matrixVelx[i * n + leftWingFaces[(k - 1) * 2]] = matrixVelx[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + leftWingFaces[(k - 1) * 2]] = matrixVely[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + leftWingFaces[(k - 1) * 2]] = matrixVelz[i * n + leftWingFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + leftWingFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + leftWingFaces[(k - 1) * 2 + 1]] - lineVel[2];
                }
            }
        }

        // Right wing
        for (k = 0; k < nSpanRightWing; k++)
        {

            for (l = 0; l < nWakeRightWing - 1; l++)
            {

                indexLine1 = rightWingGrid[k * nWakeRightWing + l] * 3;
                p1Line.x = rightWingVertices[indexLine1];
                p1Line.y = rightWingVertices[indexLine1 + 1];
                p1Line.z = rightWingVertices[indexLine1 + 2];

                indexLine2 = rightWingGrid[k * nWakeRightWing + l + 1] * 3;
                p2Line.x = rightWingVertices[indexLine2];
                p2Line.y = rightWingVertices[indexLine2 + 1];
                p2Line.z = rightWingVertices[indexLine2 + 2];

                lineFunc(p, p1Line, p2Line, lineVel);

                if (k == 0)
                {

                    matrix[i * n + rightWingFaces[k * 2]] = matrix[i * n + rightWingFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[k * 2 + 1]] = matrix[i * n + rightWingFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + rightWingFaces[k * 2]] = matrixVelx[i * n + rightWingFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + rightWingFaces[k * 2]] = matrixVely[i * n + rightWingFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + rightWingFaces[k * 2]] = matrixVelz[i * n + rightWingFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + rightWingFaces[k * 2 + 1]] = matrixVelx[i * n + rightWingFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + rightWingFaces[k * 2 + 1]] = matrixVely[i * n + rightWingFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + rightWingFaces[k * 2 + 1]] = matrixVelz[i * n + rightWingFaces[k * 2 + 1]] + lineVel[2];
                }
                else if (k == nSpanRightWing - 1)
                {

                    matrix[i * n + rightWingFaces[(k - 1) * 2]] = matrix[i * n + rightWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + rightWingFaces[(k - 1) * 2]] = matrixVelx[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + rightWingFaces[(k - 1) * 2]] = matrixVely[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + rightWingFaces[(k - 1) * 2]] = matrixVelz[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[2];
                }
                else
                {

                    matrix[i * n + rightWingFaces[(k - 1) * 2]] = matrix[i * n + rightWingFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrix[i * n + rightWingFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + rightWingFaces[(k - 1) * 2]] = matrixVelx[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + rightWingFaces[(k - 1) * 2]] = matrixVely[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + rightWingFaces[(k - 1) * 2]] = matrixVelz[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[2];

                    matrix[i * n + rightWingFaces[k * 2]] = matrix[i * n + rightWingFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + rightWingFaces[k * 2 + 1]] = matrix[i * n + rightWingFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + rightWingFaces[k * 2]] = matrixVelx[i * n + rightWingFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + rightWingFaces[k * 2]] = matrixVely[i * n + rightWingFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + rightWingFaces[k * 2]] = matrixVelz[i * n + rightWingFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + rightWingFaces[k * 2 + 1]] = matrixVelx[i * n + rightWingFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + rightWingFaces[k * 2 + 1]] = matrixVely[i * n + rightWingFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + rightWingFaces[k * 2 + 1]] = matrixVelz[i * n + rightWingFaces[k * 2 + 1]] + lineVel[2];

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

                    matrixVelx[i * n + rightWingFaces[(k - 1) * 2]] = matrixVelx[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + rightWingFaces[(k - 1) * 2]] = matrixVely[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + rightWingFaces[(k - 1) * 2]] = matrixVelz[i * n + rightWingFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + rightWingFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + rightWingFaces[(k - 1) * 2 + 1]] - lineVel[2];
                }
            }
        }

        // Tail
        for (k = 0; k < nSpanTail; k++)
        {

            for (l = 0; l < nWakeTail - 1; l++)
            {

                indexLine1 = tailGrid[k * nWakeTail + l] * 3;
                p1Line.x = tailVertices[indexLine1];
                p1Line.y = tailVertices[indexLine1 + 1];
                p1Line.z = tailVertices[indexLine1 + 2];

                indexLine2 = tailGrid[k * nWakeTail + l + 1] * 3;
                p2Line.x = tailVertices[indexLine2];
                p2Line.y = tailVertices[indexLine2 + 1];
                p2Line.z = tailVertices[indexLine2 + 2];

                lineFunc(p, p1Line, p2Line, lineVel);

                if (k == 0)
                {

                    matrix[i * n + tailFaces[k * 2]] = matrix[i * n + tailFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[k * 2 + 1]] = matrix[i * n + tailFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + tailFaces[k * 2]] = matrixVelx[i * n + tailFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + tailFaces[k * 2]] = matrixVely[i * n + tailFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + tailFaces[k * 2]] = matrixVelz[i * n + tailFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + tailFaces[k * 2 + 1]] = matrixVelx[i * n + tailFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + tailFaces[k * 2 + 1]] = matrixVely[i * n + tailFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + tailFaces[k * 2 + 1]] = matrixVelz[i * n + tailFaces[k * 2 + 1]] + lineVel[2];
                }
                else if (k == nSpanTail - 1)
                {

                    matrix[i * n + tailFaces[(k - 1) * 2]] = matrix[i * n + tailFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[(k - 1) * 2 + 1]] = matrix[i * n + tailFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + tailFaces[(k - 1) * 2]] = matrixVelx[i * n + tailFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + tailFaces[(k - 1) * 2]] = matrixVely[i * n + tailFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + tailFaces[(k - 1) * 2]] = matrixVelz[i * n + tailFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[2];
                }
                else
                {

                    matrix[i * n + tailFaces[(k - 1) * 2]] = matrix[i * n + tailFaces[(k - 1) * 2]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[(k - 1) * 2 + 1]] = matrix[i * n + tailFaces[(k - 1) * 2 + 1]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + tailFaces[(k - 1) * 2]] = matrixVelx[i * n + tailFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + tailFaces[(k - 1) * 2]] = matrixVely[i * n + tailFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + tailFaces[(k - 1) * 2]] = matrixVelz[i * n + tailFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[2];

                    matrix[i * n + tailFaces[k * 2]] = matrix[i * n + tailFaces[k * 2]] - (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);
                    matrix[i * n + tailFaces[k * 2 + 1]] = matrix[i * n + tailFaces[k * 2 + 1]] + (e3iPoint.x * lineVel[0] + e3iPoint.y * lineVel[1] + e3iPoint.z * lineVel[2]);

                    matrixVelx[i * n + tailFaces[k * 2]] = matrixVelx[i * n + tailFaces[k * 2]] - lineVel[0];
                    matrixVely[i * n + tailFaces[k * 2]] = matrixVely[i * n + tailFaces[k * 2]] - lineVel[1];
                    matrixVelz[i * n + tailFaces[k * 2]] = matrixVelz[i * n + tailFaces[k * 2]] - lineVel[2];

                    matrixVelx[i * n + tailFaces[k * 2 + 1]] = matrixVelx[i * n + tailFaces[k * 2 + 1]] + lineVel[0];
                    matrixVely[i * n + tailFaces[k * 2 + 1]] = matrixVely[i * n + tailFaces[k * 2 + 1]] + lineVel[1];
                    matrixVelz[i * n + tailFaces[k * 2 + 1]] = matrixVelz[i * n + tailFaces[k * 2 + 1]] + lineVel[2];

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

                    matrixVelx[i * n + tailFaces[(k - 1) * 2]] = matrixVelx[i * n + tailFaces[(k - 1) * 2]] + lineVel[0];
                    matrixVely[i * n + tailFaces[(k - 1) * 2]] = matrixVely[i * n + tailFaces[(k - 1) * 2]] + lineVel[1];
                    matrixVelz[i * n + tailFaces[(k - 1) * 2]] = matrixVelz[i * n + tailFaces[(k - 1) * 2]] + lineVel[2];

                    matrixVelx[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVelx[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[0];
                    matrixVely[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVely[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[1];
                    matrixVelz[i * n + tailFaces[(k - 1) * 2 + 1]] = matrixVelz[i * n + tailFaces[(k - 1) * 2 + 1]] - lineVel[2];
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
    for (i = 0; i < n; i++)
    {

        line1 = i * 3;

        // Set velocity equal 0
        velx[i] = arrayVel[line1];
        vely[i] = arrayVel[line1 + 1];
        velz[i] = arrayVel[line1 + 2];

        // Current line
        line2 = i * n;

        for (j = 0; j < n; j++)
        {

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
    BOUNDARY LAYER
#####################################################
*/
struct IntegralThickness {
    double *delta_1_ast, *delta_2_ast;
    double *phi_11, *phi_12, *phi_21, *phi_22;
    double *Ktauxx, *Ktauxy, *Ktauyx, *Ktauyy;
    double *Sx, *Sy;
    double *phi_1_ast, *phi_2_ast;
    double *theta_1_o, *theta_2_o;
    double *delta_1_line, *delta_2_line;
    double *delta_1_o, *delta_2_o;
    double *S0, *T0;
    double *C_D, *C_D_x, *C_D_o;
};

struct FaceDivergents {
    double div_phi_11_12, div_phi_21_22;
    double div_delta_1_2_ast;
    double div_phi_1_2_ast;
    double div_theta_1_2_o;
    double div_Ktau_xx_xy;
    double div_Ktau_yx_yy;
};

struct FaceGradients {
    double grad_u2_x, grad_u2_y;
    double grad_phi_x, grad_phi_y;
};

struct FaceEquations {
    double momentum_x, momentum_y;
    double kinetic_energy;
    double lateral_curvature;
    double shear_stress_x, shear_stress_y;
    double obj;
    int interaction;
};

double absValue(double a) {
    if (a < 0) {
        return - a;
    } else {
        return a;
    }
}

void trapz(int n, double *x, double *y, double *out) {
    *out = 0.0;
    for (int i = 0; i < n - 1; i++) *out = *out + 0.5 * (x[i + 1] - x[i]) * (y[i + 1] + y[i]);
}

void calculateParameters(int face, double delta, double A, double B, double Psi, double amp, double velocity, double density, double viscosity, double mach, struct IntegralThickness params) {
    
    /* Profiles */
    int i;
    double delta_eta = 1 / ((double) LAYERS - 1);
    double *eta = (double *)malloc(LAYERS * sizeof(double));
    double *U = (double *)malloc(LAYERS * sizeof(double));
    double *W = (double *)malloc(LAYERS * sizeof(double));
    double *dUdeta = (double *)malloc(LAYERS * sizeof(double));
    double *dWdeta = (double *)malloc(LAYERS * sizeof(double));
    double *S = (double *)malloc(LAYERS * sizeof(double));
    double *T = (double *)malloc(LAYERS * sizeof(double));
    double *R = (double *)malloc(LAYERS * sizeof(double));
    double *psi = (double *)malloc(LAYERS * sizeof(double));
    double *dpsideta = (double *)malloc(LAYERS * sizeof(double));

    double reynolds_delta = delta * velocity * density / viscosity;
    double ratio_mu;

    double f0, f1, f2, f3, f4;
    double df0deta, df1deta, df2deta, df3deta, df4deta;

    for (i = 0; i < LAYERS; i++) {

        // Height
        eta[i] = delta_eta * i;
        
        // Velocity
        f0 = 1 - 0.6 * (A - 3) * pow(eta[i], 3);
        f1 = 6 * pow(eta[i], 2) - 8 * pow(eta[i], 3) + 3 * pow(eta[i], 4);
        f2 = eta[i] - 3 * pow(eta[i], 2) + 3 * pow(eta[i], 3) - pow(eta[i], 4);
        f3 = (eta[i] - 4 * pow(eta[i], 2) + 6 * pow(eta[i], 3) - 4 * pow(eta[i], 4) + pow(eta[i], 5)) * pow(1 - eta[i], 2);
        f4 = (pow(eta[i], 2) - 3 * pow(eta[i], 3) + 3 * pow(eta[i], 4) - pow(eta[i], 5)) * pow(1 - eta[i], 2);

        U[i] = A * f0 * f2 + f1;
        W[i] = B * f3 + Psi * f4;

        // Angular deflection
        psi[i] = atan(division(W[i], U[i]));
        if (i == 0) {
            dpsideta[i] = 0.0;
        } else {
            dpsideta[i] = ((dUdeta[i] * W[i] - U[i] * dWdeta[i]) / (U[i] * U[i])) / (1 - pow(division(W[i], U[i]), 2));
        }

        // Density and viscosity ratio
        R[i] = 1 / (1 + 0.2 * mach * mach * (1 - pow(U[i], 2) - pow(W[i], 2)));
        ratio_mu = pow(1 / R[i], 1.5) * (1 + 1 / R[0]) / (1 / R[i] + 1 / R[0]);

        // Velocity gradient
        df0deta = - 1.8 * (A - 3) * pow(eta[i], 2);
        df1deta = 12 * eta[i] - 24 * pow(eta[i], 2) + 3 * pow(eta[i], 4);
        df2deta = 1 - 6 * eta[i] + 9 * pow(eta[i], 2) - 4 * pow(eta[i], 3);
        df3deta = (1 - 8 * eta[i] + 18 * pow(eta[i], 2) - 16 * pow(eta[i], 3) + 5 * pow(eta[i], 4)) * pow(1 - eta[i], 2) + (eta[i] - 4 * pow(eta[i], 2) + 6 * pow(eta[i], 3) - 4 * pow(eta[i], 4) + pow(eta[i], 5)) * 2 * (1 - eta[i]);
        df4deta = (2 * eta[i] - 9 * pow(eta[i], 2) + 12 * pow(eta[i], 3) - 5 * pow(eta[i], 4)) * pow(1 - eta[i], 2) + (pow(eta[i], 2) - 3 * pow(eta[i], 3) + 3 * pow(eta[i], 4) - pow(eta[i], 5)) * 2 * (1 - eta[i]);
        
        dUdeta[i] = A * (df0deta * f2 + f0 * df2deta) + df0deta;
        dWdeta[i] = B * df3deta + Psi * df4deta;
        
        // Shear stress
        S[i] = (1 / reynolds_delta) * ratio_mu * dUdeta[i];
        T[i] = (1 / reynolds_delta) * ratio_mu * dWdeta[i];

    }

    /* Integral thickness */
    double *func = (double *)malloc(LAYERS * sizeof(double));

    for (i = 0; i < LAYERS; i++) func[i] = 1 - R[i] * U[i];
    trapz(LAYERS, eta, func, &params.delta_1_ast[face]);
    params.delta_1_ast[face] = delta * params.delta_1_ast[face];

    for (i = 0; i < LAYERS; i++) func[i] = - R[i] * W[i];
    trapz(LAYERS, eta, func, &params.delta_2_ast[face]);
    params.delta_2_ast[face] = delta * params.delta_2_ast[face];

    for (i = 0; i < LAYERS; i++) func[i] = 1 - R[i] * U[i] * U[i];
    trapz(LAYERS, eta, func, &params.phi_11[face]);
    params.phi_11[face] = delta * params.phi_11[face];

    for (i = 0; i < LAYERS; i++) func[i] = - R[i] * U[i] * W[i];
    trapz(LAYERS, eta, func, &params.phi_12[face]);
    params.phi_12[face] = delta * params.phi_12[face];
    params.phi_21[face] = params.phi_12[face];

    for (i = 0; i < LAYERS; i++) func[i] = - R[i] * W[i] * W[i];
    trapz(LAYERS, eta, func, &params.phi_22[face]);
    params.phi_22[face] = delta * params.phi_22[face];

    for (i = 0; i < LAYERS; i++) func[i] = 1 - R[i] * U[i] * (U[i] * U[i] + W[i] * W[i]);
    trapz(LAYERS, eta, func, &params.phi_1_ast[face]);
    params.phi_1_ast[face] = delta * params.phi_1_ast[face];

    for (i = 0; i < LAYERS; i++) func[i] = - R[i] * W[i] * (U[i] * U[i] + W[i] * W[i]);
    trapz(LAYERS, eta, func, &params.phi_2_ast[face]);
    params.phi_2_ast[face] = delta * params.phi_2_ast[face];

    for (i = 0; i < LAYERS; i++) func[i] = 1 - U[i];
    trapz(LAYERS, eta, func, &params.delta_1_line[face]);
    params.delta_1_line[face] = delta * params.delta_1_line[face];

    for (i = 0; i < LAYERS; i++) func[i] = - W[i];
    trapz(LAYERS, eta, func, &params.delta_2_line[face]);
    params.delta_2_line[face] = delta * params.delta_2_line[face];

    for (i = 0; i < LAYERS; i++) func[i] = - Psi * R[i] * U[i] * (U[i] * U[i] + W[i] * W[i]);
    trapz(LAYERS, eta, func, &params.theta_1_o[face]);
    params.theta_1_o[face] = delta * params.theta_1_o[face];

    for (i = 0; i < LAYERS; i++) func[i] = - Psi * R[i] * W[i] * (U[i] * U[i] + W[i] * W[i]);
    trapz(LAYERS, eta, func, &params.theta_2_o[face]);
    params.theta_2_o[face] = delta * params.theta_2_o[face];

    for (i = 0; i < LAYERS; i++) func[i] = - Psi * U[i];
    trapz(LAYERS, eta, func, &params.delta_1_o[face]);
    params.delta_1_o[face] = delta * params.delta_1_o[face];

    for (i = 0; i < LAYERS; i++) func[i] = - Psi * W[i];
    trapz(LAYERS, eta, func, &params.delta_2_o[face]);
    params.delta_2_o[face] = delta * params.delta_2_o[face];

    for (i = 0; i < LAYERS; i++) func[i] = S[i] * dUdeta[i] + T[i] * dWdeta[i];
    trapz(LAYERS, eta, func, &params.C_D[face]);
    params.C_D[face] = delta * params.C_D[face];

    for (i = 0; i < LAYERS; i++) func[i] = S[i] * dWdeta[i] - T[i] * dUdeta[i];
    trapz(LAYERS, eta, func, &params.C_D_x[face]);
    params.C_D_x[face] = delta * params.C_D_x[face];

    for (i = 0; i < LAYERS; i++) func[i] = Psi * (S[i] * dUdeta[i] + T[i] * dWdeta[i]); //func[i] = S[i] * (dpsideta[i] * U[i] + psi[i] * dUdeta[i]) + T[i] * (dpsideta[i] * W[i] + psi[i] * dWdeta[i]);
    trapz(LAYERS, eta, func, &params.C_D_o[face]);
    params.C_D_o[face] = delta * params.C_D_o[face];

    params.S0[face] = S[0];
    params.T0[face] = T[0];

    for (i = 0; i < LAYERS; i++) func[i] = S[i] * U[i];
    trapz(LAYERS, eta, func, &params.Ktauxx[face]);
    params.Ktauxx[face] = delta * params.Ktauxx[face];

    for (i = 0; i < LAYERS; i++) func[i] = S[i] * W[i];
    trapz(LAYERS, eta, func, &params.Ktauxy[face]);
    params.Ktauxy[face] = delta * params.Ktauxy[face];

    for (i = 0; i < LAYERS; i++) func[i] = T[i] * U[i];
    trapz(LAYERS, eta, func, &params.Ktauyx[face]);
    params.Ktauyx[face] = delta * params.Ktauyx[face];

    for (i = 0; i < LAYERS; i++) func[i] = T[i] * W[i];
    trapz(LAYERS, eta, func, &params.Ktauyy[face]);
    params.Ktauyy[face] = delta * params.Ktauyy[face];

    double Re_theta;
    double func1;
    double func2;
    double H_k;
    double fn;
    
    H_k = (absValue(division(params.delta_1_ast[face], params.phi_11[face] - params.delta_1_ast[face])) - 0.29 * mach * mach) / (1 + 0.113 * mach * mach);

    if (H_k < 1.8) {
        
        params.Sx[face] = 0.0;

    } else {

        Re_theta = velocity * density * absValue(params.phi_11[face] - params.delta_1_ast[face]) / viscosity;
        func1 = 0.01 * sqrt(pow(2.4 * H_k - 3.7 + 2.5 * tanh(1.5 * H_k - 4.65), 2) + 0.25);
        func2 = pow(10, (1.415 / (H_k - 1) - 0.489) * tanh(20 / (H_k - 1) - 12.9) + 3.295 / (H_k - 1) + 0.44);
        fn = func1 * (Re_theta - func2);

        if (fn < 0) {
            params.Sx[face] = 0.0;
        } else {
            params.Sx[face] = amp * fn / ((params.phi_11[face] - params.delta_1_ast[face]) * density * velocity * velocity);
        }

    }

    H_k = (absValue(division(params.delta_2_ast[face], params.phi_11[face] - params.delta_1_ast[face])) - 0.29 * mach * mach) / (1 + 0.113 * mach * mach);

    if (H_k < 1.8) {
        
        params.Sy[face] = 0.0;

    } else {

        Re_theta = velocity * density * absValue(params.phi_22[face] - params.delta_2_ast[face]) / viscosity;
        func1 = 0.01 * sqrt(pow(2.4 * H_k - 3.7 + 2.5 * tanh(1.5 * H_k - 4.65), 2) + 0.25);
        func2 = pow(10, (1.415 / (H_k - 1) - 0.489) * tanh(20 / (H_k - 1) - 12.9) + 3.295 / (H_k - 1) + 0.44);
        fn = func1 * (Re_theta - func2);

        if (fn < 0) {
            params.Sy[face] = 0.0;
        } else {
            params.Sy[face] = amp * fn / ((params.phi_11[face] - params.delta_1_ast[face]) * density * velocity * velocity);
        }
        
    }

    // Free
    free(eta);
    free(U);
    free(W);
    free(dUdeta);
    free(dWdeta);
    free(S);
    free(T);
    free(R);
    free(psi);
    free(dpsideta);
    free(func);

}

void mallocIntegralThickness(int nf, struct IntegralThickness *params) {
    params->C_D = (double*)malloc(nf * sizeof(double));
    params->C_D_o = (double*)malloc(nf * sizeof(double));
    params->C_D_x = (double*)malloc(nf * sizeof(double));
    params->delta_1_ast = (double*)malloc(nf * sizeof(double));
    params->delta_1_line = (double*)malloc(nf * sizeof(double));
    params->delta_1_o = (double*)malloc(nf * sizeof(double));
    params->delta_2_ast = (double*)malloc(nf * sizeof(double));
    params->delta_2_line = (double*)malloc(nf * sizeof(double));
    params->delta_2_o = (double*)malloc(nf * sizeof(double));
    params->Ktauxx = (double*)malloc(nf * sizeof(double));
    params->Ktauxy = (double*)malloc(nf * sizeof(double));
    params->Ktauyx = (double*)malloc(nf * sizeof(double));
    params->Ktauyy = (double*)malloc(nf * sizeof(double));
    params->phi_11 = (double*)malloc(nf * sizeof(double));
    params->phi_12 = (double*)malloc(nf * sizeof(double));
    params->phi_1_ast = (double*)malloc(nf * sizeof(double));
    params->phi_21 = (double*)malloc(nf * sizeof(double));
    params->phi_22 = (double*)malloc(nf * sizeof(double));
    params->phi_2_ast = (double*)malloc(nf * sizeof(double));
    params->S0 = (double*)malloc(nf * sizeof(double));
    params->Sx = (double*)malloc(nf * sizeof(double));
    params->Sy = (double*)malloc(nf * sizeof(double));
    params->T0 = (double*)malloc(nf * sizeof(double));
    params->theta_1_o = (double*)malloc(nf * sizeof(double));
    params->theta_2_o = (double*)malloc(nf * sizeof(double));
}

void calculateParametersGradients(int nf, int *faces, double *vertices, double *faceCenters, double *velx, double *vely, double *velz, double *transpiration, double *e3, struct VerticeConnection *vertices_connection, struct FaceGradients *facesGradients) {
    
    /* Parameters */
    int face, i;
    struct Point vel, vel1, vel2, vel3;
    double scalar_vel;
    double velNorm, velNorm1, velNorm2, velNorm3;
    struct Point s1, s2, s3;
    double phi, phi1, phi2, phi3;
    double u, v;
    struct Point point1, point2, point3;
    struct Point n1, n2, n3;

    for (face = 0; face < nf; face++) {
    
        /* Face and edges velocities */
        vel.x = velx[face] - transpiration[face] * e3[3 * face];
        vel.y = vely[face] - transpiration[face] * e3[3 * face + 1];
        vel.z = velz[face] - transpiration[face] * e3[3 * face + 2];

        vel1.x = 0.0; vel1.y = 0.0; vel1.z = 0.0;
        vel2.x = 0.0; vel2.y = 0.0; vel2.z = 0.0;
        vel3.x = 0.0; vel3.y = 0.0; vel3.z = 0.0;

        // Vertice 1
        for(i = 0; i < vertices_connection[faces[3 * face]].n; i++) {
            vel1.x = vel1.x + velx[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i];
            vel1.y = vel1.y + vely[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i];
            vel1.z = vel1.z + velz[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i];
        }

        scalar_vel = vel1.x * e3[3 * face] + vel1.y * e3[3 * face + 1] + vel1.z * e3[3 * face + 2];

        vel1.x = vel1.x - scalar_vel * e3[3 * face];
        vel1.y = vel1.y - scalar_vel * e3[3 * face + 1];
        vel1.z = vel1.z - scalar_vel * e3[3 * face + 2];

        // Vertice 2
        for(i = 0; i < vertices_connection[faces[3 * face + 1]].n; i++) {
            vel2.x = vel2.x + velx[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i];
            vel2.y = vel2.y + vely[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i];
            vel2.z = vel2.z + velz[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i];
        }

        scalar_vel = vel2.x * e3[3 * face] + vel2.y * e3[3 * face + 1] + vel2.z * e3[3 * face + 2];

        vel2.x = vel2.x - scalar_vel * e3[3 * face];
        vel2.y = vel2.y - scalar_vel * e3[3 * face + 1];
        vel2.z = vel2.z - scalar_vel * e3[3 * face + 2];

        // Vertice 3
        for(i = 0; i < vertices_connection[faces[3 * face + 2]].n; i++) {
            vel3.x = vel3.x + velx[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i];
            vel3.y = vel3.y + vely[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i];
            vel3.z = vel3.z + velz[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i];
        }

        scalar_vel = vel3.x * e3[3 * face] + vel3.y * e3[3 * face + 1] + vel3.z * e3[3 * face + 2];

        vel3.x = vel3.x - scalar_vel * e3[3 * face];
        vel3.y = vel3.y - scalar_vel * e3[3 * face + 1];
        vel3.z = vel3.z - scalar_vel * e3[3 * face + 2];

        /* Vel norm */
        

        velNorm = (vel.x * vel.x + vel.y * vel.y + vel.z * vel.z);
        velNorm1 = (vel1.x * vel1.x + vel1.y * vel1.y + vel1.z * vel1.z) / velNorm;
        velNorm2 = (vel2.x * vel2.x + vel2.y * vel2.y + vel2.z * vel2.z) / velNorm;
        velNorm3 = (vel3.x * vel3.x + vel3.y * vel3.y + vel3.z * vel3.z) / velNorm;

        /* Angles */
        

        /* Face and edges angles */
        scalar_vel = sqrt(vel.x * vel.x + vel.y * vel.y + vel.z * vel.z);
        s1.x = vel.x / scalar_vel; s1.y = vel.y / scalar_vel; s1.z = vel.z / scalar_vel;
        s3.x = e3[3 * face]; s3.y = e3[3 * face + 1]; s3.z = e3[3 * face + 2];
        s2 = cross(s3, s1);

        phi = 0.0;

        // Vertice 1
        phi1 = atan2(s2.x * vel1.x + s2.y * vel1.y + s2.z * vel1.z, s1.x * vel1.x + s1.y * vel1.y + s1.z * vel1.z);

        // Vertice 2
        phi2 = atan2(s2.x * vel2.x + s2.y * vel2.y + s2.z * vel2.z, s1.x * vel2.x + s1.y * vel2.y + s1.z * vel2.z);

        // Vertice 3
        phi3 = atan2(s2.x * vel3.x + s2.y * vel3.y + s2.z * vel3.z, s1.x * vel3.x + s1.y * vel3.y + s1.z * vel3.z);

        /* Gradients */
        

        // grad_q2
        point1.x = (vertices[3 * faces[3 * face]] - faceCenters[3 * face]) * s1.x + (vertices[3 * faces[3 * face] + 1] - faceCenters[3 * face + 1]) * s1.y + (vertices[3 * faces[3 * face] + 2] - faceCenters[3 * face + 2]) * s1.z;
        point1.y = (vertices[3 * faces[3 * face]] - faceCenters[3 * face]) * s2.x + (vertices[3 * faces[3 * face] + 1] - faceCenters[3 * face + 1]) * s2.y + (vertices[3 * faces[3 * face] + 2] - faceCenters[3 * face + 2]) * s2.z;
        point1.z = velNorm1 - 1.0;

        point2.x = (vertices[3 * faces[3 * face + 1]] - faceCenters[3 * face]) * s1.x + (vertices[3 * faces[3 * face + 1] + 1] - faceCenters[3 * face + 1]) * s1.y + (vertices[3 * faces[3 * face + 1] + 2] - faceCenters[3 * face + 2]) * s1.z;
        point2.y = (vertices[3 * faces[3 * face + 1]] - faceCenters[3 * face]) * s2.x + (vertices[3 * faces[3 * face + 1] + 1] - faceCenters[3 * face + 1]) * s2.y + (vertices[3 * faces[3 * face + 1] + 2] - faceCenters[3 * face + 2]) * s2.z;
        point2.z = velNorm2 - 1.0;

        point3.x = (vertices[3 * faces[3 * face + 2]] - faceCenters[3 * face]) * s1.x + (vertices[3 * faces[3 * face + 2] + 1] - faceCenters[3 * face + 1]) * s1.y + (vertices[3 * faces[3 * face + 2] + 2] - faceCenters[3 * face + 2]) * s1.z;
        point3.y = (vertices[3 * faces[3 * face + 2]] - faceCenters[3 * face]) * s2.x + (vertices[3 * faces[3 * face + 2] + 1] - faceCenters[3 * face + 1]) * s2.y + (vertices[3 * faces[3 * face + 2] + 2] - faceCenters[3 * face + 2]) * s2.z;
        point3.z = velNorm3 - 1.0;

        n1 = cross(point1, point2);
        n2 = cross(point2, point3);
        n3 = cross(point3, point1);

        facesGradients[face].grad_u2_x = - (n1.x / n1.z + n2.x / n2.z + n3.x / n3.z);
        facesGradients[face].grad_u2_y = - (n1.y / n1.z + n2.y / n2.z + n3.y / n3.z);

        // grad_q2
        point1.z = phi1;
        point2.z = phi2;
        point3.z = phi3;

        n1 = cross(point1, point2);
        n2 = cross(point2, point3);
        n3 = cross(point3, point1);

        facesGradients[face].grad_phi_x = - (n1.x / n1.z + n2.x / n2.z + n3.x / n3.z);
        facesGradients[face].grad_phi_y = - (n1.y / n1.z + n2.y / n2.z + n3.y / n3.z);

    }

}

void calculateParametersDivergents(int face, double *velNorm, int *faces, double *vertices, double *velx, double *vely, double *velz, double *transpiration, double *e3, struct VerticeConnection *vertices_connection, struct IntegralThickness params, struct FaceDivergents *faceDivergents) {

    /* Face base vectors */
    struct Point s12, s23, s31, l1, l2, l3, l1_aux, l2_aux, l3_aux, s1, s2, s3;
    double s12_norm, s23_norm, s31_norm, vec_norm;

    // Edges
    s12.x = vertices[3 * faces[3 * face + 1]] - vertices[3 * faces[3 * face]];
    s12.y = vertices[3 * faces[3 * face + 1] + 1] - vertices[3 * faces[3 * face] + 1];
    s12.z = vertices[3 * faces[3 * face + 1] + 2] - vertices[3 * faces[3 * face] + 2];

    s12_norm = norm(s12); s12.x = s12.x / s12_norm; s12.y = s12.y / s12_norm; s12.z = s12.z / s12_norm;

    s23.x = vertices[3 * faces[3 * face + 2]] - vertices[3 * faces[3 * face + 1]];
    s23.y = vertices[3 * faces[3 * face + 2] + 1] - vertices[3 * faces[3 * face + 1] + 1];
    s23.z = vertices[3 * faces[3 * face + 2] + 2] - vertices[3 * faces[3 * face + 1] + 2];

    s23_norm = norm(s23); s23.x = s23.x / s23_norm; s23.y = s23.y / s23_norm; s23.z = s23.z / s23_norm;

    s31.x = vertices[3 * faces[3 * face]] - vertices[3 * faces[3 * face + 2]];
    s31.y = vertices[3 * faces[3 * face] + 1] - vertices[3 * faces[3 * face + 2] + 1];
    s31.z = vertices[3 * faces[3 * face] + 2] - vertices[3 * faces[3 * face + 2] + 2];

    s31_norm = norm(s31); s31.x = s31.x / s31_norm; s31.y = s31.y / s31_norm; s31.z = s31.z / s31_norm;

    // Velocity base
    s3.x = e3[3 * face]; s3.y = e3[3 * face + 1]; s3.z = e3[3 * face + 2];
    s1.x = velx[face] - transpiration[face] * s3.x; s1.y = vely[face] - transpiration[face] * s3.y; s1.z = velz[face] - transpiration[face] * s3.z; vec_norm = norm(s1); s1.x = s1.x / vec_norm; s1.y = s1.y / vec_norm; s1.z = s1.z / vec_norm;
    s2 = cross(s3, s1);

    // Lateral vectors
    l1_aux = cross(s12, s3);
    l2_aux = cross(s23, s3);
    l3_aux = cross(s31, s3);

    l1.x = dot(s1, l1_aux); l1.y = dot(s2, l1_aux); l1.z = 0.0;
    l2.x = dot(s1, l2_aux); l2.y = dot(s2, l2_aux); l2.z = 0.0;
    l3.x = dot(s1, l3_aux); l3.y = dot(s2, l3_aux); l3.z = 0.0;

    // Vertices values
    int i;
    double a, b;
    double delta_1_ast_1 = 0.0; double delta_2_ast_1 = 0.0; double delta_1_ast_2 = 0.0; double delta_2_ast_2 = 0.0; double delta_1_ast_3 = 0.0; double delta_2_ast_3 = 0.0;
    double phi_11_1 = 0.0; double phi_12_1 = 0.0; double phi_11_2 = 0.0; double phi_12_2 = 0.0; double phi_11_3 = 0.0; double phi_12_3 = 0.0;
    double phi_21_1 = 0.0; double phi_22_1 = 0.0; double phi_21_2 = 0.0; double phi_22_2 = 0.0; double phi_21_3 = 0.0; double phi_22_3 = 0.0;
    double phi_1_ast_1 = 0.0; double phi_2_ast_1 = 0.0; double phi_1_ast_2 = 0.0; double phi_2_ast_2 = 0.0; double phi_1_ast_3 = 0.0; double phi_2_ast_3 = 0.0;
    double theta_1_o_1 = 0.0; double theta_2_o_1 = 0.0; double theta_1_o_2 = 0.0; double theta_2_o_2 = 0.0; double theta_1_o_3 = 0.0; double theta_2_o_3 = 0.0;
    double ktau_xx_1 = 0.0; double ktau_xy_1 = 0.0; double ktau_xx_2 = 0.0; double ktau_xy_2 = 0.0; double ktau_xx_3 = 0.0; double ktau_xy_3 = 0.0;
    double ktau_yx_1 = 0.0; double ktau_yy_1 = 0.0; double ktau_yx_2 = 0.0; double ktau_yy_2 = 0.0; double ktau_yx_3 = 0.0; double ktau_yy_3 = 0.0;
    
    double factor_1;
    double factor_2;
    double factor_3;

    for(i = 0; i < vertices_connection[faces[3 * face]].n; i++) {

        factor_1 = division(velNorm[vertices_connection[faces[3 * face]].faces[i]], velNorm[face]);
        factor_2 = pow(division(velNorm[vertices_connection[faces[3 * face]].faces[i]], velNorm[face]), 2);
        factor_3 = pow(division(velNorm[vertices_connection[faces[3 * face]].faces[i]], velNorm[face]), 3);

        // delta_1_ast / delta_2_ast
        delta_1_ast_1 = delta_1_ast_1 + params.delta_1_ast[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_1;
        delta_2_ast_1 = delta_2_ast_1 + params.delta_2_ast[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_1;

        // phi_11 / phi_12
        phi_11_1 = phi_11_1 + params.phi_11[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_2;
        phi_12_1 = phi_12_1 + params.phi_12[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_2;

        // phi_21 / phi_22
        phi_21_1 = phi_21_1 + params.phi_21[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_2;
        phi_22_1 = phi_22_1 + params.phi_22[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_2;

        // phi_1_ast / phi_2_ast
        phi_1_ast_1 = phi_1_ast_1 + params.phi_1_ast[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_3;
        phi_2_ast_1 = phi_2_ast_1 + params.phi_2_ast[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_3;

        // theta_1_o / theta_2_o
        theta_1_o_1 = theta_1_o_1 + params.theta_1_o[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_3;
        theta_2_o_1 = theta_2_o_1 + params.theta_2_o[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_3;

        // ktau_xx / ktau_xy
        ktau_xx_1 = ktau_xx_1 + params.Ktauxx[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_3;
        ktau_xy_1 = ktau_xy_1 + params.Ktauxy[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_3;

        ktau_yx_1 = ktau_yx_1 + params.Ktauyx[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_3;
        ktau_yy_1 = ktau_yy_1 + params.Ktauyy[vertices_connection[faces[3 * face]].faces[i]] * vertices_connection[faces[3 * face]].coeffs[i] * factor_3;

    }

    for(i = 0; i < vertices_connection[faces[3 * face + 1]].n; i++) {
        
        factor_1 = division(velNorm[vertices_connection[faces[3 * face + 1]].faces[i]], velNorm[face]);
        factor_2 = pow(division(velNorm[vertices_connection[faces[3 * face + 1]].faces[i]], velNorm[face]), 2);
        factor_3 = pow(division(velNorm[vertices_connection[faces[3 * face + 1]].faces[i]], velNorm[face]), 3);

        // delta_1_ast / delta_2_ast
        delta_1_ast_2 = delta_1_ast_2 + params.delta_1_ast[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_1;
        delta_2_ast_2 = delta_2_ast_2 + params.delta_2_ast[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_1;

        // phi_11 / phi_12
        phi_11_2 = phi_11_2 + params.phi_11[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_2;
        phi_12_2 = phi_12_2 + params.phi_12[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_2;

        // phi_21 / phi_22
        phi_21_2 = phi_21_2 + params.phi_21[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_2;
        phi_22_2 = phi_22_2 + params.phi_22[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_2;

        // phi_1_ast / phi_2_ast
        phi_1_ast_2 = phi_1_ast_2 + params.phi_1_ast[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_3;
        phi_2_ast_2 = phi_2_ast_2 + params.phi_2_ast[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_3;

        // theta_1_o / theta_2_o
        theta_1_o_2 = theta_1_o_2 + params.theta_1_o[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_3;
        theta_2_o_2 = theta_2_o_2 + params.theta_2_o[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_3;

        // ktau_xx / ktau_xy
        ktau_xx_2 = ktau_xx_2 + params.Ktauxx[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_3;
        ktau_xy_2 = ktau_xy_2 + params.Ktauxy[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_3;

        ktau_yx_2 = ktau_yx_2 + params.Ktauyx[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_3;
        ktau_yy_2 = ktau_yy_2 + params.Ktauyy[vertices_connection[faces[3 * face + 1]].faces[i]] * vertices_connection[faces[3 * face + 1]].coeffs[i] * factor_3;

    }

    for(i = 0; i < vertices_connection[faces[3 * face + 2]].n; i++) {

        factor_1 = division(velNorm[vertices_connection[faces[3 * face + 2]].faces[i]], velNorm[face]);
        factor_2 = pow(division(velNorm[vertices_connection[faces[3 * face + 2]].faces[i]], velNorm[face]), 2);
        factor_3 = pow(division(velNorm[vertices_connection[faces[3 * face + 2]].faces[i]], velNorm[face]), 3);

        // delta_1_ast / delta_2_ast
        delta_1_ast_3 = delta_1_ast_3 + params.delta_1_ast[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_1;
        delta_2_ast_3 = delta_2_ast_3 + params.delta_2_ast[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_1;

        // phi_11 / phi_12
        phi_11_3 = phi_11_3 + params.phi_11[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_2;
        phi_12_3 = phi_12_3 + params.phi_12[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_2;

        // phi_21 / phi_22
        phi_21_3 = phi_21_3 + params.phi_21[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_2;
        phi_22_3 = phi_22_3 + params.phi_22[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_2;

        // phi_1_ast / phi_2_ast
        phi_1_ast_3 = phi_1_ast_3 + params.phi_1_ast[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_3;
        phi_2_ast_3 = phi_2_ast_3 + params.phi_2_ast[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_3;

        // theta_1_o / theta_2_o
        theta_1_o_3 = theta_1_o_3 + params.theta_1_o[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_3;
        theta_2_o_3 = theta_2_o_3 + params.theta_2_o[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_3;

        // ktau_xx / ktau_xy
        ktau_xx_3 = ktau_xx_3 + params.Ktauxx[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_3;
        ktau_xy_3 = ktau_xy_3 + params.Ktauxy[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_3;

        ktau_yx_3 = ktau_yx_3 + params.Ktauyx[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_3;
        ktau_yy_3 = ktau_yy_3 + params.Ktauyy[vertices_connection[faces[3 * face + 2]].faces[i]] * vertices_connection[faces[3 * face + 2]].coeffs[i] * factor_3;

    }

    // M
    faceDivergents->div_delta_1_2_ast = s12_norm * 0.5 * ((delta_1_ast_1 * l1.x + delta_2_ast_1 * l1.y) + (delta_1_ast_2 * l1.x + delta_2_ast_2 * l1.y));
    faceDivergents->div_delta_1_2_ast = faceDivergents->div_delta_1_2_ast + s23_norm * 0.5 * ((delta_1_ast_2 * l2.x + delta_2_ast_2 * l2.y) + (delta_1_ast_3 * l2.x + delta_2_ast_3 * l2.y));
    faceDivergents->div_delta_1_2_ast = faceDivergents->div_delta_1_2_ast + s31_norm * 0.5 * ((delta_1_ast_3 * l3.x + delta_2_ast_3 * l3.y) + (delta_1_ast_1 * l3.x + delta_2_ast_1 * l3.y));

    // J
    faceDivergents->div_phi_11_12 = s12_norm * 0.5 * ((phi_11_1 * l1.x + phi_12_1 * l1.y) + (phi_11_2 * l1.x + phi_12_2 * l1.y));
    faceDivergents->div_phi_11_12 = faceDivergents->div_phi_11_12 + s23_norm * 0.5 * ((phi_11_2 * l2.x + phi_12_2 * l2.y) + (phi_11_3 * l2.x + phi_12_3 * l2.y));
    faceDivergents->div_phi_11_12 = faceDivergents->div_phi_11_12 + s31_norm * 0.5 * ((phi_11_3 * l3.x + phi_12_3 * l3.y) + (phi_11_1 * l3.x + phi_12_1 * l3.y));

    faceDivergents->div_phi_21_22 = s12_norm * 0.5 * ((phi_21_1 * l1.x + phi_22_1 * l1.y) + (phi_21_2 * l1.x + phi_22_2 * l1.y));
    faceDivergents->div_phi_21_22 = faceDivergents->div_phi_21_22 + s23_norm * 0.5 * ((phi_21_2 * l2.x + phi_22_2 * l2.y) + (phi_21_3 * l2.x + phi_22_3 * l2.y));
    faceDivergents->div_phi_21_22 = faceDivergents->div_phi_21_22 + s31_norm * 0.5 * ((phi_21_3 * l3.x + phi_22_3 * l3.y) + (phi_21_1 * l3.x + phi_22_1 * l3.y));

    // E
    faceDivergents->div_phi_1_2_ast = s12_norm * 0.5 * ((phi_1_ast_1 * l1.x + phi_2_ast_1 * l1.y) + (phi_1_ast_2 * l1.x + phi_2_ast_2 * l1.y));
    faceDivergents->div_phi_1_2_ast = faceDivergents->div_phi_1_2_ast + s23_norm * 0.5 * ((phi_1_ast_2 * l2.x + phi_2_ast_2 * l2.y) + (phi_1_ast_3 * l2.x + phi_2_ast_3 * l2.y));
    faceDivergents->div_phi_1_2_ast = faceDivergents->div_phi_1_2_ast + s31_norm * 0.5 * ((phi_1_ast_3 * l3.x + phi_2_ast_3 * l3.y) + (phi_1_ast_1 * l3.x + phi_2_ast_1 * l3.y));

    // Ko
    faceDivergents->div_theta_1_2_o = s12_norm * 0.5 * ((theta_1_o_1 * l1.x + theta_2_o_1 * l1.y) + (theta_1_o_2 * l1.x + theta_2_o_2 * l1.y));
    faceDivergents->div_theta_1_2_o = faceDivergents->div_theta_1_2_o + s23_norm * 0.5 * ((theta_1_o_2 * l2.x + theta_2_o_2 * l2.y) + (theta_1_o_3 * l2.x + theta_2_o_3 * l2.y));
    faceDivergents->div_theta_1_2_o = faceDivergents->div_theta_1_2_o + s31_norm * 0.5 * ((theta_1_o_3 * l3.x + theta_2_o_3 * l3.y) + (theta_1_o_1 * l3.x + theta_2_o_1 * l3.y));

    // Ktau
    faceDivergents->div_Ktau_xx_xy = s12_norm * 0.5 * ((ktau_xx_1 * l1.x + ktau_xy_1 * l1.y) + (ktau_xx_2 * l1.x + ktau_xy_2 * l1.y));
    faceDivergents->div_Ktau_xx_xy = faceDivergents->div_Ktau_xx_xy + s23_norm * 0.5 * ((ktau_xx_2 * l2.x + ktau_xy_2 * l2.y) + (ktau_xx_3 * l2.x + ktau_xy_3 * l2.y));
    faceDivergents->div_Ktau_xx_xy = faceDivergents->div_Ktau_xx_xy + s31_norm * 0.5 * ((ktau_xx_3 * l3.x + ktau_xy_3 * l3.y) + (ktau_xx_1 * l3.x + ktau_xy_1 * l3.y));

    faceDivergents->div_Ktau_yx_yy = s12_norm * 0.5 * ((ktau_yx_1 * l1.x + ktau_yy_1 * l1.y) + (ktau_yx_2 * l1.x + ktau_yy_2 * l1.y));
    faceDivergents->div_Ktau_yx_yy = faceDivergents->div_Ktau_yx_yy + s23_norm * 0.5 * ((ktau_yx_2 * l2.x + ktau_yy_2 * l2.y) + (ktau_yx_3 * l2.x + ktau_yy_3 * l2.y));
    faceDivergents->div_Ktau_yx_yy = faceDivergents->div_Ktau_yx_yy + s31_norm * 0.5 * ((ktau_yx_3 * l3.x + ktau_yy_3 * l3.y) + (ktau_yx_1 * l3.x + ktau_yy_1 * l3.y));
}

void getIntegralThicknessParams(int face, struct IntegralThickness params, double *out) {

    out[0] = params.C_D[face];
    out[1] = params.C_D_o[face];
    out[2] = params.C_D_x[face];
    out[3] = params.delta_1_ast[face];
    out[4] = params.delta_1_line[face];
    out[5] = params.delta_1_o[face];
    out[6] = params.delta_2_ast[face];
    out[7] = params.delta_2_line[face];
    out[8] = params.delta_2_o[face];
    out[9] = params.Ktauxx[face];
    out[10] = params.Ktauxy[face];
    out[11] = params.Ktauyx[face];
    out[12] = params.Ktauyy[face];
    out[13] = params.phi_11[face];
    out[14] = params.phi_12[face];
    out[15] = params.phi_1_ast[face];
    out[16] = params.phi_21[face];
    out[17] = params.phi_22[face];
    out[18] = params.phi_2_ast[face];
    out[19] = params.S0[face];
    out[20] = params.Sx[face];
    out[21] = params.Sy[face];
    out[22] = params.T0[face];
    out[23] = params.theta_1_o[face];
    out[24] = params.theta_2_o[face];

}

void setIntegralThicknessParams(int face, struct IntegralThickness params, double *out) {

    params.C_D[face] = out[0];
    params.C_D_o[face] = out[1];
    params.C_D_x[face] = out[2];
    params.delta_1_ast[face] = out[3];
    params.delta_1_line[face] = out[4];
    params.delta_1_o[face] = out[5];
    params.delta_2_ast[face] = out[6];
    params.delta_2_line[face] = out[7];
    params.delta_2_o[face] = out[8];
    params.Ktauxx[face] = out[9];
    params.Ktauxy[face] = out[10];
    params.Ktauyx[face] = out[11];
    params.Ktauyy[face] = out[12];
    params.phi_11[face] = out[13];
    params.phi_12[face] = out[14];
    params.phi_1_ast[face] = out[15];
    params.phi_21[face] = out[16];
    params.phi_22[face] = out[17];
    params.phi_2_ast[face] = out[18];
    params.S0[face] = out[19];
    params.Sx[face] = out[20];
    params.Sy[face] = out[21];
    params.T0[face] = out[22];
    params.theta_1_o[face] = out[23];
    params.theta_2_o[face] = out[24];

}

void calculateGradient(int face, int inte, double *tau_w_x, double *tau_w_y, double density, struct FaceEquations *equations, double eps, double *velNorm, int *faces, double *vertices, double *velx, double *vely, double *velz, double *transpiration, double *e3, struct VerticeConnection *vertices_connection, struct FacesConnection *faces_connection, struct FaceGradients *facesGradients, double *gradient, struct IntegralThickness params, struct IntegralThickness params_delta, struct IntegralThickness params_A, struct IntegralThickness params_B, struct IntegralThickness params_Psi, struct IntegralThickness params_amp) {

    /* Parameters */
    double transpiration_aux;
    int i;
    struct FaceDivergents faceDivergents;
    double *params_1 = (double*)malloc(25 * sizeof(double));
    double *params_2 = (double*)malloc(25 * sizeof(double));
    double *obj_func_eps = (double*)malloc(5 * (1 + faces_connection[face].n) * sizeof(double));
    double momentum_x, momentum_y, kinetic_energy, lateral_curvature, shear_stress_x, shear_stress_y;

    double f1 = 1.0;
    double f2 = 1.0;
    double f3 = 5.0;
    double f4 = 2.0;
    double f5 = 1.0;
    double f6 = 1.0;

    /* Reference values */
    calculateParametersDivergents(face, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, params, &faceDivergents);

    transpiration_aux = density * velNorm[face] * faceDivergents.div_delta_1_2_ast;
    tau_w_x[face] = density * velNorm[face] * velNorm[face] * params.S0[face];
    tau_w_y[face] = density * velNorm[face] * velNorm[face] * params.T0[face];

    equations->momentum_x = faceDivergents.div_phi_11_12 - faceDivergents.div_delta_1_2_ast - params.S0[face];
    equations->momentum_y = faceDivergents.div_phi_21_22 - params.T0[face];
    equations->kinetic_energy = faceDivergents.div_phi_1_2_ast - faceDivergents.div_delta_1_2_ast - (params.delta_1_line[face] * facesGradients[face].grad_u2_x + params.delta_2_line[face] * facesGradients[face].grad_u2_y) - 2 * params.C_D[face];
    equations->lateral_curvature = faceDivergents.div_theta_1_2_o + (params.phi_1_ast[face] * facesGradients[face].grad_phi_x + params.phi_2_ast[face] * facesGradients[face].grad_phi_y) + 0.5 * (params.delta_1_line[face] * facesGradients[face].grad_u2_y - params.delta_2_line[face] * facesGradients[face].grad_u2_x) - (params.delta_1_o[face] * facesGradients[face].grad_u2_x + params.delta_2_o[face] * facesGradients[face].grad_u2_y) + params.C_D_x[face] - 2 * params.C_D_o[face];
    equations->shear_stress_x = faceDivergents.div_Ktau_xx_xy - params.Sx[face];
    equations->shear_stress_y = faceDivergents.div_Ktau_yx_yy - params.Sy[face];
    equations->obj = f1 * pow(equations->momentum_x, 2) + f2 * pow(equations->momentum_y, 2) + f3 * pow(equations->kinetic_energy, 2) + f4 * pow(equations->lateral_curvature, 2) + f5 * pow(equations->shear_stress_x, 2) + f6 * pow(equations->shear_stress_y, 2);

    /* Current face derivatives */

    // delta
    getIntegralThicknessParams(face, params, params_1);
    getIntegralThicknessParams(face, params_delta, params_2);
    setIntegralThicknessParams(face, params, params_2);

    calculateParametersDivergents(face, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, params, &faceDivergents);

    momentum_x = faceDivergents.div_phi_11_12 - faceDivergents.div_delta_1_2_ast - params.S0[face];
    momentum_y = faceDivergents.div_phi_21_22 - params.T0[face];
    kinetic_energy = faceDivergents.div_phi_1_2_ast - faceDivergents.div_delta_1_2_ast - (params.delta_1_line[face] * facesGradients[face].grad_u2_x + params.delta_2_line[face] * facesGradients[face].grad_u2_y) - 2 * params.C_D[face];
    lateral_curvature = faceDivergents.div_theta_1_2_o + (params.phi_1_ast[face] * facesGradients[face].grad_phi_x + params.phi_2_ast[face] * facesGradients[face].grad_phi_y) + 0.5 * (params.delta_1_line[face] * facesGradients[face].grad_u2_y - params.delta_2_line[face] * facesGradients[face].grad_u2_x) - (params.delta_1_o[face] * facesGradients[face].grad_u2_x + params.delta_2_o[face] * facesGradients[face].grad_u2_y) + params.C_D_x[face] - 2 * params.C_D_o[face];
    shear_stress_x = faceDivergents.div_Ktau_xx_xy - params.Sx[face];
    shear_stress_y = faceDivergents.div_Ktau_yx_yy - params.Sy[face];
    obj_func_eps[0] = f1 * pow(momentum_x, 2) + f2 * pow(momentum_y, 2) + f3 * pow(kinetic_energy, 2) + f4 * pow(lateral_curvature, 2) + f5 * pow(shear_stress_x, 2) + f6 * pow(shear_stress_y, 2);

    // A
    getIntegralThicknessParams(face, params_A, params_2);
    setIntegralThicknessParams(face, params, params_2);

    calculateParametersDivergents(face, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, params, &faceDivergents);

    momentum_x = faceDivergents.div_phi_11_12 - faceDivergents.div_delta_1_2_ast - params.S0[face];
    momentum_y = faceDivergents.div_phi_21_22 - params.T0[face];
    kinetic_energy = faceDivergents.div_phi_1_2_ast - faceDivergents.div_delta_1_2_ast - (params.delta_1_line[face] * facesGradients[face].grad_u2_x + params.delta_2_line[face] * facesGradients[face].grad_u2_y) - 2 * params.C_D[face];
    lateral_curvature = faceDivergents.div_theta_1_2_o + (params.phi_1_ast[face] * facesGradients[face].grad_phi_x + params.phi_2_ast[face] * facesGradients[face].grad_phi_y) + 0.5 * (params.delta_1_line[face] * facesGradients[face].grad_u2_y - params.delta_2_line[face] * facesGradients[face].grad_u2_x) - (params.delta_1_o[face] * facesGradients[face].grad_u2_x + params.delta_2_o[face] * facesGradients[face].grad_u2_y) + params.C_D_x[face] - 2 * params.C_D_o[face];
    shear_stress_x = faceDivergents.div_Ktau_xx_xy - params.Sx[face];
    shear_stress_y = faceDivergents.div_Ktau_yx_yy - params.Sy[face];
    obj_func_eps[1] = f1 * pow(momentum_x, 2) + f2 * pow(momentum_y, 2) + f3 * pow(kinetic_energy, 2) + f4 * pow(lateral_curvature, 2) + f5 * pow(shear_stress_x, 2) + f6 * pow(shear_stress_y, 2);

    // B
    getIntegralThicknessParams(face, params_B, params_2);
    setIntegralThicknessParams(face, params, params_2);

    calculateParametersDivergents(face, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, params, &faceDivergents);

    momentum_x = faceDivergents.div_phi_11_12 - faceDivergents.div_delta_1_2_ast - params.S0[face];
    momentum_y = faceDivergents.div_phi_21_22 - params.T0[face];
    kinetic_energy = faceDivergents.div_phi_1_2_ast - faceDivergents.div_delta_1_2_ast - (params.delta_1_line[face] * facesGradients[face].grad_u2_x + params.delta_2_line[face] * facesGradients[face].grad_u2_y) - 2 * params.C_D[face];
    lateral_curvature = faceDivergents.div_theta_1_2_o + (params.phi_1_ast[face] * facesGradients[face].grad_phi_x + params.phi_2_ast[face] * facesGradients[face].grad_phi_y) + 0.5 * (params.delta_1_line[face] * facesGradients[face].grad_u2_y - params.delta_2_line[face] * facesGradients[face].grad_u2_x) - (params.delta_1_o[face] * facesGradients[face].grad_u2_x + params.delta_2_o[face] * facesGradients[face].grad_u2_y) + params.C_D_x[face] - 2 * params.C_D_o[face];
    shear_stress_x = faceDivergents.div_Ktau_xx_xy - params.Sx[face];
    shear_stress_y = faceDivergents.div_Ktau_yx_yy - params.Sy[face];
    obj_func_eps[2] = f1 * pow(momentum_x, 2) + f2 * pow(momentum_y, 2) + f3 * pow(kinetic_energy, 2) + f4 * pow(lateral_curvature, 2) + f5 * pow(shear_stress_x, 2) + f6 * pow(shear_stress_y, 2);

    // Psi
    getIntegralThicknessParams(face, params_Psi, params_2);
    setIntegralThicknessParams(face, params, params_2);

    calculateParametersDivergents(face, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, params, &faceDivergents);

    momentum_x = faceDivergents.div_phi_11_12 - faceDivergents.div_delta_1_2_ast - params.S0[face];
    momentum_y = faceDivergents.div_phi_21_22 - params.T0[face];
    kinetic_energy = faceDivergents.div_phi_1_2_ast - faceDivergents.div_delta_1_2_ast - (params.delta_1_line[face] * facesGradients[face].grad_u2_x + params.delta_2_line[face] * facesGradients[face].grad_u2_y) - 2 * params.C_D[face];
    lateral_curvature = faceDivergents.div_theta_1_2_o + (params.phi_1_ast[face] * facesGradients[face].grad_phi_x + params.phi_2_ast[face] * facesGradients[face].grad_phi_y) + 0.5 * (params.delta_1_line[face] * facesGradients[face].grad_u2_y - params.delta_2_line[face] * facesGradients[face].grad_u2_x) - (params.delta_1_o[face] * facesGradients[face].grad_u2_x + params.delta_2_o[face] * facesGradients[face].grad_u2_y) + params.C_D_x[face] - 2 * params.C_D_o[face];
    shear_stress_x = faceDivergents.div_Ktau_xx_xy - params.Sx[face];
    shear_stress_y = faceDivergents.div_Ktau_yx_yy - params.Sy[face];
    obj_func_eps[3] = f1 * pow(momentum_x, 2) + f2 * pow(momentum_y, 2) + f3 * pow(kinetic_energy, 2) + f4 * pow(lateral_curvature, 2) + f5 * pow(shear_stress_x, 2) + f6 * pow(shear_stress_y, 2);

    // amp
    getIntegralThicknessParams(face, params_amp, params_2);
    setIntegralThicknessParams(face, params, params_2);

    calculateParametersDivergents(face, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, params, &faceDivergents);

    momentum_x = faceDivergents.div_phi_11_12 - faceDivergents.div_delta_1_2_ast - params.S0[face];
    momentum_y = faceDivergents.div_phi_21_22 - params.T0[face];
    kinetic_energy = faceDivergents.div_phi_1_2_ast - faceDivergents.div_delta_1_2_ast - (params.delta_1_line[face] * facesGradients[face].grad_u2_x + params.delta_2_line[face] * facesGradients[face].grad_u2_y) - 2 * params.C_D[face];
    lateral_curvature = faceDivergents.div_theta_1_2_o + (params.phi_1_ast[face] * facesGradients[face].grad_phi_x + params.phi_2_ast[face] * facesGradients[face].grad_phi_y) + 0.5 * (params.delta_1_line[face] * facesGradients[face].grad_u2_y - params.delta_2_line[face] * facesGradients[face].grad_u2_x) - (params.delta_1_o[face] * facesGradients[face].grad_u2_x + params.delta_2_o[face] * facesGradients[face].grad_u2_y) + params.C_D_x[face] - 2 * params.C_D_o[face];
    shear_stress_x = faceDivergents.div_Ktau_xx_xy - params.Sx[face];
    shear_stress_y = faceDivergents.div_Ktau_yx_yy - params.Sy[face];
    obj_func_eps[4] = f1 * pow(momentum_x, 2) + f2 * pow(momentum_y, 2) + f3 * pow(kinetic_energy, 2) + f4 * pow(lateral_curvature, 2) + f5 * pow(shear_stress_x, 2) + f6 * pow(shear_stress_y, 2);

    setIntegralThicknessParams(face, params, params_1);

    for (i = 0; i < faces_connection[face].n; i++) {

        // delta
        getIntegralThicknessParams(faces_connection[face].faces[i], params, params_1);
        getIntegralThicknessParams(faces_connection[face].faces[i], params_delta, params_2);
        setIntegralThicknessParams(faces_connection[face].faces[i], params, params_2);

        calculateParametersDivergents(face, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, params, &faceDivergents);

        momentum_x = faceDivergents.div_phi_11_12 - faceDivergents.div_delta_1_2_ast - params.S0[face];
        momentum_y = faceDivergents.div_phi_21_22 - params.T0[face];
        kinetic_energy = faceDivergents.div_phi_1_2_ast - faceDivergents.div_delta_1_2_ast - (params.delta_1_line[face] * facesGradients[face].grad_u2_x + params.delta_2_line[face] * facesGradients[face].grad_u2_y) - 2 * params.C_D[face];
        lateral_curvature = faceDivergents.div_theta_1_2_o + (params.phi_1_ast[face] * facesGradients[face].grad_phi_x + params.phi_2_ast[face] * facesGradients[face].grad_phi_y) + 0.5 * (params.delta_1_line[face] * facesGradients[face].grad_u2_y - params.delta_2_line[face] * facesGradients[face].grad_u2_x) - (params.delta_1_o[face] * facesGradients[face].grad_u2_x + params.delta_2_o[face] * facesGradients[face].grad_u2_y) + params.C_D_x[face] - 2 * params.C_D_o[face];
        shear_stress_x = faceDivergents.div_Ktau_xx_xy - params.Sx[face];
        shear_stress_y = faceDivergents.div_Ktau_yx_yy - params.Sy[face];
        obj_func_eps[5 + 5 * i] = f1 * pow(momentum_x, 2) + f2 * pow(momentum_y, 2) + f3 * pow(kinetic_energy, 2) + f4 * pow(lateral_curvature, 2) + f5 * pow(shear_stress_x, 2) + f6 * pow(shear_stress_y, 2);

        // A
        getIntegralThicknessParams(faces_connection[face].faces[i], params_A, params_2);
        setIntegralThicknessParams(faces_connection[face].faces[i], params, params_2);

        calculateParametersDivergents(face, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, params, &faceDivergents);

        momentum_x = faceDivergents.div_phi_11_12 - faceDivergents.div_delta_1_2_ast - params.S0[face];
        momentum_y = faceDivergents.div_phi_21_22 - params.T0[face];
        kinetic_energy = faceDivergents.div_phi_1_2_ast - faceDivergents.div_delta_1_2_ast - (params.delta_1_line[face] * facesGradients[face].grad_u2_x + params.delta_2_line[face] * facesGradients[face].grad_u2_y) - 2 * params.C_D[face];
        lateral_curvature = faceDivergents.div_theta_1_2_o + (params.phi_1_ast[face] * facesGradients[face].grad_phi_x + params.phi_2_ast[face] * facesGradients[face].grad_phi_y) + 0.5 * (params.delta_1_line[face] * facesGradients[face].grad_u2_y - params.delta_2_line[face] * facesGradients[face].grad_u2_x) - (params.delta_1_o[face] * facesGradients[face].grad_u2_x + params.delta_2_o[face] * facesGradients[face].grad_u2_y) + params.C_D_x[face] - 2 * params.C_D_o[face];
        shear_stress_x = faceDivergents.div_Ktau_xx_xy - params.Sx[face];
        shear_stress_y = faceDivergents.div_Ktau_yx_yy - params.Sy[face];
        obj_func_eps[5 + 5 * i + 1] = f1 * pow(momentum_x, 2) + f2 * pow(momentum_y, 2) + f3 * pow(kinetic_energy, 2) + f4 * pow(lateral_curvature, 2) + f5 * pow(shear_stress_x, 2) + f6 * pow(shear_stress_y, 2);

        // B
        getIntegralThicknessParams(faces_connection[face].faces[i], params_B, params_2);
        setIntegralThicknessParams(faces_connection[face].faces[i], params, params_2);

        calculateParametersDivergents(face, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, params, &faceDivergents);

        momentum_x = faceDivergents.div_phi_11_12 - faceDivergents.div_delta_1_2_ast - params.S0[face];
        momentum_y = faceDivergents.div_phi_21_22 - params.T0[face];
        kinetic_energy = faceDivergents.div_phi_1_2_ast - faceDivergents.div_delta_1_2_ast - (params.delta_1_line[face] * facesGradients[face].grad_u2_x + params.delta_2_line[face] * facesGradients[face].grad_u2_y) - 2 * params.C_D[face];
        lateral_curvature = faceDivergents.div_theta_1_2_o + (params.phi_1_ast[face] * facesGradients[face].grad_phi_x + params.phi_2_ast[face] * facesGradients[face].grad_phi_y) + 0.5 * (params.delta_1_line[face] * facesGradients[face].grad_u2_y - params.delta_2_line[face] * facesGradients[face].grad_u2_x) - (params.delta_1_o[face] * facesGradients[face].grad_u2_x + params.delta_2_o[face] * facesGradients[face].grad_u2_y) + params.C_D_x[face] - 2 * params.C_D_o[face];
        shear_stress_x = faceDivergents.div_Ktau_xx_xy - params.Sx[face];
        shear_stress_y = faceDivergents.div_Ktau_yx_yy - params.Sy[face];
        obj_func_eps[5 + 5 * i + 2] = f1 * pow(momentum_x, 2) + f2 * pow(momentum_y, 2) + f3 * pow(kinetic_energy, 2) + f4 * pow(lateral_curvature, 2) + f5 * pow(shear_stress_x, 2) + f6 * pow(shear_stress_y, 2);

        // Psi
        getIntegralThicknessParams(faces_connection[face].faces[i], params_Psi, params_2);
        setIntegralThicknessParams(faces_connection[face].faces[i], params, params_2);

        calculateParametersDivergents(face, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, params, &faceDivergents);

        momentum_x = faceDivergents.div_phi_11_12 - faceDivergents.div_delta_1_2_ast - params.S0[face];
        momentum_y = faceDivergents.div_phi_21_22 - params.T0[face];
        kinetic_energy = faceDivergents.div_phi_1_2_ast - faceDivergents.div_delta_1_2_ast - (params.delta_1_line[face] * facesGradients[face].grad_u2_x + params.delta_2_line[face] * facesGradients[face].grad_u2_y) - 2 * params.C_D[face];
        lateral_curvature = faceDivergents.div_theta_1_2_o + (params.phi_1_ast[face] * facesGradients[face].grad_phi_x + params.phi_2_ast[face] * facesGradients[face].grad_phi_y) + 0.5 * (params.delta_1_line[face] * facesGradients[face].grad_u2_y - params.delta_2_line[face] * facesGradients[face].grad_u2_x) - (params.delta_1_o[face] * facesGradients[face].grad_u2_x + params.delta_2_o[face] * facesGradients[face].grad_u2_y) + params.C_D_x[face] - 2 * params.C_D_o[face];
        shear_stress_x = faceDivergents.div_Ktau_xx_xy - params.Sx[face];
        shear_stress_y = faceDivergents.div_Ktau_yx_yy - params.Sy[face];
        obj_func_eps[5 + 5 * i + 3] = f1 * pow(momentum_x, 2) + f2 * pow(momentum_y, 2) + f3 * pow(kinetic_energy, 2) + f4 * pow(lateral_curvature, 2) + f5 * pow(shear_stress_x, 2) + f6 * pow(shear_stress_y, 2);

        // amp
        getIntegralThicknessParams(faces_connection[face].faces[i], params_amp, params_2);
        setIntegralThicknessParams(faces_connection[face].faces[i], params, params_2);

        calculateParametersDivergents(face, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, params, &faceDivergents);

        momentum_x = faceDivergents.div_phi_11_12 - faceDivergents.div_delta_1_2_ast - params.S0[face];
        momentum_y = faceDivergents.div_phi_21_22 - params.T0[face];
        kinetic_energy = faceDivergents.div_phi_1_2_ast - faceDivergents.div_delta_1_2_ast - (params.delta_1_line[face] * facesGradients[face].grad_u2_x + params.delta_2_line[face] * facesGradients[face].grad_u2_y) - 2 * params.C_D[face];
        lateral_curvature = faceDivergents.div_theta_1_2_o + (params.phi_1_ast[face] * facesGradients[face].grad_phi_x + params.phi_2_ast[face] * facesGradients[face].grad_phi_y) + 0.5 * (params.delta_1_line[face] * facesGradients[face].grad_u2_y - params.delta_2_line[face] * facesGradients[face].grad_u2_x) - (params.delta_1_o[face] * facesGradients[face].grad_u2_x + params.delta_2_o[face] * facesGradients[face].grad_u2_y) + params.C_D_x[face] - 2 * params.C_D_o[face];
        shear_stress_x = faceDivergents.div_Ktau_xx_xy - params.Sx[face];
        shear_stress_y = faceDivergents.div_Ktau_yx_yy - params.Sy[face];
        obj_func_eps[5 + 5 * i + 4] = f1 * pow(momentum_x, 2) + f2 * pow(momentum_y, 2) + f3 * pow(kinetic_energy, 2) + f4 * pow(lateral_curvature, 2) + f5 * pow(shear_stress_x, 2) + f6 * pow(shear_stress_y, 2);

        setIntegralThicknessParams(faces_connection[face].faces[i], params, params_1);
        
    }

    // Calculate grandient
    gradient[0] = (obj_func_eps[0] - equations->obj) / eps;
    gradient[1] = (obj_func_eps[1] - equations->obj) / eps;
    gradient[2] = (obj_func_eps[2] - equations->obj) / eps;
    gradient[3] = (obj_func_eps[3] - equations->obj) / eps;
    gradient[4] = (obj_func_eps[4] - equations->obj) / eps;

    for (i = 0; i < faces_connection[face].n; i++) {
        gradient[5 + 5 * i] = (obj_func_eps[5 + 5 * i] - equations->obj) / eps;
        gradient[5 + 5 * i + 1] = (obj_func_eps[5 + 5 * i + 1] - equations->obj) / eps;
        gradient[5 + 5 * i + 2] = (obj_func_eps[5 + 5 * i + 2] - equations->obj) / eps;
        gradient[5 + 5 * i + 3] = (obj_func_eps[5 + 5 * i + 3] - equations->obj) / eps;
        gradient[5 + 5 * i + 4] = (obj_func_eps[5 + 5 * i + 4] - equations->obj) / eps;
    }

    transpiration[face] = absValue(transpiration_aux);

    /* Obj func to compare */
    equations->obj = absValue(equations->momentum_x) + absValue(equations->momentum_y) + absValue(equations->kinetic_energy) + absValue(equations->lateral_curvature) + absValue(equations->shear_stress_x) + absValue(equations->shear_stress_y);

}

double signDouble(double a) {
    if (a < 0.0) {
        return - 1.0;
    } else {
        return 1.0;
    }
}

void solveBoundaryLayer(int nf,
                        int nv,
                        struct Point freestream,
                        struct VerticeConnection *vertices_connection,
                        double *vertices,
                        int *faces,
                        double *facesCenter,
                        double *facesAreas,
                        double *e1, double *e2, double *e3,
                        double *p1, double *p2, double *p3,
                        double *transpiration,
                        double *delta, double *A, double *B, double *Psi, double *Ctau1, double *Ctau2,
                        double *tau_wall_x, double *tau_wall_y,
                        double *velNorm,
                        double *velx, double *vely, double *velz,
                        double *mach,
                        double density, double viscosity,
                        double *cp,
                        double *matrix, double *array,
                        double *matrixVelx, double *matrixVely, double *matrixVelz, double *arrayVel,
                        double *doublet,
                        double freestreamNorm) {
    
    /* Loop */
    int i, j, k;

    /* Equations parameters */
    struct IntegralThickness integralThickness;
    struct IntegralThickness integralThickness_delta;
    struct IntegralThickness integralThickness_A;
    struct IntegralThickness integralThickness_B;
    struct IntegralThickness integralThickness_Psi;
    struct IntegralThickness integralThickness_amp;

    mallocIntegralThickness(nf, &integralThickness);
    mallocIntegralThickness(nf, &integralThickness_delta);
    mallocIntegralThickness(nf, &integralThickness_A);
    mallocIntegralThickness(nf, &integralThickness_B);
    mallocIntegralThickness(nf, &integralThickness_Psi);
    mallocIntegralThickness(nf, &integralThickness_amp);

    struct FaceGradients *facesGradients = (struct FaceGradients*)malloc(nf * sizeof(struct FaceGradients));

    /* Parameters */
    double eps = 1e-8;
    double *delta_list, *A_list, *B_list, *Psi_list, *amp_list;
    double delta_norm, A_norm, B_norm, Psi_norm, amp_norm;

    delta_list = (double*)malloc(nf * sizeof(double));
    A_list = (double*)malloc(nf * sizeof(double));
    B_list = (double*)malloc(nf * sizeof(double));
    Psi_list = (double*)malloc(nf * sizeof(double));
    amp_list = (double*)malloc(nf * sizeof(double));

    delta_norm = 1e-3; A_norm = 1.0; B_norm = 0.01; Psi_norm = 0.01; amp_norm = 1e-6;

    double max_x_value = -10.0;
    for (i = 0; i < nv; i++) if (vertices[3 * i] > max_x_value) max_x_value = vertices[3 * i];

    double delta_factor;

    for (i = 0; i < nf; i++) {

        delta_factor = (freestream.x * e3[3 * i] + freestream.y * e3[3 * i + 1] + freestream.z * e3[3 * i + 2]) / freestreamNorm;

        if (delta_factor < 0) {
            delta_factor = pow(1.0 + delta_factor, 2);
        } else {
            delta_factor = 1.0;
        }

        delta_list[i] = (1 / delta_norm) * (delta_factor + 1e-8) * (5.0 * (max_x_value - facesCenter[3 * i]) / sqrt(density * freestreamNorm * (max_x_value - facesCenter[3 * i]) / viscosity));// (1 / delta_norm) * (0.000001 + 5.0 * (max_x_value - facesCenter[3 * i]) / sqrt(density * freestreamNorm * (max_x_value - facesCenter[3 * i]) / viscosity));
        A_list[i] = 1e-5;
        B_list[i] = 0.0;
        Psi_list[i] = 0.0;
        amp_list[i] = 1e-15;

    }

    /* Gradient */
    double *gradient;
    struct FaceEquations faceEquations;

    /* Faces connection */
    struct FacesConnection *faces_connection = (struct FacesConnection *)malloc(nf * sizeof(struct FacesConnection));
    calculateFacesConnection(nv, nf, faces, vertices_connection, faces_connection);

    /* Convergence */
    int max_interactions = 500;
    struct FaceEquations maxFaceEquations;
    struct FaceEquations previousMaxFaceEquations;
    double step_ref = 1e-1;
    double step;
    double step_aux;
    struct Point drag;
    double *max_inc = (double*)malloc(5 * sizeof(double)); max_inc[0] = 1e-4 * delta_norm; max_inc[1] = 1e-5 * A_norm; max_inc[2] = 1e-5 * B_norm; max_inc[3] = 1e-5 * Psi_norm; max_inc[4] = 1e-3 * amp_norm;
    double *steps = (double*)malloc(5 * sizeof(double));
    double *max_ratio = (double*)malloc(5 * sizeof(double));
    int check_max;
    double *delta_list_aux = (double*)malloc(nf * sizeof(double));
    double *A_list_aux = (double*)malloc(nf * sizeof(double));
    double *B_list_aux = (double*)malloc(nf * sizeof(double));
    double *Psi_list_aux = (double*)malloc(nf * sizeof(double));
    double *amp_list_aux = (double*)malloc(nf * sizeof(double));
    double alpha = 0.5;
    double grad_val;
    double **gradient_matrix = (double**)malloc(nf * sizeof(double*));
    
    for (i = 0; i < nf; i++) {
        gradient_matrix[i] = (double*)malloc(5 * (1 + faces_connection[i].n) * sizeof(double));
    }

    /* Surface base vectors */
    struct Point *s1 = (struct Point*)malloc(nf * sizeof(struct Point));
    struct Point *s2 = (struct Point*)malloc(nf * sizeof(struct Point));
    struct Point s3;
    double aux;

    for (i = 0; i < nf; i++) {

        aux = sqrt(velx[i] * velx[i] + vely[i] * vely[i] + velz[i] * velz[i]);
        
        s1[i].x = velx[i] / aux;
        s1[i].y = vely[i] / aux;
        s1[i].z = velz[i] / aux;

        s3.x = e3[3 * i];
        s3.y = e3[3 * i + 1];
        s3.z = e3[3 * i + 2];

        s2[i] = cross(s3, s1[i]);

    }

    printf("\n                                    Interactions\n");
    printf("--------------------------------------------------------------------------------------------\n");

    /* Loop */
    for (i = 1; i <= max_interactions; i++) {

        /* Calculate parameters gradients */
        calculateParametersGradients(nf, faces, vertices, facesCenter, velx, vely, velz, transpiration, e3, vertices_connection, facesGradients);

        /* Calculate parameters */
        for (j = 0; j < nf; j++) {
            calculateParameters(j, delta_norm * delta_list[j], A_norm * A_list[j], B_norm * B_list[j], Psi_norm * Psi_list[j], amp_norm * amp_list[j], velNorm[j], density, viscosity, mach[j], integralThickness);
            calculateParameters(j, delta_norm * (delta_list[j] + eps), A_norm * A_list[j], B_norm * B_list[j], Psi_norm * Psi_list[j], amp_norm * amp_list[j], velNorm[j], density, viscosity, mach[j], integralThickness_delta);
            calculateParameters(j, delta_norm * delta_list[j], A_norm * (A_list[j] + eps), B_norm * B_list[j], Psi_norm * Psi_list[j], amp_norm * amp_list[j], velNorm[j], density, viscosity, mach[j], integralThickness_A);
            calculateParameters(j, delta_norm * delta_list[j], A_norm * A_list[j], B_norm * (B_list[j] + eps), Psi_norm * Psi_list[j], amp_norm * amp_list[j], velNorm[j], density, viscosity, mach[j], integralThickness_B);
            calculateParameters(j, delta_norm * delta_list[j], A_norm * A_list[j], B_norm * B_list[j], Psi_norm * (Psi_list[j] + eps), amp_norm * amp_list[j], velNorm[j], density, viscosity, mach[j], integralThickness_Psi);
            calculateParameters(j, delta_norm * delta_list[j], A_norm * A_list[j], B_norm * B_list[j], Psi_norm * Psi_list[j], amp_norm * (amp_list[j] + eps), velNorm[j], density, viscosity, mach[j], integralThickness_amp);
        }

        /* Assign next values */
        for (j = 0; j < nf; j++) {
            delta_list_aux[j] = delta_list[j];
            A_list_aux[j] =  A_list[j];
            B_list_aux[j] = B_list[j];
            Psi_list_aux[j] = Psi_list[j];
            amp_list_aux[j] = amp_list[j];
        }

        /* Loop over faces */
        for (j = 0; j < nf; j++) {

            /* Calculate gradient */
            gradient = (double*)malloc(5 * (1 + faces_connection[j].n) * sizeof(double));

            calculateGradient(j, i, tau_wall_x, tau_wall_y, density, &faceEquations, eps, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, faces_connection, facesGradients, gradient, integralThickness, integralThickness_delta, integralThickness_A, integralThickness_B, integralThickness_Psi, integralThickness_amp);

            /* Step */
            step_aux = max_inc[0] / absValue(gradient[0]);
            if (max_inc[1] / absValue(gradient[1]) > step_aux) step_aux = max_inc[1] / absValue(gradient[1]);
            if (max_inc[2] / absValue(gradient[2]) > step_aux) step_aux = max_inc[2] / absValue(gradient[2]);
            if (max_inc[3] / absValue(gradient[3]) > step_aux) step_aux = max_inc[3] / absValue(gradient[3]);
            if (max_inc[4] / absValue(gradient[4]) > step_aux) step_aux = max_inc[4] / absValue(gradient[4]);

            for (k = 0; k < faces_connection[j].n; k++) {
                if (max_inc[0] / absValue(gradient[5 + 5 * k]) > step_aux) step_aux = max_inc[0] / absValue(gradient[5 + 5 * k]);
                if (max_inc[1] / absValue(gradient[5 + 5 * k + 1]) > step_aux) step_aux = max_inc[1] / absValue(gradient[5 + 5 * k + 1]);
                if (max_inc[2] / absValue(gradient[5 + 5 * k + 2]) > step_aux) step_aux = max_inc[2] / absValue(gradient[5 + 5 * k + 2]);
                if (max_inc[3] / absValue(gradient[5 + 5 * k + 3]) > step_aux) step_aux = max_inc[3] / absValue(gradient[5 + 5 * k + 3]);
                if (max_inc[4] / absValue(gradient[5 + 5 * k + 4]) > step_aux) step_aux = max_inc[4] / absValue(gradient[5 + 5 * k + 4]);
            }

            if (step_aux < step_ref) {
                step = step_aux;
                check_max = 0;
            } else {
                step = step_ref;
                check_max = 1;
            }

            /* Increase */
            if (i == 1) {

                delta_list_aux[j] = delta_list_aux[j] - step * gradient[0];
                A_list_aux[j] = A_list_aux[j] - step * gradient[1];
                B_list_aux[j] = B_list_aux[j] - step * gradient[2];
                Psi_list_aux[j] = Psi_list_aux[j] - step * gradient[3];
                amp_list_aux[j] = amp_list_aux[j] - step * gradient[4];

                if (delta_list_aux[j] < 1e-8) delta_list_aux[j] = 1e-10;
                if (amp_list_aux[j] < 0.0) amp_list_aux[j] = 0.0;

                for (k = 0; k < faces_connection[j].n; k++) {
                    delta_list_aux[faces_connection[j].faces[k]] = delta_list_aux[faces_connection[j].faces[k]] - step * gradient[5 + 5 * k];
                    A_list_aux[faces_connection[j].faces[k]] = A_list_aux[faces_connection[j].faces[k]] - step * gradient[5 + 5 * k + 1];
                    B_list_aux[faces_connection[j].faces[k]] = B_list_aux[faces_connection[j].faces[k]] - step * gradient[5 + 5 * k + 2];
                    Psi_list_aux[faces_connection[j].faces[k]] = Psi_list_aux[faces_connection[j].faces[k]] - step * gradient[5 + 5 * k + 3];
                    amp_list_aux[faces_connection[j].faces[k]] = amp_list_aux[faces_connection[j].faces[k]] - step * gradient[5 + 5 * k + 4];

                    if (delta_list_aux[faces_connection[j].faces[k]] < 1e-8) delta_list_aux[faces_connection[j].faces[k]] = 1e-10;
                    if (amp_list_aux[faces_connection[j].faces[k]] < 0.0) amp_list_aux[faces_connection[j].faces[k]] = 0.0;

                }

            } else {
                
                if (division(absValue(gradient[0]), absValue(gradient_matrix[j][0])) > 1.05) gradient[0] = 1.05 * absValue(gradient_matrix[j][0]) * signDouble(gradient[0]);
                delta_list_aux[j] = delta_list_aux[j] - step * (alpha * gradient[0] + (1 - alpha) * gradient_matrix[j][0]);
                
                if (division(absValue(gradient[1]), absValue(gradient_matrix[j][1])) > 1.05) gradient[1] = 1.05 * absValue(gradient_matrix[j][1]) * signDouble(gradient[1]);
                A_list_aux[j] = A_list_aux[j] - step * (alpha * gradient[1] + (1 - alpha) * gradient_matrix[j][1]);
                
                if (division(absValue(gradient[2]), absValue(gradient_matrix[j][2])) > 1.05) gradient[2] = 1.05 * absValue(gradient_matrix[j][2]) * signDouble(gradient[2]);
                B_list_aux[j] = B_list_aux[j] - step * (alpha * gradient[2] + (1 - alpha) * gradient_matrix[j][2]);
                
                if (division(absValue(gradient[3]), absValue(gradient_matrix[j][3])) > 1.05) gradient[3] = 1.05 * absValue(gradient_matrix[j][3]) * signDouble(gradient[3]);
                Psi_list_aux[j] = Psi_list_aux[j] - step * (alpha * gradient[3] + (1 - alpha) * gradient_matrix[j][3]);
                
                if (division(absValue(gradient[4]), absValue(gradient_matrix[j][4])) > 1.05) gradient[4] = 1.05 * absValue(gradient_matrix[j][4]) * signDouble(gradient[4]);
                amp_list_aux[j] = amp_list_aux[j] - step * (alpha * gradient[4] + (1 - alpha) * gradient_matrix[j][4]);

                if (delta_list_aux[j] < 1e-8) delta_list_aux[j] = 1e-10;
                if (amp_list_aux[j] < 0.0) amp_list_aux[j] = 0.0;

                // for (k = 0; k < faces_connection[j].n; k++) {

                //     if (division(absValue(gradient[5 + 5 * k]), absValue(gradient_matrix[j][5 + 5 * k])) > 1.05) gradient[5 + 5 * k] = 1.05 * absValue(gradient_matrix[j][5 + 5 * k]) * signDouble(gradient[5 + 5 * k]);
                //     delta_list_aux[faces_connection[j].faces[k]] = delta_list_aux[faces_connection[j].faces[k]] - step * (alpha * gradient[5 + 5 * k] + (1 - alpha) * gradient_matrix[j][5 + 5 * k]);

                //     if (division(absValue(gradient[5 + 5 * k + 1]), absValue(gradient_matrix[j][5 + 5 * k + 1])) > 1.05) gradient[5 + 5 * k + 1] = 1.05 * absValue(gradient_matrix[j][5 + 5 * k + 1]) * signDouble(gradient[5 + 5 * k + 1]);
                //     A_list_aux[faces_connection[j].faces[k]] = A_list_aux[faces_connection[j].faces[k]] - step * (alpha * gradient[5 + 5 * k + 1] + (1 - alpha) * gradient_matrix[j][5 + 5 * k + 1]);

                //     if (division(absValue(gradient[5 + 5 * k + 2]), absValue(gradient_matrix[j][5 + 5 * k + 2])) > 1.05) gradient[5 + 5 * k + 2] = 1.05 * absValue(gradient_matrix[j][5 + 5 * k + 2]) * signDouble(gradient[5 + 5 * k + 2]);
                //     B_list_aux[faces_connection[j].faces[k]] = B_list_aux[faces_connection[j].faces[k]] - step * (alpha * gradient[5 + 5 * k + 2] + (1 - alpha) * gradient_matrix[j][5 + 5 * k + 2]);
                    
                //     if (division(absValue(gradient[5 + 5 * k + 3]), absValue(gradient_matrix[j][5 + 5 * k + 3])) > 1.05) gradient[5 + 5 * k + 3] = 1.05 * absValue(gradient_matrix[j][5 + 5 * k + 3]) * signDouble(gradient[5 + 5 * k + 3]);
                //     Psi_list_aux[faces_connection[j].faces[k]] = Psi_list_aux[faces_connection[j].faces[k]] - step * (alpha * gradient[5 + 5 * k + 3] + (1 - alpha) * gradient_matrix[j][5 + 5 * k + 3]);
                    
                //     if (division(absValue(gradient[5 + 5 * k + 4]), absValue(gradient_matrix[j][5 + 5 * k + 4])) > 1.05) gradient[5 + 5 * k + 4] = 1.05 * absValue(gradient_matrix[j][5 + 5 * k + 4]) * signDouble(gradient[5 + 5 * k + 4]);
                //     amp_list_aux[faces_connection[j].faces[k]] = amp_list_aux[faces_connection[j].faces[k]] - step * (alpha * gradient[5 + 5 * k + 4] + (1 - alpha) * gradient_matrix[j][5 + 5 * k + 4]);

                //     if (delta_list_aux[faces_connection[j].faces[k]] < 1e-8) delta_list_aux[faces_connection[j].faces[k]] = 1e-10;
                //     if (amp_list_aux[faces_connection[j].faces[k]] < 0.0) amp_list_aux[faces_connection[j].faces[k]] = 0.0;

                // }
            
            }

            /* Error */
            if (j == 0) {
                previousMaxFaceEquations = faceEquations;
            } else if (faceEquations.obj > previousMaxFaceEquations.obj) {
                previousMaxFaceEquations = faceEquations;
            }

            /* Assign gradient */
            gradient_matrix[j][0] = gradient[0];
            gradient_matrix[j][1] = gradient[1];
            gradient_matrix[j][2] = gradient[2];
            gradient_matrix[j][3] = gradient[3];
            gradient_matrix[j][4] = gradient[4];

            for (k = 0; k < faces_connection[j].n; k++) {
                gradient_matrix[j][5 + 5 * k] = gradient[5 + 5 * k];
                gradient_matrix[j][5 + 5 * k + 1] = gradient[5 + 5 * k + 1];
                gradient_matrix[j][5 + 5 * k + 2] = gradient[5 + 5 * k + 2];
                gradient_matrix[j][5 + 5 * k + 3] = gradient[5 + 5 * k + 3];
                gradient_matrix[j][5 + 5 * k + 4] = gradient[5 + 5 * k + 4];
            }

            /* Free */
            free(gradient);

        }

        /* Calculate new error */
        // for (j = 0; j < nf; j++) {
        //     calculateParameters(j, delta_norm * delta_list_aux[j], A_norm * A_list_aux[j], B_norm * B_list_aux[j], Psi_norm * Psi_list_aux[j], amp_norm * amp_list_aux[j], velNorm[j], density, viscosity, mach[j], integralThickness);
        //     calculateParameters(j, delta_norm * (delta_list_aux[j] + eps), A_norm * A_list_aux[j], B_norm * B_list_aux[j], Psi_norm * Psi_list_aux[j], amp_norm * amp_list_aux[j], velNorm[j], density, viscosity, mach[j], integralThickness_delta);
        //     calculateParameters(j, delta_norm * delta_list_aux[j], A_norm * (A_list_aux[j] + eps), B_norm * B_list_aux[j], Psi_norm * Psi_list_aux[j], amp_norm * amp_list_aux[j], velNorm[j], density, viscosity, mach[j], integralThickness_A);
        //     calculateParameters(j, delta_norm * delta_list_aux[j], A_norm * A_list_aux[j], B_norm * (B_list_aux[j] + eps), Psi_norm * Psi_list_aux[j], amp_norm * amp_list_aux[j], velNorm[j], density, viscosity, mach[j], integralThickness_B);
        //     calculateParameters(j, delta_norm * delta_list_aux[j], A_norm * A_list_aux[j], B_norm * B_list_aux[j], Psi_norm * (Psi_list_aux[j] + eps), amp_norm * amp_list_aux[j], velNorm[j], density, viscosity, mach[j], integralThickness_Psi);
        //     calculateParameters(j, delta_norm * delta_list_aux[j], A_norm * A_list_aux[j], B_norm * B_list_aux[j], Psi_norm * Psi_list_aux[j], amp_norm * (amp_list_aux[j] + eps), velNorm[j], density, viscosity, mach[j], integralThickness_amp);
        // }

        // for (j = 0; j < nf; j++) {

        //     /* Calculate gradient */
        //     gradient = (double*)malloc(5 * (1 + faces_connection[j].n) * sizeof(double));

        //     calculateGradient(j, i, tau_wall_x, tau_wall_y, density, &faceEquations, eps, velNorm, faces, vertices, velx, vely, velz, transpiration, e3, vertices_connection, faces_connection, facesGradients, gradient, integralThickness, integralThickness_delta, integralThickness_A, integralThickness_B, integralThickness_Psi, integralThickness_amp);

        //     /* Error */
        //     if (j == 0) {
        //         maxFaceEquations = faceEquations;
        //     } else if (faceEquations.obj > maxFaceEquations.obj) {
        //         maxFaceEquations = faceEquations;
        //     }

        //     /* Free */
        //     free(gradient);

        // }

        /* Step */
        // if (division(maxFaceEquations.obj, previousMaxFaceEquations.obj) < 1.1) {
            
        //     for (j = 0; j < nf; j++) {
        //         delta_list[j] = delta_list_aux[j];
        //         A_list[j] =  A_list_aux[j];
        //         B_list[j] = B_list_aux[j];
        //         Psi_list[j] = Psi_list_aux[j];
        //         amp_list[j] = amp_list_aux[j];
        //     }

        //     if (division(maxFaceEquations.obj, previousMaxFaceEquations.obj) < 1.0) step_ref = step_ref * 1.005;

        // } else {

        //     step_ref = step_ref / 1.1;

        //     printf(" [%.3e]", step_ref);

        // }

        if (division(maxFaceEquations.obj, previousMaxFaceEquations.obj) < 1.0) step_ref = step_ref * 1.005;
        for (j = 0; j < nf; j++) {
            delta_list[j] = delta_list_aux[j];
            A_list[j] =  A_list_aux[j];
            B_list[j] = B_list_aux[j];
            Psi_list[j] = Psi_list_aux[j];
            amp_list[j] = amp_list_aux[j];
        }

        /* Calculate drag */
        drag.x = 0.0; drag.y = 0.0; drag.z = 0.0;

        for (j = 0; j < nf; j++) {
            drag.x = drag.x + facesAreas[j] * (s1[j].x * tau_wall_x[j] + s2[j].x * tau_wall_y[j]);
            drag.y = drag.y + facesAreas[j] * (s1[j].y * tau_wall_x[j] + s2[j].y * tau_wall_y[j]);
            drag.z = drag.z + facesAreas[j] * (s1[j].z * tau_wall_x[j] + s2[j].z * tau_wall_y[j]);
        }

        /* Print error */
        printf("    %d    %.3e %.3e %.3e %.3e %.3e %.3e %.3e %.3e\n", i, previousMaxFaceEquations.obj, freestream.x * drag.x + freestream.y * drag.y + freestream.z * drag.z, previousMaxFaceEquations.momentum_x, previousMaxFaceEquations.momentum_y, previousMaxFaceEquations.kinetic_energy, previousMaxFaceEquations.lateral_curvature, previousMaxFaceEquations.shear_stress_x, previousMaxFaceEquations.shear_stress_y);

        if (previousMaxFaceEquations.obj < 1e-8) break;
        if (previousMaxFaceEquations.obj > 1e2) break;

    }

    /* Assing values */
    for (i = 0; i < nf; i++) {
        delta[i] = delta_list[i] * delta_norm;
        A[i] = A_list[i] * A_norm;
        B[i] = B_list[i] * B_norm;
        Psi[i] = Psi_list[i] * Psi_norm;
        Ctau1[i] = amp_list[i] * amp_norm / sqrt(1.4);
        Ctau2[i] = Ctau1[i];
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

    printf("Aerodynamic solver\n");

    printf("  - Potential flow\n");

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
    printf("    > Creating linear system\n");
    createLinearSystem(nf, facesAreas, facesMaxDistance, facesCenter, controlPoints, p1, p2, p3, e1, e2, e3, freestream, sigma, nSpanLeftWing, nWakeLeftWing, leftWingGrid, leftWingVertices, leftWingFaces, nSpanRightWing, nWakeRightWing, rightWingGrid, rightWingVertices, rightWingFaces, nSpanTail, nWakeTail, tailGrid, tailVertices, tailFaces, matrix, array, matrixVelx, matrixVely, matrixVelz, arrayVel);

    /* Solve linear system with zero transpiration */
    printf("    > Solving linear system\n");
    calculateDoubletDistribution(nf, matrix, array, transpiration, doublet);

    /* Calculate potential surface parameters */
    calculateSurfaceParameters(nf, matrixVelx, matrixVely, matrixVelz, arrayVel, doublet, freestreamNorm, velx, vely, velz, velNorm, cp, mach, sound_speed);

    /* Boundary layer correction */

    /* Boundary layer parameters */
    double *tau_wall_x, *tau_wall_y;
    struct VerticeConnection *vertices_connection;
    struct Point freestream_point;

    /* Initialize */
    tau_wall_x = (double *)calloc(nf, sizeof(double));
    tau_wall_y = (double *)calloc(nf, sizeof(double));
    vertices_connection = (struct VerticeConnection *)malloc(nv * sizeof(struct VerticeConnection));
    freestream_point.x = freestream[0];
    freestream_point.y = freestream[1];
    freestream_point.z = freestream[2];

    calculateVerticesConnection(nv, nf, vertices, faces, vertices_connection);

    if (type == 1)
    {
        printf("  - Boundary layer correction\n");
        solveBoundaryLayer(nf, nv, freestream_point, vertices_connection, vertices, faces, facesCenter, facesAreas, e1, e2, e3, p1, p2, p3, transpiration, delta, A, B, Psi, Ctau1, Ctau2, tau_wall_x, tau_wall_y, velNorm, velx, vely, velz, mach, density, viscosity, cp, matrix, array, matrixVelx, matrixVely, matrixVelz, arrayVel, doublet, freestreamNorm);
        printf("\n");
    }

    /* Parameters */
    int i, j;
    struct Point e3_point, vel_point, s1, s2;
    double aux;

    // Surface shear stress
    for (i = 0; i < nf; i++) {

        e3_point.x = e3[3 * i];
        e3_point.y = e3[3 * i + 1];
        e3_point.z = e3[3 * i + 2];

        vel_point.x = velx[i] - transpiration[i] * e3_point.x;
        vel_point.y = vely[i] - transpiration[i] * e3_point.y;
        vel_point.z = velz[i] - transpiration[i] * e3_point.z;

        aux = norm(vel_point);

        s1.x = vel_point.x / aux;
        s1.y = vel_point.y / aux;
        s1.z = vel_point.z / aux;

        s2 = cross(s1, e3_point);

        tau_x[i] = s1.x * tau_wall_x[i] + s2.x * tau_wall_y[i];
        tau_y[i] = s1.y * tau_wall_x[i] + s2.y * tau_wall_y[i];
        tau_z[i] = s1.z * tau_wall_x[i] + s2.z * tau_wall_y[i];
    }

    // Vertices values
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
            delta_v[i] = delta_v[i] + delta[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            A_v[i] = A_v[i] + A[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            B_v[i] = B_v[i] + B[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            Psi_v[i] = Psi_v[i] + Psi[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            Ctau1_v[i] = Ctau1_v[i] + Ctau1[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            Ctau2_v[i] = Ctau2_v[i] + Ctau2[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            tau_x_v[i] = tau_x_v[i] + tau_x[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            tau_y_v[i] = tau_y_v[i] + tau_y[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
            tau_z_v[i] = tau_z_v[i] + tau_z[vertices_connection[i].faces[j]] * vertices_connection[i].coeffs[j];
        }
    
    }

    free(matrix);
    free(array);
    free(matrixVelx);
    free(matrixVely);
    free(matrixVelz);
    free(arrayVel);
}