#include "verticesConnection.h"
#include "structs.h"
#include "customMath.h"

void getVerticesConnection(struct Input input, struct VerticesConnection *verticesConnection)
{

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
    for (i = 0; i < input.mesh.surface.nv; i++) {

        // Reset the number of faces and angles
        n = 0;
        sum = 0.0;

        /* Loop over faces */
        for (j = 0; j < input.mesh.surface.nf; j++) {

            faceLine = j * 3;

            /* Check if the face contain the vertice */
            if ((input.mesh.surface.faces[faceLine] == i) || (input.mesh.surface.faces[faceLine + 1] == i) || (input.mesh.surface.faces[faceLine + 2] == i)) {

                /* Calculate the angle */
                if (input.mesh.surface.faces[faceLine] == i)
                {
                    verticeLine1 = 3 * input.mesh.surface.faces[faceLine];
                    verticeLine2 = 3 * input.mesh.surface.faces[faceLine + 1];
                    verticeLine3 = 3 * input.mesh.surface.faces[faceLine + 2];
                }
                else if (input.mesh.surface.faces[faceLine + 1] == i)
                {
                    verticeLine3 = 3 * input.mesh.surface.faces[faceLine];
                    verticeLine1 = 3 * input.mesh.surface.faces[faceLine + 1];
                    verticeLine2 = 3 * input.mesh.surface.faces[faceLine + 2];
                }
                else
                {
                    verticeLine2 = 3 * input.mesh.surface.faces[faceLine];
                    verticeLine3 = 3 * input.mesh.surface.faces[faceLine + 1];
                    verticeLine1 = 3 * input.mesh.surface.faces[faceLine + 2];
                }

                point1.x = input.mesh.surface.vertices[verticeLine2] - input.mesh.surface.vertices[verticeLine1];
                point1.y = input.mesh.surface.vertices[verticeLine2 + 1] - input.mesh.surface.vertices[verticeLine1 + 1];
                point1.z = input.mesh.surface.vertices[verticeLine2 + 2] - input.mesh.surface.vertices[verticeLine1 + 2];

                point2.x = input.mesh.surface.vertices[verticeLine3] - input.mesh.surface.vertices[verticeLine1];
                point2.y = input.mesh.surface.vertices[verticeLine3 + 1] - input.mesh.surface.vertices[verticeLine1 + 1];
                point2.z = input.mesh.surface.vertices[verticeLine3 + 2] - input.mesh.surface.vertices[verticeLine1 + 2];

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
        verticesConnection[i].n = n;
        verticesConnection[i].coeffs = (double *)malloc(n * sizeof(double));
        verticesConnection[i].faces = (int *)malloc(n * sizeof(int));

        for (k = 0; k < n; k++)
        {
            verticesConnection[i].coeffs[k] = angles[k] / sum;
            verticesConnection[i].faces[k] = facesIds[k];
        }

    }

    free(facesIds);
    free(angles);

}

void getVerticesValues(struct Input input, struct VerticesConnection *verticesConnection, double *facesvalues, double *verticesValues)
{

    int i, j;

    for (i = 0; i < input.mesh.surface.nv; i++)
    {
        
        verticesValues[i] = 0.0;
        
        for (j = 0; j < verticesConnection[i].n; j++)
        {
            verticesValues[i] = verticesValues[i] + facesvalues[verticesConnection[i].faces[j]] * verticesConnection[i].coeffs[j];
        }

    }

}
