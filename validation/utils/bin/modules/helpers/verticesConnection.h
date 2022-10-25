#ifndef VERTICES_CONNECTION_H
#define VERTICES_CONNECTION_H

#include "structs.h"

void getVerticesConnection(struct Input input, struct VerticesConnection *verticesConnection);
void getVerticesValues(struct Input input, struct VerticesConnection *verticesConnection, double *facesvalues, double *verticesValues);

#include "verticesConnection.c"

#endif