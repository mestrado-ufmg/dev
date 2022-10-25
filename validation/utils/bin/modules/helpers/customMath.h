#ifndef CUSTOM_MATH_H
#define CUSTOM_MATH_H

#include "constants.h"
#include "structs.h"
#include <math.h>

double division(double a, double b);
double norm(struct Point p);
struct Point cross(struct Point p1, struct Point p2);
double dot(struct Point p1, struct Point p2);
double angleBetweenVectors(struct Point p1, struct Point p2);

#include "customMath.c"

#endif