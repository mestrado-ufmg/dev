#include "customMath.h"
#include "structs.h"
#include <math.h>

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