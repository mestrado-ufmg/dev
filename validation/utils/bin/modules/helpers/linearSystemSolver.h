#ifndef LINEAR_SYSTEM_SOLVER_H
#define LINEAR_SYSTEM_SOLVER_H

/*
    Solves a sparse linear system using least square error
    in the form of a triplet form representation.

    Parameters
    ----------
    - n = matrix size (length of rhs)
    - na = entries size (length of a)
    - a(k) = value of entry [left hand side]
    - ia(k) = row of entry
    - ja(k) = column of entry
    - rhs = right hand side of the equation
    - x = solution
*/
void solveGMRES(int n, int na, double *a, int *ia, int *ja, double *rhs, double *x);

#include "linearSystemSolver.c"

#endif