#include "../helpers/linearSystemSolver.h"
#include "data.h"

void getDoubleDistributionImp(struct PotentialFlowData data)
{
    for (int i = 0; i < data.n; i++) data.doublet[i] = 0.0;
    solveGMRES(data.n, data.na, data.a, data.ia, data.ja, data.rhs, data.doublet);
}