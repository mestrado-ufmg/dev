#include "potentialFlow.h"
#include "getLinearSystemImp.c"
#include "getDoubletDistributionImp.c"
#include "getSurfaceParametersImp.c"

void getLinearSystem(struct Input input, struct PotentialFlowData data)
{
    getLinearSystemImp(input, data);
}

void getDoubleDistribution(struct PotentialFlowData data)
{
    getDoubleDistributionImp(data);
}

void getSurfaceParameters(struct Input input, struct PotentialFlowData data)
{
    getSurfaceParametersImp(input, data);
}