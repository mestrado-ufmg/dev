#ifndef POTENTIAL_FLOW_H
#define POTENTIAL_FLOW_H

#include "../helpers/structs.h"
#include "data.h"

void getLinearSystem(struct Input input, struct PotentialFlowData data);
void getDoubleDistribution(struct PotentialFlowData data);
void getSurfaceParameters(struct Input input, struct PotentialFlowData data);

#include "potentialFlow.c"

#endif