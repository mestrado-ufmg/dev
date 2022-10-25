#include "./modules/helpers/warnings.h"
#include "./modules/helpers/structs.h"
#include "./modules/helpers/verticesConnection.h"
#include "./modules/potentialFlow/potentialFlow.h"
#include "./modules/potentialFlow/data.h"
#include "./modules/posproc/posproc.h"

void solve(struct Input input, double *vel_x_v, double *vel_y_v, double *vel_z_v, double *transpiration_v, double *sigma_v, double *doublet_v, double *forces)
{   

    /* Parameters */
    struct PotentialFlowData potentialFlowData = getPotentialFlowData(input.mesh.surface.nf, input.mesh.surface.e3, input.environment.vel_x, input.environment.vel_y, input.environment.vel_z);
    struct VerticesConnection *verticesConnetion = getVerticesConnectionData(input.mesh.surface.nv);

    /* Potential flow */
    warnings(1);
    warnings(2);
    getLinearSystem(input, potentialFlowData);
    warnings(3);
    getDoubleDistribution(potentialFlowData);
    warnings(4);
    getSurfaceParameters(input, potentialFlowData);
    getForces(input, potentialFlowData.cp, forces);

    /* Vertices values */
    warnings(5);
    getVerticesConnection(input, verticesConnetion);
    
    getVerticesValues(input, verticesConnetion, potentialFlowData.vel_x, vel_x_v);
    getVerticesValues(input, verticesConnetion, potentialFlowData.vel_y, vel_y_v);
    getVerticesValues(input, verticesConnetion, potentialFlowData.vel_z, vel_z_v);
    getVerticesValues(input, verticesConnetion, potentialFlowData.sigma, sigma_v);
    getVerticesValues(input, verticesConnetion, potentialFlowData.doublet, doublet_v);
    getVerticesValues(input, verticesConnetion, potentialFlowData.transpiration, transpiration_v);

    /* Free */
    // freePotentialFlowData(potentialFlowData);

}