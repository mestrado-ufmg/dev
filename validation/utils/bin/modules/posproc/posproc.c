#include "posproc.h"
#include "../helpers/structs.h"

void getForces(struct Input input, double *cp, double *forces)
{
    
    forces[0] = 0.0;
    forces[1] = 0.0;
    forces[2] = 0.0;

    for (int i = 0; i < input.mesh.surface.nf; i++)
    {
        forces[0] = forces[0] - cp[i] * input.mesh.surface.facesAreas[i] * input.mesh.surface.e3[3 * i];
        forces[1] = forces[1] - cp[i] * input.mesh.surface.facesAreas[i] * input.mesh.surface.e3[3 * i + 1];
        forces[2] = forces[2] - cp[i] * input.mesh.surface.facesAreas[i] * input.mesh.surface.e3[3 * i + 2];
    }

    forces[0] = forces[0] * 0.5 * input.environment.density * input.environment.velNorm * input.environment.velNorm;
    forces[1] = forces[1] * 0.5 * input.environment.density * input.environment.velNorm * input.environment.velNorm;
    forces[2] = forces[2] * 0.5 * input.environment.density * input.environment.velNorm * input.environment.velNorm;

}