#include "../helpers/structs.h"
#include "data.h"

void getSurfaceParametersImp(struct Input input, struct PotentialFlowData data)
{
    
    int i, j;

    for (i = 0; i < input.mesh.surface.nf; i++)
    {

        data.vel_x[i] = data.rhs_vel_x[i];
        data.vel_y[i] = data.rhs_vel_y[i];
        data.vel_z[i] = data.rhs_vel_z[i];

        for (j = 0; j < input.mesh.surface.nf; j++)
        {
            data.vel_x[i] = data.vel_x[i] + data.a_vel_x[i * input.mesh.surface.nf + j] * data.doublet[j];
            data.vel_y[i] = data.vel_y[i] + data.a_vel_y[i * input.mesh.surface.nf + j] * data.doublet[j];
            data.vel_z[i] = data.vel_z[i] + data.a_vel_z[i * input.mesh.surface.nf + j] * data.doublet[j];
        }

        data.cp[i] = 1 - (pow(data.vel_x[i], 2) + pow(data.vel_y[i], 2) + pow(data.vel_z[i], 2)) / pow(input.environment.velNorm, 2);
        
        data.transpiration[i] = data.vel_x[i] * input.mesh.surface.e3[3 * i] + data.vel_y[i] * input.mesh.surface.e3[3 * i + 1] + data.vel_z[i] * input.mesh.surface.e3[3 * i + 2];
    
    }

}