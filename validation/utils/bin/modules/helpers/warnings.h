#ifndef WARNINGS_H
#define WARNINGS_H

#include <stdio.h>

void warnings(int step)
{
    if (step == 1) printf("  > Potential flow\n");
    if (step == 2) printf("    * Creating linear system\n");
    if (step == 3) printf("    * Solving linear system\n");
    if (step == 4) printf("    * Pos processing\n");
    if (step == 5) printf("  > Calculating vertices values\n");
}

#endif