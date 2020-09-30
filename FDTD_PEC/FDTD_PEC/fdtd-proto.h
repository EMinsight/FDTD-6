#include "fdtd-grid.h"

double ezInc(double time, double location);
void gridInit(Grid* g);
void updateE(Grid* g);
void updateH(Grid* g);
void snapshot3d(Grid* g);
