
#include "fdtd-grid.h"

void gridInit(Grid* g, int isSpherePresent, bool isAid);
void gridInit1d(Grid* g);
void gridInit2d(Grid* g);

void tfsfInit(Grid* g);
void tfsf(Grid* g);
void addSource(Grid* g);
double ezInc(double time, double location);
double getDist(int x, int y, int z, int inDirection);
void update_tfsf_H(Grid* g);
void update_tfsf_E(Grid* g);

double ezInc(double time, double location);

void updateE(Grid* g, bool isAid);
void updateH(Grid* g, bool isAid);

void snapshot3d(Grid* g);
