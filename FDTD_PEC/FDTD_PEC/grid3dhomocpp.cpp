/***********************************************************/
/*                对网格参数进行初始化                     */
/***********************************************************/
#include "fdtd-macro.h"
#include "fdtd-alloc.h"

#include <math.h>

extern int sizeX;
extern int sizeY;
extern int sizeZ;
extern int maxTime;

extern double dx;
extern double dy;
extern double dz;
extern double dt;
extern double c0;
extern double mu0;
extern double epsR;
extern double eps0;

void gridInit(Grid* g) {
    int mm, nn, pp;
    int ii, jj, kk;

    SizeX = sizeX;
    SizeY = sizeY;
    SizeZ = sizeZ;
    MaxTime = maxTime;

    ALLOC_3D(g->hx, SizeX, SizeY - 1, SizeZ - 1, double);
    ALLOC_3D(g->chxh, SizeX, SizeY, SizeZ, double);
    ALLOC_3D(g->chxe, SizeX, SizeY, SizeZ, double);
    ALLOC_3D(g->hy, SizeX - 1, SizeY, SizeZ - 1, double);
    ALLOC_3D(g->chyh, SizeX, SizeY, SizeZ, double);
    ALLOC_3D(g->chye, SizeX, SizeY, SizeZ, double);
    ALLOC_3D(g->hz, SizeX - 1, SizeY - 1, SizeZ, double);
    ALLOC_3D(g->chzh, SizeX, SizeY, SizeZ, double);
    ALLOC_3D(g->chze, SizeX, SizeY, SizeZ, double);

    ALLOC_3D(g->ex, SizeX - 1, SizeY, SizeZ, double);
    ALLOC_3D(g->cexe, SizeX, SizeY, SizeZ, double);
    ALLOC_3D(g->cexh, SizeX, SizeY, SizeZ, double);
    ALLOC_3D(g->ey, SizeX, SizeY - 1, SizeZ, double);
    ALLOC_3D(g->ceye, SizeX, SizeY, SizeZ, double);
    ALLOC_3D(g->ceyh, SizeX, SizeY, SizeZ, double);
    ALLOC_3D(g->ez, SizeX, SizeY, SizeZ - 1, double);
    ALLOC_3D(g->ceze, SizeX, SizeY, SizeZ, double);
    ALLOC_3D(g->cezh, SizeX, SizeY, SizeZ, double);

    ALLOC_3D(g->sig, SizeX, SizeY, SizeZ, double);
    ALLOC_3D(g->eps, SizeX, SizeY, SizeZ, double);

    for (mm = 0; mm < SizeX; mm++) {
        for (nn = 0; nn < SizeY; nn++) {
            for (pp = 0; pp < SizeZ; pp++) {
                Sig(mm, nn, pp) = 0.0;
                Eps(mm, nn, pp) = epsR * eps0;
            }
        }
    }
    for (mm = 0; mm < SizeX; mm++) {
        for (nn = 0; nn < SizeY; nn++) {
            for (pp = 0; pp < SizeZ; pp++) {
                Cexe(mm, nn, pp) = (1.0 - Sig(mm, nn, pp) * dt / (2.0 * Eps(mm, nn, pp))) / (1.0 + Sig(mm, nn, pp) * dt / (2.0 * Eps(mm, nn, pp)));
                Cexh(mm, nn, pp) = (dt / Eps(mm, nn, pp)) / (1.0 + Sig(mm, nn, pp) * dt / (2.0 * Eps(mm, nn, pp)));
                Ceye(mm, nn, pp) = (1.0 - Sig(mm, nn, pp) * dt / (2.0 * Eps(mm, nn, pp))) / (1.0 + Sig(mm, nn, pp) * dt / (2.0 * Eps(mm, nn, pp)));
                Ceyh(mm, nn, pp) = (dt / Eps(mm, nn, pp)) / (1.0 + Sig(mm, nn, pp) * dt / (2.0 * Eps(mm, nn, pp)));
                Ceze(mm, nn, pp) = (1.0 - Sig(mm, nn, pp) * dt / (2.0 * Eps(mm, nn, pp))) / (1.0 + Sig(mm, nn, pp) * dt / (2.0 * Eps(mm, nn, pp)));
                Cezh(mm, nn, pp) = (dt / Eps(mm, nn, pp)) / (1.0 + Sig(mm, nn, pp) * dt / (2.0 * Eps(mm, nn, pp)));
                Chxh(mm, nn, pp) = 1.0;
                Chxe(mm, nn, pp) = dt / mu0;
                Chyh(mm, nn, pp) = 1.0;
                Chye(mm, nn, pp) = dt / mu0;
                Chzh(mm, nn, pp) = 1.0;
                Chze(mm, nn, pp) = dt / mu0;
            }
        }
    }
}


