#include "fdtd-macro.h"
#include "fdtd-alloc.h"
#include <math.h>

#define NLOSS     20
#define MAX_LOSS  0.35

extern int m;
extern int ma;

extern int xPML1;
extern int xPML2;
extern int yPML1;
extern int yPML2;
extern int zPML1;
extern int zPML2;

extern int sizeX;
extern int sizeY;
extern int sizeZ;
extern int maxTime;
extern double radius;

extern double sig_x;
extern double sig_y;
extern double sig_z;
extern double alp_x;
extern double alp_y;
extern double alp_z;
extern double kap_x;
extern double kap_y;
extern double kap_z;

extern double dx;
extern double dy;
extern double dz;
extern double dt;
extern double c0;
extern double mu0;
extern double epsR;
extern double eps0;
void gridInit1d(Grid* g)
{
    double imp0 = 377.0, depthInLayer, lossFactor;
    int mm;

    SizeX = 2 * (int)floor(sqrt(sizeX * sizeX + sizeY * sizeY + sizeZ * sizeZ)) + 1;
    SizeX += NLOSS;
    Type = oneDGrid;

    ALLOC_1D(g->hy, SizeX, double);
    ALLOC_1D(g->chyh, SizeX, double);
    ALLOC_1D(g->chye, SizeX, double);
    ALLOC_1D(g->ez, SizeX, double);
    ALLOC_1D(g->ceze, SizeX, double);
    ALLOC_1D(g->cezh, SizeX, double);


    for (mm = 0; mm < SizeX; mm++) {
        if (mm < SizeX - NLOSS) {
            Ceze1(mm) = 1.0;
            Cezh1(mm) = Cdtds * imp0;
            Chyh1(mm) = 1.0;
            Chye1(mm) = Cdtds / imp0;
        }
        else {
            depthInLayer = mm - (SizeX - NLOSS) + 0.5;
            lossFactor = MAX_LOSS * pow(depthInLayer / NLOSS, 2);
            Ceze1(mm) = (1.0 - lossFactor) / (1.0 + lossFactor);
            Cezh1(mm) = Cdtds * imp0 / (1.0 + lossFactor);
            depthInLayer += 0.5;
            lossFactor = MAX_LOSS * pow(depthInLayer / NLOSS, 2);
            Chyh1(mm) = (1.0 - lossFactor) / (1.0 + lossFactor);
            Chye1(mm) = Cdtds / imp0 / (1.0 + lossFactor);
        }
    }
    return;
}

void gridInit2d(Grid* g)
{
    double imp0 = 377.0;
    int mm, nn;


    Type = tmZGrid;

    ALLOC_2D(g->hx, SizeX, SizeY - 1, double);
    ALLOC_2D(g->chxh, SizeX, SizeY - 1, double);
    ALLOC_2D(g->chxe, SizeX, SizeY - 1, double);
    ALLOC_2D(g->hy, SizeX - 1, SizeY, double);
    ALLOC_2D(g->chyh, SizeX - 1, SizeY, double);
    ALLOC_2D(g->chye, SizeX - 1, SizeY, double);
    ALLOC_2D(g->ez, SizeX, SizeY, double);
    ALLOC_2D(g->ceze, SizeX, SizeY, double);
    ALLOC_2D(g->cezh, SizeX, SizeY, double);


    for (mm = 0; mm < SizeX; mm++) {
        for (nn = 0; nn < SizeY; nn++) {
            Ceze2(mm, nn) = 1.0;
            Cezh2(mm, nn) = Cdtds * imp0;
        }
    }
    for (mm = 0; mm < SizeX; mm++) {
        for (nn = 0; nn < SizeY - 1; nn++) {
            Chxh2(mm, nn) = 1.0;
            Chxe2(mm, nn) = Cdtds / imp0;
        }
    }
    for (mm = 0; mm < SizeX - 1; mm++) {
        for (nn = 0; nn < SizeY; nn++) {
            Chyh2(mm, nn) = 1.0;
            Chye2(mm, nn) = Cdtds / imp0;
        }
    }
    return;
}
void gridInit(Grid* g, int isSpherePresent, bool isAid) {
    double imp0 = 377.0;
    int mm, nn, pp;
    int ii, jj, kk;

    int m_c = 60, n_c = 60, p_c = 60;
    double m2, n2, p2, r2;

    Type = threeGrid;
    SizeX = sizeX;
    SizeY = sizeY;
    SizeZ = sizeZ;
    MaxTime = maxTime;
    Cdtds = 1.0 / sqrt(3.0);

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

    ALLOC_1D(g->den_ex, SizeX - 1, double);
    ALLOC_1D(g->den_hx, SizeX - 1, double);
    ALLOC_1D(g->den_ey, SizeY - 1, double);
    ALLOC_1D(g->den_hy, SizeY - 1, double);
    ALLOC_1D(g->den_ez, SizeZ - 1, double);
    ALLOC_1D(g->den_hz, SizeZ - 1, double);

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
    ALLOC_3D(g->p_Exy1, SizeX - 1, yPML1, SizeZ, double);
    ALLOC_3D(g->p_Exy2, SizeX - 1, yPML2, SizeZ, double);
    ALLOC_3D(g->p_Exz1, SizeX - 1, SizeY, zPML1, double);
    ALLOC_3D(g->p_Exz2, SizeX - 1, SizeY, zPML2, double);
    ALLOC_3D(g->p_Eyz1, SizeX, SizeY - 1, zPML1, double);
    ALLOC_3D(g->p_Eyz2, SizeX, SizeY - 1, zPML2, double);
    ALLOC_3D(g->p_Eyx1, xPML1, SizeY - 1, SizeZ, double);
    ALLOC_3D(g->p_Eyx2, xPML2, SizeY - 1, SizeZ, double);
    ALLOC_3D(g->p_Ezx1, xPML1, SizeY, SizeZ - 1, double);
    ALLOC_3D(g->p_Ezx2, xPML2, SizeY, SizeZ - 1, double);
    ALLOC_3D(g->p_Ezy1, SizeX, yPML1, SizeZ - 1, double);
    ALLOC_3D(g->p_Ezy2, SizeX, yPML2, SizeZ - 1, double);

    ALLOC_3D(g->p_Hxy1, SizeX, yPML1 - 1, SizeZ - 1, double);
    ALLOC_3D(g->p_Hxy2, SizeX, yPML2 - 1, SizeZ - 1, double);
    ALLOC_3D(g->p_Hxz1, SizeX, SizeY - 1, zPML1 - 1, double);
    ALLOC_3D(g->p_Hxz2, SizeX, SizeY - 1, zPML2 - 1, double);
    ALLOC_3D(g->p_Hyz1, SizeX - 1, SizeY, zPML1 - 1, double);
    ALLOC_3D(g->p_Hyz2, SizeX - 1, SizeY, zPML2 - 1, double);
    ALLOC_3D(g->p_Hyx1, xPML1 - 1, SizeY, SizeZ - 1, double);
    ALLOC_3D(g->p_Hyx2, xPML2 - 1, SizeY, SizeZ - 1, double);
    ALLOC_3D(g->p_Hzx1, xPML1 - 1, SizeY - 1, SizeZ, double);
    ALLOC_3D(g->p_Hzx2, xPML2 - 1, SizeY - 1, SizeZ, double);
    ALLOC_3D(g->p_Hzy1, SizeX - 1, yPML1 - 1, SizeZ, double);
    ALLOC_3D(g->p_Hzy2, SizeX - 1, yPML2 - 1, SizeZ, double);

    ALLOC_1D(g->be_x1, xPML1, double);
    ALLOC_1D(g->ce_x1, xPML1, double);
    ALLOC_1D(g->alpe_x1, xPML1, double);
    ALLOC_1D(g->sige_x1, xPML1, double);
    ALLOC_1D(g->kape_x1, xPML1, double);
    ALLOC_1D(g->be_x2, xPML2, double);
    ALLOC_1D(g->ce_x2, xPML2, double);
    ALLOC_1D(g->alpe_x2, xPML2, double);
    ALLOC_1D(g->sige_x2, xPML2, double);
    ALLOC_1D(g->kape_x2, xPML2, double);

    ALLOC_1D(g->be_y1, yPML1, double);
    ALLOC_1D(g->ce_y1, yPML1, double);
    ALLOC_1D(g->alpe_y1, yPML1, double);
    ALLOC_1D(g->sige_y1, yPML1, double);
    ALLOC_1D(g->kape_y1, yPML1, double);
    ALLOC_1D(g->be_y2, yPML2, double);
    ALLOC_1D(g->ce_y2, yPML2, double);
    ALLOC_1D(g->alpe_y2, yPML2, double);
    ALLOC_1D(g->sige_y2, yPML2, double);
    ALLOC_1D(g->kape_y2, yPML2, double);

    ALLOC_1D(g->be_z1, zPML1, double);
    ALLOC_1D(g->ce_z1, zPML1, double);
    ALLOC_1D(g->alpe_z1, zPML1, double);
    ALLOC_1D(g->sige_z1, zPML1, double);
    ALLOC_1D(g->kape_z1, zPML1, double);
    ALLOC_1D(g->be_z2, zPML2, double);
    ALLOC_1D(g->ce_z2, zPML2, double);
    ALLOC_1D(g->alpe_z2, zPML2, double);
    ALLOC_1D(g->sige_z2, zPML2, double);
    ALLOC_1D(g->kape_z2, zPML2, double);

    ALLOC_1D(g->bh_x1, xPML1 - 1, double);
    ALLOC_1D(g->ch_x1, xPML1 - 1, double);
    ALLOC_1D(g->alph_x1, xPML1 - 1, double);
    ALLOC_1D(g->sigh_x1, xPML1 - 1, double);
    ALLOC_1D(g->kaph_x1, xPML1 - 1, double);
    ALLOC_1D(g->bh_x2, xPML2 - 1, double);
    ALLOC_1D(g->ch_x2, xPML2 - 1, double);
    ALLOC_1D(g->alph_x2, xPML2 - 1, double);
    ALLOC_1D(g->sigh_x2, xPML2 - 1, double);
    ALLOC_1D(g->kaph_x2, xPML2 - 1, double);

    ALLOC_1D(g->bh_y1, yPML1 - 1, double);
    ALLOC_1D(g->ch_y1, yPML1 - 1, double);
    ALLOC_1D(g->alph_y1, yPML1 - 1, double);
    ALLOC_1D(g->sigh_y1, yPML1 - 1, double);
    ALLOC_1D(g->kaph_y1, yPML1 - 1, double);
    ALLOC_1D(g->bh_y2, yPML2 - 1, double);
    ALLOC_1D(g->ch_y2, yPML2 - 1, double);
    ALLOC_1D(g->alph_y2, yPML2 - 1, double);
    ALLOC_1D(g->sigh_y2, yPML2 - 1, double);
    ALLOC_1D(g->kaph_y2, yPML2 - 1, double);

    ALLOC_1D(g->bh_z1, zPML1 - 1, double);
    ALLOC_1D(g->ch_z1, zPML1 - 1, double);
    ALLOC_1D(g->alph_z1, zPML1 - 1, double);
    ALLOC_1D(g->sigh_z1, zPML1 - 1, double);
    ALLOC_1D(g->kaph_z1, zPML1 - 1, double);
    ALLOC_1D(g->bh_z2, zPML2 - 1, double);
    ALLOC_1D(g->ch_z2, zPML2 - 1, double);
    ALLOC_1D(g->alph_z2, zPML2 - 1, double);
    ALLOC_1D(g->sigh_z2, zPML2 - 1, double);
    ALLOC_1D(g->kaph_z2, zPML2 - 1, double);



    //x dir
    for (mm = 0; mm < xPML1; mm++) {
        g->sige_x1[mm] = sig_x * pow((xPML1 - 1.0 - mm) / (xPML1 - 1.0), m);
        g->alpe_x1[mm] = alp_x * pow(mm / (xPML1 - 1.0), ma);
        g->kape_x1[mm] = 1.0 + (kap_x - 1.0) * pow((xPML1 - 1.0 - mm) / (xPML1 - 1.0), m);
        g->be_x1[mm] = exp(-(g->sige_x1[mm] / g->kape_x1[mm] + g->alpe_x1[mm]) * dt / eps0);
        if ((g->sige_x1[mm] == 0.0) && (g->alpe_x1[mm] == 0.0) && (mm + 1 == xPML1))
            g->ce_x1[mm] = 0.0;
        else
            g->ce_x1[mm] = g->sige_x1[mm] * (g->be_x1[mm] - 1.0) / (g->sige_x1[mm] + g->kape_x1[mm] * g->alpe_x1[mm]) / g->kape_x1[mm];
    }
    for (mm = 0; mm < xPML1 - 1; mm++) {
        g->sigh_x1[mm] = sig_x * pow((xPML1 - 1.5 - mm) / (xPML1 - 1.0), m);
        g->alph_x1[mm] = alp_x * pow((mm + 0.5) / (xPML1 - 1.0), ma);
        g->kaph_x1[mm] = 1.0 + (kap_x - 1.0) * pow((xPML1 - 1.5 - mm) / (xPML1 - 1.0), m);
        g->bh_x1[mm] = exp(-(g->sigh_x1[mm] / g->kaph_x1[mm] + g->alph_x1[mm]) / eps0 * dt);
        g->ch_x1[mm] = g->sigh_x1[mm] * (g->bh_x1[mm] - 1.0) / (g->sigh_x1[mm] + g->kaph_x1[mm] * g->alph_x1[mm]) / g->kaph_x1[mm];
    }
    for (mm = 0; mm < xPML2; mm++) {
        g->sige_x2[mm] = sig_x * pow((xPML2 - 1.0 - mm) / (xPML2 - 1.0), m);
        g->alpe_x2[mm] = alp_x * pow(mm / (xPML2 - 1.0), ma);
        g->kape_x2[mm] = 1.0 + (kap_x - 1.0) * pow((xPML2 - 1.0 - mm) / (xPML2 - 1.0), m);
        g->be_x2[mm] = exp(-(g->sige_x2[mm] / g->kape_x2[mm] + g->alpe_x2[mm]) / eps0 * dt);
        if ((g->sige_x2[mm] == 0.0) && (g->alpe_x2[mm] == 0.0) && (mm + 1 == xPML2))
            g->ce_x2[mm] = 0.0;
        else
            g->ce_x2[mm] = g->sige_x2[mm] * (g->be_x2[mm] - 1.0) / (g->sige_x2[mm] + g->kape_x2[mm] * g->alpe_x2[mm]) / g->kape_x2[mm];
    }
    for (mm = 0; mm < xPML2 - 1; mm++) {
        g->sigh_x2[mm] = sig_x * pow((xPML2 - 1.5 - mm) / (xPML2 - 1.0), m);
        g->alph_x2[mm] = alp_x * pow((mm + 0.5) / (xPML2 - 1.0), ma);
        g->kaph_x2[mm] = 1.0 + (kap_x - 1) * pow((xPML2 - 1.5 - mm) / (xPML2 - 1.0), m);
        g->bh_x2[mm] = exp(-(g->sigh_x2[mm] / g->kaph_x2[mm] + g->alph_x2[mm]) / eps0 * dt);
        g->ch_x2[mm] = g->sigh_x2[mm] * (g->bh_x2[mm] - 1.0) / (g->sigh_x2[mm] + g->kaph_x2[mm] * g->alph_x2[mm]) / g->kaph_x2[mm];
    }
    //y dir
    for (nn = 0; nn < yPML1; nn++) {
        g->sige_y1[nn] = sig_y * pow((yPML1 - 1.0 - nn) / (yPML1 - 1.0), m);
        g->alpe_y1[nn] = alp_y * pow(nn / (yPML1 - 1.0), ma);
        g->kape_y1[nn] = 1.0 + (kap_y - 1.0) * pow((yPML1 - 1.0 - nn) / (yPML1 - 1.0), m);
        g->be_y1[nn] = exp(-(g->sige_y1[nn] / g->kape_y1[nn] + g->alpe_y1[nn]) / eps0 * dt);
        if ((g->sige_y1[nn] == 0.0) && (g->alpe_y1[nn] == 0.0) && (nn + 1 == yPML1))
            g->ce_y1[nn] = 0.0;
        else
            g->ce_y1[nn] = g->sige_y1[nn] * (g->be_y1[nn] - 1.0) / (g->sige_y1[nn] + g->kape_y1[nn] * g->alpe_y1[nn]) / g->kape_y1[nn];
    }
    for (nn = 0; nn < yPML1 - 1; nn++) {
        g->sigh_y1[nn] = sig_y * pow((yPML1 - 1.5 - nn) / (yPML1 - 1.0), m);
        g->alph_y1[nn] = alp_y * pow((nn + 0.5) / (yPML1 - 1.0), ma);
        g->kaph_y1[nn] = 1.0 + (kap_y - 1.0) * pow((yPML1 - 1.5 - nn) / (yPML1 - 1.0), m);
        g->bh_y1[nn] = exp(-(g->sigh_y1[nn] / g->kaph_y1[nn] + g->alph_y1[nn]) / eps0 * dt);
        g->ch_y1[nn] = g->sigh_y1[nn] * (g->bh_y1[nn] - 1.0) / (g->sigh_y1[nn] + g->kaph_y1[nn] * g->alph_y1[nn]) / g->kaph_y1[nn];
    }
    for (nn = 0; nn < yPML2; nn++) {
        g->sige_y2[nn] = sig_y * pow((yPML2 - 1.0 - nn) / (yPML2 - 1.0), m);
        g->alpe_y2[nn] = alp_y * pow(nn / (yPML2 - 1.0), ma);
        g->kape_y2[nn] = 1.0 + (kap_y - 1.0) * pow((yPML2 - 1.0 - nn) / (yPML2 - 1.0), m);
        g->be_y2[nn] = exp(-(g->sige_y2[nn] / g->kape_y2[nn] + g->alpe_y2[nn]) / eps0 * dt);
        if ((g->sige_y2[nn] == 0.0) && (g->alpe_y2[nn] == 0.0) && (nn + 1 == yPML2))
            g->ce_y2[nn] = 0.0;
        else
            g->ce_y2[nn] = g->sige_y2[nn] * (g->be_y2[nn] - 1.0) / (g->sige_y2[nn] + g->kape_y2[nn] * g->alpe_y2[nn]) / g->kape_y2[nn];
    }
    for (nn = 0; nn < yPML2 - 1; nn++) {
        g->sigh_y2[nn] = sig_y * pow((yPML2 - 1.5 - nn) / (yPML2 - 1.0), m);
        g->alph_y2[nn] = alp_y * pow((nn + 0.5) / (yPML2 - 1.0), ma);
        g->kaph_y2[nn] = 1.0 + (kap_y - 1) * pow((yPML2 - 1.5 - nn) / (yPML2 - 1.0), m);
        g->bh_y2[nn] = exp(-(g->sigh_y2[nn] / g->kaph_y2[nn] + g->alph_y2[nn]) / eps0 * dt);
        g->ch_y2[nn] = g->sigh_y2[nn] * (g->bh_y2[nn] - 1.0) / (g->sigh_y2[nn] + g->kaph_y2[nn] * g->alph_y2[nn]) / g->kaph_y2[nn];
    }
    //z dir
    for (pp = 0; pp < zPML1; pp++) {
        g->sige_z1[pp] = sig_z * pow((zPML1 - 1.0 - pp) / (zPML1 - 1.0), m);
        g->alpe_z1[pp] = alp_z * pow(pp / (zPML1 - 1.0), ma);
        g->kape_z1[pp] = 1.0 + (kap_z - 1.0) * pow((zPML1 - 1.0 - pp) / (zPML1 - 1.0), m);
        g->be_z1[pp] = exp(-(g->sige_z1[pp] / g->kape_z1[pp] + g->alpe_z1[pp]) / eps0 * dt);
        if ((g->sige_z1[pp] == 0.0) && (g->alpe_z1[pp] == 0.0) && (pp + 1 == zPML1))
            g->ce_z1[pp] = 0.0;
        else
            g->ce_z1[pp] = g->sige_z1[pp] * (g->be_z1[pp] - 1.0) / (g->sige_z1[pp] + g->kape_z1[pp] * g->alpe_z1[pp]) / g->kape_z1[pp];
    }
    for (pp = 0; pp < zPML1 - 1; pp++) {
        g->sigh_z1[pp] = sig_z * pow((zPML1 - 1.5 - pp) / (zPML1 - 1.0), m);
        g->alph_z1[pp] = alp_z * pow((pp + 0.5) / (zPML1 - 1.0), ma);
        g->kaph_z1[pp] = 1.0 + (kap_z - 1) * pow((zPML1 - 1.5 - pp) / (zPML1 - 1.0), m);
        g->bh_z1[pp] = exp(-(g->sigh_z1[pp] / g->kaph_z1[pp] + g->alph_z1[pp]) / eps0 * dt);
        g->ch_z1[pp] = g->sigh_z1[pp] * (g->bh_z1[pp] - 1.0) / (g->sigh_z1[pp] + g->kaph_z1[pp] * g->alph_z1[pp]) / g->kaph_z1[pp];
    }
    for (pp = 0; pp < zPML2; pp++) {
        g->sige_z2[pp] = sig_z * pow((zPML2 - 1.0 - pp) / (zPML2 - 1.0), m);
        g->alpe_z2[pp] = alp_z * pow(pp / (zPML2 - 1.0), ma);
        g->kape_z2[pp] = 1.0 + (kap_z - 1) * pow((zPML2 - 1.0 - pp) / (zPML2 - 1.0), m);
        g->be_z2[pp] = exp(-(g->sige_z2[pp] / g->kape_z2[pp] + g->alpe_z2[pp]) / eps0 * dt);
        if ((g->sige_z2[pp] == 0.0) && (g->alpe_z2[pp] == 0.0) && (pp + 1 == zPML2))
            g->ce_z2[pp] = 0.0;
        else
            g->ce_z2[pp] = g->sige_z2[pp] * (g->be_z2[pp] - 1.0) / (g->sige_z2[pp] + g->kape_z2[pp] * g->alpe_z2[pp]) / g->kape_z2[pp];
    }
    for (pp = 0; pp < zPML2 - 1; pp++) {
        g->sigh_z2[pp] = sig_z * pow((zPML2 - 1.5 - pp) / (zPML2 - 1.0), m);
        g->alph_z2[pp] = alp_z * pow((pp + 0.5) / (zPML2 - 1.0), ma);
        g->kaph_z2[pp] = 1.0 + (kap_z - 1.0) * pow((zPML2 - 1.5 - pp) / (zPML2 - 1.0), m);
        g->bh_z2[pp] = exp(-(g->sigh_z2[pp] / g->kaph_z2[pp] + g->alph_z2[pp]) / eps0 * dt);
        g->ch_z2[pp] = g->sigh_z2[pp] * (g->bh_z2[pp] - 1.0) / (g->sigh_z2[pp] + g->kaph_z2[pp] * g->alph_z2[pp]) / g->kaph_z2[pp];
    }
    ii = xPML2 - 2;
    for (mm = 0; mm < SizeX - 1; mm++) {
        if (mm < xPML1 - 1) {
            g->den_hx[mm] = 1.0 / g->kaph_x1[mm] / dx;
        }
        else if (mm >= SizeX - xPML2) {
            g->den_hx[mm] = 1.0 / g->kaph_x2[ii] / dx;
            ii -= 1;
        }
        else {
            g->den_hx[mm] = 1.0 / dx;
        }
    }
    jj = yPML2 - 2;
    for (nn = 0; nn < SizeY - 1; nn++) {
        if (nn < yPML1 - 1) {
            g->den_hy[nn] = 1.0 / g->kaph_y1[nn] / dy;
        }
        else if (nn >= SizeY - yPML2) {
            g->den_hy[nn] = 1.0 / g->kaph_y2[jj] / dy;
            jj -= 1;
        }
        else {
            g->den_hy[nn] = 1.0 / dy;
        }
    }
    kk = zPML2 - 2;
    for (pp = 1; pp < SizeZ - 1; pp++) {
        if (pp < zPML1) {
            g->den_hz[pp] = 1.0 / g->kaph_z1[pp - 1] / dz;
        }
        else if (pp >= SizeZ - zPML2) {
            g->den_hz[pp] = 1.0 / g->kaph_z2[kk] / dz;
            kk -= 1;
        }
        else {
            g->den_hz[pp] = 1.0 / dz;
        }
    }
    ii = xPML2 - 1;
    for (mm = 0; mm < SizeX - 1; mm++) {
        if (mm < xPML1) {
            g->den_ex[mm] = 1.0 / g->kape_x1[mm] / dx;
        }
        else if (mm >= SizeX - xPML2) {
            g->den_ex[mm] = 1.0 / g->kape_x2[ii] / dx;
            ii -= 1;
        }
        else {
            g->den_ex[mm] = 1.0 / dx;
        }
    }
    jj = yPML2 - 1;
    for (nn = 0; nn < SizeY - 1; nn++) {
        if (nn < yPML1) {
            g->den_ey[nn] = 1.0 / g->kape_y1[nn] / dy;
        }
        else if (nn >= SizeY - yPML2) {
            g->den_ey[nn] = 1.0 / g->kape_y2[jj] / dy;
            jj -= 1;
        }
        else {
            g->den_ey[nn] = 1.0 / dy;
        }
    }
    kk = zPML2 - 1;
    for (pp = 0; pp < SizeZ - 1; pp++) {
        if (pp < zPML1) {
            g->den_ez[pp] = 1.0 / g->kape_z1[pp] / dz;
        }
        else if (pp >= SizeZ - zPML2 - 1) {
            g->den_ez[pp] = 1.0 / g->kape_z2[kk] / dz;
            kk -= 1;
        }
        else {
            g->den_ez[pp] = 1.0 / dz;
        }
    }
    /*if (!isAid) {
        

        
    }
    else {
        for (mm = 0; mm < SizeX - 1; mm++) {
            g->den_hx[mm] = 1.0 / dx;
        }
        for (nn = 0; nn < SizeY - 1; nn++) {
            g->den_hy[nn] = 1.0 / dy;
        }
        for (pp = 1; pp < SizeZ - 1; pp++) {
            g->den_hz[pp] = 1.0 / dz;
        }
        for (mm = 0; mm < SizeX - 1; mm++) {
            g->den_ex[mm] = 1.0 / dx;
        }
        for (nn = 0; nn < SizeY - 1; nn++) {
            g->den_ey[nn] = 1.0 / dy;
        }
        for (pp = 0; pp < SizeZ - 1; pp++) {
            g->den_ez[pp] = 1.0 / dz;
        }
    }*/
    
    
    
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

    if (isSpherePresent) {
        r2 = radius * radius;
        for (mm = 2; mm < SizeX - 2; mm++) {
            m2 = (mm + 0.5 - m_c) * (mm + 0.5 - m_c);
            for (nn = 2; nn < SizeY - 2; nn++) {
                n2 = (nn + 0.5 - n_c) * (nn + 0.5 - n_c);
                for (pp = 2; pp < SizeZ - 2; pp++) {
                    p2 = (pp + 0.5 - p_c) * (pp + 0.5 - p_c);

                    if (m2 + n2 + p2 < r2) {
                        Cexe(mm, nn, pp) = 0.0;
                        Cexe(mm, nn + 1, pp) = 0.0;
                        Cexe(mm, nn, pp + 1) = 0.0;
                        Cexe(mm, nn + 1, pp + 1) = 0.0;
                        Cexh(mm, nn, pp) = 0.0;
                        Cexh(mm, nn + 1, pp) = 0.0;
                        Cexh(mm, nn, pp + 1) = 0.0;
                        Cexh(mm, nn + 1, pp + 1) = 0.0;

                        Ceye(mm, nn, pp) = 0.0;
                        Ceye(mm + 1, nn, pp) = 0.0;
                        Ceye(mm, nn, pp + 1) = 0.0;
                        Ceye(mm + 1, nn, pp + 1) = 0.0;
                        Ceyh(mm, nn, pp) = 0.0;
                        Ceyh(mm + 1, nn, pp) = 0.0;
                        Ceyh(mm, nn, pp + 1) = 0.0;
                        Ceyh(mm + 1, nn, pp + 1) = 0.0;

                        Ceze(mm, nn, pp) = 0.0;
                        Ceze(mm + 1, nn, pp) = 0.0;
                        Ceze(mm, nn + 1, pp) = 0.0;
                        Ceze(mm + 1, nn + 1, pp) = 0.0;
                        Cezh(mm, nn, pp) = 0.0;
                        Cezh(mm + 1, nn, pp) = 0.0;
                        Cezh(mm, nn + 1, pp) = 0.0;
                        Cezh(mm + 1, nn + 1, pp) = 0.0;
                    }
                }
            }
        }
    }
    return;
}
