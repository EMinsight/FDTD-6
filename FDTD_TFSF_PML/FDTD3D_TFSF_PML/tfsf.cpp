
#include <string.h>
#include <cmath>
#include "fdtd-macro.h"
#include "fdtd-proto.h"
#include "fdtd-alloc.h"

extern int    firstX, lastX, firstY, lastY, firstZ, lastZ, sizeX, sizeY, sizeZ;
extern int    inDirection;
extern double ths;
extern double phs;
extern double er[3];
extern double PI;

static Grid* g1, * g2;

void tfsfInit(Grid* g)
{

    ALLOC_1D(g1, 1, Grid);
    memcpy(g1, g, sizeof(Grid));
    gridInit1d(g1);

    ALLOC_1D(g2, 1, Grid);
    memcpy(g2, g, sizeof(Grid));
    gridInit(g2, 0, true);

    return;
}

void tfsf(Grid* g)
{
    update_tfsf_H(g);

    updateH(g1, false);
    updateE(g1, false);
                
    updateH(g2, false);
    updateE(g2, false);

    Ez1G(g1, 0) += ezInc(TimeG(g1), 0.0);
    TimeG(g1)++;

    addSource(g2);
    TimeG(g2)++;

    update_tfsf_E(g);
    //snapshot3d(g2);
    
    return;
}
double ezInc(double time, double location)
{
    double ppw = 20.0;
    double arg;

    arg = PI * ((time - location) / ppw - 1.0);
    arg = arg * arg;
    double ans = (1.0 - 2.0 * arg) * exp(-arg);
    return ans;

}
void addSource(Grid* g) {
    int mm, nn, pp;
    int p;
    double r, w;

    if (inDirection == 1) {
        for (pp = 0; pp < SizeZ; pp++) {
            nn = 20;
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
            mm = 20;
            for (nn = 0; nn < SizeY; nn++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
        pp = 20;
        for (nn = 0; nn < SizeY; nn++) {
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
    }
    else if (inDirection == 2) {
        for (pp = 0; pp < SizeZ; pp++) {
            nn = 20;
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
            mm = SizeX - 22;
            for (nn = 0; nn < SizeY; nn++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
        pp = 20;
        for (nn = 0; nn < SizeY; nn++) {
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
    }
    else if (inDirection == 3) {
        for (pp = 0; pp < SizeZ; pp++) {
            nn = SizeY - 22;
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
            mm = SizeX - 22;
            for (nn = 0; nn < SizeY; nn++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
        pp = 20;
        for (nn = 0; nn < SizeY; nn++) {
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
    }
    else if (inDirection == 4) {
        for (pp = 0; pp < SizeZ; pp++) {
            nn = SizeY - 22;
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                  ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                  EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
           mm = 20;
           for (nn = 0; nn < SizeY; nn++) {
               r = getDist(mm, nn, pp, inDirection);
               p = (int)floor(r);
               w = r - p;
               ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
               if (nn < SizeY - 1)
                   EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
               if (pp < SizeZ - 1)
                   EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
           }
        }
        pp = 20;
        for (nn = 0; nn < SizeY; nn++) {
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
    }
    else if (inDirection == 5) {
        for (pp = 0; pp < SizeZ; pp++) {
            nn = 20;
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
            mm = 20;
            for (nn = 0; nn < SizeY; nn++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
        pp = SizeZ - 22;
        for (nn = 0; nn < SizeY; nn++) {
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
    }
    else if (inDirection == 6) {
        for (pp = 0; pp < SizeZ; pp++) {
            nn = 20;
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
            mm = SizeX - 22;
            for (nn = 0; nn < SizeY; nn++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
        pp = SizeZ - 22;
        for (nn = 0; nn < SizeY; nn++) {
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
    }
    else if (inDirection == 7) {
        for (pp = 0; pp < SizeZ; pp++) {
            nn = SizeY - 22;
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
            mm = SizeX - 22;
            for (nn = 0; nn < SizeY; nn++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
        pp = SizeZ - 22;
        for (nn = 0; nn < SizeY; nn++) {
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
    }
    else if (inDirection == 8) {
        for (pp = 0; pp < SizeZ; pp++) {
            nn = SizeY - 22;
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
            mm = 20;
            for (nn = 0; nn < SizeY; nn++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                if (pp < SizeZ - 1)
                    EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
        pp = SizeZ - 22;
        for (nn = 0; nn < SizeY; nn++) {
            for (mm = 0; mm < SizeX; mm++) {
                r = getDist(mm, nn, pp, inDirection);
                p = (int)floor(r);
                w = r - p;
                if (mm < SizeX - 1)
                    ExG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[0];
                if (nn < SizeY - 1)
                    EyG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[1];
                EzG(g, mm, nn, pp) = ((1 - w) * Ez1G(g1, p) + w * Ez1G(g1, p + 1)) * er[2];
            }
        }
    }
    
}
double getDist(int x, int y, int z, int inDirection) {
    double x1 = 0;
    double y1 = 0;
    double z1 = 0;
    double x2 = 0;
    double y2 = 0;
    double z2 = 0;
    if (inDirection == 1) {
        x1 = 0;
        y1 = 0;
        z1 = 0;
        x2 = sin(ths) * cos(phs);
        y2 = sin(ths) * sin(phs);
        z2 = cos(ths);
    }
    else if (inDirection == 2) {
        x1 = sizeX;
        y1 = 0;
        z1 = 0;
        x2 = sizeX + sin(ths) * cos(phs);
        y2 = sin(ths) * sin(phs);
        z2 = cos(ths);
    }
    else if (inDirection == 3) {
        x1 = sizeX;
        y1 = sizeY;
        z1 = 0;
        x2 = sizeX + sin(ths) * cos(phs);
        y2 = sizeY + sin(ths) * sin(phs);
        z2 = cos(ths);
    }
    else if (inDirection == 4) {
        x1 = 0;
        y1 = sizeY;
        z1 = 0;
        x2 = sin(ths) * cos(phs);
        y2 = sizeY + sin(ths) * sin(phs);
        z2 = cos(ths);
    }
    else if (inDirection == 5) {
        x1 = 0;
        y1 = 0;
        z1 = sizeZ;
        x2 = sin(ths) * cos(phs);
        y2 = sin(ths) * sin(phs);
        z2 = sizeZ + cos(ths);
    }
    else if (inDirection == 6) {
        x1 = sizeX;
        y1 = 0;
        z1 = sizeZ;
        x2 = sizeX + sin(ths) * cos(phs);
        y2 = sin(ths) * sin(phs);
        z2 = sizeZ + cos(ths);
    }
    else if (inDirection == 7) {
        x1 = sizeX;
        y1 = sizeY;
        z1 = sizeZ;
        x2 = sizeX + sin(ths) * cos(phs);
        y2 = sizeY + sin(ths) * sin(phs);
        z2 = sizeZ + cos(ths);
    }
    else if (inDirection == 8) {
        x1 = 0;
        y1 = sizeY;
        z1 = sizeZ;
        x2 = sin(ths) * cos(phs);
        y2 = sizeY + sin(ths) * sin(phs);
        z2 = sizeZ + cos(ths);
    }
    double ljx = sin(ths) * cos(phs);
    double ljy = sin(ths) * sin(phs);
    double ljz = cos(ths);
    double m = (x - x1) * ljx + (y - y1) * ljy + (z - z1) * ljz;
    double X = x1 + ljx * m;
    double Y = y1 + ljy * m;
    double Z = z1 + ljz * m;
    return sqrt((X - x1) * (X - x1) + (Y - y1) * (Y - y1) + (Z - z1) * (Z - z1));
}
void update_tfsf_H(Grid* g) 
{
    int mm, nn, pp;

    //在x平面上修正散射场
    mm = firstX;
    for (nn = firstY; nn <= lastY; nn++) {

        for (pp = firstZ; pp < lastZ; pp++) {
            Hy(mm - 1, nn, pp) -= Chye(mm - 1, nn, pp) * EzG(g2, mm, nn, pp) * g->den_hx[mm];

        }
    }
    for (nn = firstY; nn < lastY; nn++) {

        for (pp = firstZ; pp <= lastZ; pp++) {
            Hz(mm - 1, nn, pp) += Chze(mm - 1, nn, pp) * EyG(g2, mm, nn, pp) * g->den_hx[mm];

        }
    }
    mm = lastX;
    for (nn = firstY; nn <= lastY; nn++) {

        for (pp = firstZ; pp < lastZ; pp++) {
            Hy(mm, nn, pp) += Chye(mm, nn, pp) * EzG(g2, mm, nn, pp) * g->den_hx[mm];

        }
    }
    for (nn = firstY; nn < lastY; nn++) {

        for (pp = firstZ; pp <= lastZ; pp++) {
            Hz(mm, nn, pp) -= Chze(mm, nn, pp) * EyG(g2, mm, nn, pp) * g->den_hx[mm];

        }
    }
    //在y平面上修正散射场
    nn = firstY;
    for (mm = firstX; mm <= lastX; mm++) {
        for (pp = firstZ; pp < lastZ; pp++) {
            Hx(mm, nn - 1, pp) += Chxe(mm, nn - 1, pp) * EzG(g2, mm, nn, pp) * g->den_hy[nn];
        }
    }
    for (mm = firstX; mm < lastX; mm++) {
        for (pp = firstZ; pp <= lastZ; pp++) {
            Hz(mm, nn - 1, pp) -= Chze(mm, nn - 1, pp) * ExG(g2, mm, nn, pp) * g->den_hy[nn];
        }
    }
    nn = lastY;
    for (mm = firstX; mm <= lastX; mm++) {
        for (pp = firstZ; pp < lastZ; pp++) {
            Hx(mm, nn, pp) -= Chxe(mm, nn, pp) * EzG(g2, mm, nn, pp) * g->den_hy[nn];
        }
    }
    for (mm = firstX; mm < lastX; mm++) {
        for (pp = firstZ; pp <= lastZ; pp++) {
            Hz(mm, nn, pp) += Chze(mm, nn, pp) * ExG(g2, mm, nn, pp) * g->den_hy[nn];
        }
    }
    //在z平面上修正散射场
    pp = firstZ;
    for (mm = firstX; mm <= lastX; mm++) {
        for (nn = firstY; nn < lastY; nn++) {
            Hx(mm, nn, pp - 1) -= Chxe(mm, nn, pp - 1) * EyG(g2, mm, nn, pp) * g->den_hz[pp];
        }
    }
    for (mm = firstX; mm < lastX; mm++) {
        for (nn = firstY; nn <= lastY; nn++) {
            Hy(mm, nn, pp - 1) += Chye(mm, nn, pp - 1) * ExG(g2, mm, nn, pp) * g->den_hz[pp];
        }
    }
    pp = lastZ;
    for (mm = firstX; mm <= lastX; mm++) {
        for (nn = firstY; nn < lastY; nn++) {
            Hx(mm, nn, pp) += Chxe(mm, nn, pp) * EyG(g2, mm, nn, pp) * g->den_hz[pp];
        }
    }
    for (mm = firstX; mm < lastX; mm++) {
        for (nn = firstY; nn <= lastY; nn++) {
            Hy(mm, nn, pp) -= Chye(mm, nn, pp) * ExG(g2, mm, nn, pp) * g->den_hz[pp];
        }
    }
}
void update_tfsf_E(Grid* g) {
    int mm, nn, pp;
    //在x平面上修正总场
    mm = firstX;
    for (nn = firstY; nn <= lastY; nn++) {
        for (pp = firstZ; pp < lastZ; pp++) {
            Ez(mm, nn, pp) -= Cezh(mm, nn, pp) * HyG(g2, mm - 1, nn, pp) * g->den_ex[mm];

        }
    }
    for (nn = firstY; nn < lastY; nn++) {
        for (pp = firstZ; pp <= lastZ; pp++) {
            Ey(mm, nn, pp) += Ceyh(mm, nn, pp) * HzG(g2, mm - 1, nn, pp) * g->den_ex[mm];

        }
    }
    mm = lastX;
    for (nn = firstY; nn <= lastY; nn++) {
        for (pp = firstZ; pp < lastZ; pp++) {
            Ez(mm, nn, pp) += Cezh(mm, nn, pp) * HyG(g2, mm, nn, pp) * g->den_ex[mm];

        }
    }
    for (nn = firstY; nn < lastY; nn++) {
        for (pp = firstZ; pp <= lastZ; pp++) {
            Ey(mm, nn, pp) -= Ceyh(mm, nn, pp) * HzG(g2, mm, nn, pp) * g->den_ex[mm];

        }
    }
    //在y平面上修正总场
    nn = firstY;
    for (mm = firstX; mm <= lastX; mm++) {
        for (pp = firstZ; pp < lastZ; pp++) {
            Ez(mm, nn, pp) += Cezh(mm, nn, pp) * HxG(g2, mm, nn - 1, pp) * g->den_ey[nn];
        }
    }
    for (mm = firstX; mm < lastX; mm++) {
        for (pp = firstZ; pp <= lastZ; pp++) {
            Ex(mm, nn, pp) -= Cexh(mm, nn, pp) * HzG(g2, mm, nn - 1, pp) * g->den_ey[nn];
        }
    }
    nn = lastY;
    for (mm = firstX; mm <= lastX; mm++) {
        for (pp = firstZ; pp < lastZ; pp++) {
            Ez(mm, nn, pp) -= Cezh(mm, nn, pp) * HxG(g2, mm, nn, pp) * g->den_ey[nn];
        }
    }
    for (mm = firstX; mm < lastX; mm++) {
        for (pp = firstZ; pp <= lastZ; pp++) {
            Ex(mm, nn, pp) += Cexh(mm, nn, pp) * HzG(g2, mm, nn, pp) * g->den_ey[nn];
        }
    }
    //在z平面上修正总场
    pp = firstZ;
    for (mm = firstX; mm <= lastX; mm++) {
        for (nn = firstY; nn < lastY; nn++) {
            Ey(mm, nn, pp) -= Ceyh(mm, nn, pp) * HxG(g2, mm, nn, pp - 1) * g->den_ez[pp];
        }
    }
    for (mm = firstX; mm < lastX; mm++) {
        for (nn = firstY; nn <= lastY; nn++) {
            Ex(mm, nn, pp) += Cexh(mm, nn, pp) * HyG(g2, mm, nn, pp - 1) * g->den_ez[pp];
        }
    }
    pp = lastZ;
    for (mm = firstX; mm <= lastX; mm++) {
        for (nn = firstY; nn < lastY; nn++) {
            Ey(mm, nn, pp) += Ceyh(mm, nn, pp) * HxG(g2, mm, nn, pp) * g->den_ez[pp];
        }
    }
    for (mm = firstX; mm < lastX; mm++) {
        for (nn = firstY; nn <= lastY; nn++) {
            Ex(mm, nn, pp) -= Cexh(mm, nn, pp) * HyG(g2, mm, nn, pp) * g->den_ez[pp];
        }
    }
}