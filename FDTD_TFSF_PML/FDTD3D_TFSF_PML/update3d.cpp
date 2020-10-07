
#include "fdtd-macro.h"
#include <stdio.h>

extern int xPML1;
extern int xPML2;
extern int yPML1;
extern int yPML2;
extern int zPML1;
extern int zPML2;

extern double dx;
extern double dy;
extern double dz;
extern double dt;

void updateH(Grid* g, bool isAid)
{
    int mm, nn, pp;
    int ii, jj, kk;

    if (Type == oneDGrid) {
        for (mm = 0; mm < SizeX - 1; mm++) {
            Hy1(mm) = Chyh1(mm) * Hy1(mm) + Chye1(mm) * (Ez1(mm + 1) - Ez1(mm));
        }
    }
    
    if (Type == tmZGrid) {
        for (mm = 0; mm < SizeX; mm++) {
            for (nn = 0; nn < SizeY - 1; nn++) {
                Hx2(mm, nn) = Chxh2(mm, nn) * Hx2(mm, nn) - Chxe2(mm, nn) * (Ez2(mm, nn + 1) - Ez2(mm, nn));
            }
        }
        for (mm = 0; mm < SizeX - 1; mm++) {
            for (nn = 0; nn < SizeY; nn++) {
                Hy2(mm, nn) = Chyh2(mm, nn) * Hy2(mm, nn) + Chye2(mm, nn) * (Ez2(mm + 1, nn) - Ez2(mm, nn));
            }
        }
    }
    if (Type == teZGrid) {
        for (mm = 0; mm < SizeX - 1; mm++) {
            for (nn = 0; nn < SizeY - 1; nn++) {
                Hz2(mm, nn) = Chzh2(mm, nn) * Hz2(mm, nn) - Chze2(mm, nn) * ((Ey2(mm + 1, nn) - Ey2(mm, nn)) - (Ex2(mm, nn + 1) - Ex2(mm, nn)));
            }
        }
    }
    if (Type == threeGrid) {
        for (mm = 0; mm < SizeX - 1; mm++) {
            for (nn = 0; nn < SizeY - 1; nn++) {
                for (pp = 0; pp < SizeZ - 1; pp++) {
                    Hx(mm, nn, pp) = Chxh(mm, nn, pp) * Hx(mm, nn, pp) + Chxe(mm, nn, pp) * ((Ey(mm, nn, pp + 1) - Ey(mm, nn, pp)) * g->den_hz[pp] - (Ez(mm, nn + 1, pp) - Ez(mm, nn, pp)) * g->den_hy[nn]);
                }
            }
        }
        
        for (mm = 0; mm < SizeX - 1; mm++) {
            for (nn = 0; nn < SizeY - 1; nn++) {
                for (pp = 0; pp < SizeZ - 1; pp++) {
                    Hy(mm, nn, pp) = Chyh(mm, nn, pp) * Hy(mm, nn, pp) + Chye(mm, nn, pp) * ((Ez(mm + 1, nn, pp) - Ez(mm, nn, pp)) * g->den_hx[mm] - (Ex(mm, nn, pp + 1) - Ex(mm, nn, pp)) * g->den_hz[pp]);
                }
            }
        }

        for (mm = 0; mm < SizeX - 1; mm++) {
            for (nn = 0; nn < SizeY - 1; nn++) {
                for (pp = 0; pp < SizeZ - 1; pp++) {

                    Hz(mm, nn, pp) = Chzh(mm, nn, pp) * Hz(mm, nn, pp) + Chze(mm, nn, pp) * ((Ex(mm, nn + 1, pp) - Ex(mm, nn, pp)) * g->den_hy[nn] - (Ey(mm + 1, nn, pp) - Ey(mm, nn, pp)) * g->den_hx[mm]);
                }
            }
        }
        if (!isAid) {
            for (pp = 0; pp < SizeZ - 1; pp++) {
                for (mm = 0; mm < SizeX - 1; mm++) {
                    for (nn = 0; nn < yPML1 - 1; nn++) {
                        Hxy1(mm, nn, pp) = g->bh_y1[nn] * Hxy1(mm, nn, pp) - g->ch_y1[nn] * (Ez(mm, nn + 1, pp) - Ez(mm, nn, pp)) / dy;
                        Hx(mm, nn, pp) = Hx(mm, nn, pp) + Chxe(mm, nn, pp) * Hxy1(mm, nn, pp);
                    }
                    jj = yPML2 - 2;
                    for (nn = SizeY - yPML2; nn < SizeY - 1; nn++) {
                        Hxy2(mm, jj, pp) = g->bh_y2[jj] * Hxy2(mm, jj, pp) - g->ch_y2[jj] * (Ez(mm, nn + 1, pp) - Ez(mm, nn, pp)) / dy;
                        Hx(mm, nn, pp) = Hx(mm, nn, pp) + Chxe(mm, nn, pp) * Hxy2(mm, jj, pp);
                        jj -= 1;
                    }
                }
            }
            for (mm = 0; mm < SizeX - 1; mm++) {
                for (nn = 0; nn < SizeY - 1; nn++) {
                    for (pp = 0; pp < zPML1 - 1; pp++) {
                        Hxz1(mm, nn, pp + 1) = g->bh_z1[pp + 1] * Hxz1(mm, nn, pp + 1) + g->ch_z1[pp + 1] * (Ey(mm, nn, pp + 1) - Ey(mm, nn, pp)) / dz;
                        Hx(mm, nn, pp) = Hx(mm, nn, pp) + Chxe(mm, nn, pp) * Hxz1(mm, nn, pp + 1);
                    }
                    kk = zPML2 - 2;
                    for (pp = SizeZ - zPML2; pp < SizeZ - 1; pp++) {
                        Hxz2(mm, nn, kk) = g->bh_z2[kk] * Hxz2(mm, nn, kk) + g->ch_z2[kk] * (Ey(mm, nn, pp) - Ey(mm, nn, pp - 1)) / dz;
                        Hx(mm, nn, pp) = Hx(mm, nn, pp) + Chxe(mm, nn, pp) * Hxz2(mm, nn, kk);
                        kk -= 1;
                    }
                }
            }
            for (pp = 0; pp < SizeZ - 1; pp++) {
                for (nn = 0; nn < SizeY - 1; nn++) {
                    for (mm = 0; mm < xPML1 - 1; mm++) {
                        Hyx1(mm, nn, pp) = g->bh_x1[mm] * Hyx1(mm, nn, pp) + g->ch_x1[mm] * (Ez(mm + 1, nn, pp) - Ez(mm, nn, pp)) / dx;
                        Hy(mm, nn, pp) = Hy(mm, nn, pp) + Chye(mm, nn, pp) * Hyx1(mm, nn, pp);
                    }
                    ii = xPML2 - 2;
                    for (mm = SizeX - xPML2; mm < SizeX - 1; mm++) {
                        Hyx2(ii, nn, pp) = g->bh_x2[ii] * Hyx2(ii, nn, pp) + g->ch_x2[ii] * (Ez(mm + 1, nn, pp) - Ez(mm, nn, pp)) / dx;
                        Hy(mm, nn, pp) = Hy(mm, nn, pp) + Chye(mm, nn, pp) * Hyx2(ii, nn, pp);
                        ii -= 1;
                    }
                }
            }
            for (mm = 0; mm < SizeX - 1; mm++) {
                for (nn = 0; nn < SizeY; nn++) {
                    for (pp = 0; pp < zPML1 - 1; pp++) {
                        Hyz1(mm, nn, pp - 1) = g->bh_z1[pp - 1] * Hyz1(mm, nn, pp - 1) - g->ch_z1[pp - 1] * (Ex(mm, nn, pp + 1) - Ex(mm, nn, pp)) / dz;
                        Hy(mm, nn, pp) = Hy(mm, nn, pp) + Chye(mm, nn, pp) * Hyz1(mm, nn, pp - 1);
                    }
                    kk = zPML2 - 2;
                    for (pp = SizeZ - zPML2; pp < SizeZ - 1; pp++) {
                        Hyz2(mm, nn, kk) = g->bh_z2[kk] * Hyz2(mm, nn, kk) - g->ch_z2[kk] * (Ex(mm, nn, pp + 1) - Ex(mm, nn, pp)) / dz;
                        Hy(mm, nn, pp) = Hy(mm, nn, pp) + Chye(mm, nn, pp) * Hyz2(mm, nn, kk);
                        kk -= 1;
                    }
                }
            }
            for (pp = 0; pp < SizeZ - 1; pp++) {
                for (nn = 0; nn < SizeY - 1; nn++) {
                    for (mm = 0; mm < xPML1 - 1; mm++) {
                        Hzx1(mm, nn, pp) = g->bh_x1[mm] * Hzx1(mm, nn, pp) - g->ch_x1[mm] * (Ey(mm + 1, nn, pp) - Ey(mm, nn, pp)) / dx;
                        Hz(mm, nn, pp) = Hz(mm, nn, pp) + Chze(mm, nn, pp) * Hzx1(mm, nn, pp);
                    }
                    ii = xPML2 - 2;
                    for (mm = SizeX - xPML2; mm < SizeX - 1; mm++) {
                        Hzx2(ii, nn, pp) = g->bh_x2[ii] * Hzx2(ii, nn, pp) - g->ch_x2[ii] * (Ey(mm + 1, nn, pp) - Ey(mm, nn, pp)) / dx;
                        Hz(mm, nn, pp) = Hz(mm, nn, pp) + Chze(mm, nn, pp) * Hzx2(ii, nn, pp);
                        ii -= 1;
                    }
                }
                for (mm = 0; mm < SizeX - 1; mm++) {
                    for (nn = 0; nn < yPML1 - 1; nn++) {
                        Hzy1(mm, nn, pp) = g->bh_y1[nn] * Hzy1(mm, nn, pp) + g->ch_y1[nn] * (Ex(mm, nn + 1, pp) - Ex(mm, nn, pp)) / dy;
                        Hz(mm, nn, pp) = Hz(mm, nn, pp) + Chze(mm, nn, pp) * Hzy1(mm, nn, pp);
                    }
                    jj = yPML2 - 2;
                    for (nn = SizeY - yPML2; nn < SizeY - 1; nn++) {
                        Hzy2(mm, jj, pp) = g->bh_y2[jj] * Hzy2(mm, jj, pp) + g->ch_y2[jj] * (Ex(mm, nn + 1, pp) - Ex(mm, nn, pp)) / dy;
                        Hz(mm, nn, pp) = Hz(mm, nn, pp) + Chze(mm, nn, pp) * Hzy2(mm, jj, pp);
                        jj -= 1;
                    }
                }
            }
        }
    }

    return;
}


void updateE(Grid* g, bool isAid)
{
    int mm, nn, pp;
    int ii, jj, kk;

    if (Type == oneDGrid) {
        for (mm = 1; mm < SizeX - 1; mm++) {
            Ez1(mm) = Ceze1(mm) * Ez1(mm) + Cezh1(mm) * (Hy1(mm) - Hy1(mm - 1));
        }
    }
    
    if (Type == tmZGrid) {
        for (mm = 1; mm < SizeX - 1; mm++) {
            for (nn = 1; nn < SizeY - 1; nn++) {
                Ez2(mm, nn) = Ceze2(mm, nn) * Ez2(mm, nn) + Cezh2(mm, nn) * ((Hy2(mm, nn) - Hy2(mm - 1, nn)) - (Hx2(mm, nn) - Hx2(mm, nn - 1)));
            }
        }
    }
    if (Type == teZGrid) {
        for (mm = 1; mm < SizeX - 1; mm++) {
            for (nn = 1; nn < SizeY - 1; nn++) {
                Ex2(mm, nn) = Cexe2(mm, nn) * Ex2(mm, nn) + Cexh2(mm, nn) * (Hz2(mm, nn) - Hz2(mm, nn - 1));
            }
        }
        for (mm = 1; mm < SizeX - 1; mm++) {
            for (nn = 1; nn < SizeY - 1; nn++) {
                Ey2(mm, nn) = Ceye2(mm, nn) * Ey2(mm, nn) - Ceyh2(mm, nn) * (Hz2(mm, nn) - Hz2(mm - 1, nn));
            }
        }
    }
    if (Type == threeGrid) {
        for (mm = 0; mm < SizeX - 1; mm++) {
            for (nn = 1; nn < SizeY - 1; nn++) {
                for (pp = 1; pp < SizeZ - 1; pp++) {

                    Ex(mm, nn, pp) = Cexe(mm, nn, pp) * Ex(mm, nn, pp) + Cexh(mm, nn, pp) * ((Hz(mm, nn, pp) - Hz(mm, nn - 1, pp)) * g->den_ey[nn] - (Hy(mm, nn, pp) - Hy(mm, nn, pp - 1)) * g->den_ez[pp]);

                }
            }
        }
        
        for (mm = 1; mm < SizeX - 1; mm++) {
            for (nn = 0; nn < SizeY - 1; nn++) {
                for (pp = 1; pp < SizeZ - 1; pp++) {
                    Ey(mm, nn, pp) = Ceye(mm, nn, pp) * Ey(mm, nn, pp) + Ceyh(mm, nn, pp) * ((Hx(mm, nn, pp) - Hx(mm, nn, pp - 1)) * g->den_ez[pp] - (Hz(mm, nn, pp) - Hz(mm - 1, nn, pp)) * g->den_ex[mm]);
                }
            }
        }
        
        for (mm = 1; mm < SizeX - 1; mm++) {
            for (nn = 1; nn < SizeY - 1; nn++) {
                for (pp = 0; pp < SizeZ - 1; pp++) {
                    Ez(mm, nn, pp) = Ceze(mm, nn, pp) * Ez(mm, nn, pp) + Cezh(mm, nn, pp) * ((Hy(mm, nn, pp) - Hy(mm - 1, nn, pp)) * g->den_ex[mm] - (Hx(mm, nn, pp) - Hx(mm, nn - 1, pp)) * g->den_ey[nn]);
                }
            }
        }
        if (!isAid) {
            for (pp = 1; pp < SizeZ - 1; pp++) {
                for (mm = 0; mm < SizeX - 1; mm++) {
                    for (nn = 1; nn < yPML1; nn++) {
                        Exy1(mm, nn, pp) = g->be_y1[nn] * Exy1(mm, nn, pp) + g->ce_y1[nn] * (Hz(mm, nn, pp) - Hz(mm, nn - 1, pp)) / dy;
                        Ex(mm, nn, pp) = Ex(mm, nn, pp) + Cexh(mm, nn, pp) * Exy1(mm, nn, pp);
                    }
                    jj = yPML2 - 1;
                    for (nn = SizeY - yPML2; nn < SizeY - 1; nn++) {
                        Exy2(mm, jj, pp) = g->be_y2[jj] * Exy2(mm, jj, pp) + g->ce_y2[jj] * (Hz(mm, nn, pp) - Hz(mm, nn - 1, pp)) / dy;
                        Ex(mm, nn, pp) = Ex(mm, nn, pp) + Cexh(mm, nn, pp) * Exy2(mm, jj, pp);
                        jj -= 1;
                    }
                }
            }
            for (mm = 0; mm < SizeX - 1; mm++) {
                for (nn = 1; nn < SizeY - 1; nn++) {
                    for (pp = 1; pp < zPML1; pp++) {
                        Exz1(mm, nn, pp) = g->be_z1[pp] * Exz1(mm, nn, pp) - g->ce_z1[pp] * (Hy(mm, nn, pp) - Hy(mm, nn, pp - 1)) / dz;
                        Ex(mm, nn, pp) = Ex(mm, nn, pp) + Cexh(mm, nn, pp) * Exz1(mm, nn, pp);
                    }
                    kk = zPML2 - 1;
                    for (pp = SizeZ - zPML2 - 1; pp < SizeZ - 1; pp++) {
                        Exz2(mm, nn, kk) = g->be_z2[kk] * Exz2(mm, nn, kk) - g->ce_z2[kk] * (Hy(mm, nn, pp) - Hy(mm, nn, pp - 1)) / dz;
                        Ex(mm, nn, pp) = Ex(mm, nn, pp) + Cexh(mm, nn, pp) * Exz2(mm, nn, kk);
                        kk -= 1;
                    }
                }
            }
            for (pp = 1; pp < SizeZ - 1; pp++) {
                for (nn = 0; nn < SizeY - 1; nn++) {
                    for (mm = 1; mm < xPML1; mm++) {
                        Eyx1(mm, nn, pp) = g->be_x1[mm] * Eyx1(mm, nn, pp) - g->ce_x1[mm] * (Hz(mm, nn, pp) - Hz(mm - 1, nn, pp)) / dx;
                        Ey(mm, nn, pp) = Ey(mm, nn, pp) + Ceyh(mm, nn, pp) * Eyx1(mm, nn, pp);
                    }
                    ii = xPML2 - 1;
                    for (mm = SizeX - xPML2; mm < SizeX - 1; mm++) {
                        Eyx2(ii, nn, pp) = g->be_x2[ii] * Eyx2(ii, nn, pp) - g->ce_x2[ii] * (Hz(mm, nn, pp) - Hz(mm - 1, nn, pp)) / dx;
                        Ey(mm, nn, pp) = Ey(mm, nn, pp) + Ceyh(mm, nn, pp) * Eyx2(ii, nn, pp);
                        ii -= 1;
                    }
                }
            }
            for (mm = 1; mm < SizeX - 1; mm++) {
                for (nn = 0; nn < SizeY - 1; nn++) {
                    for (pp = 1; pp < zPML1; pp++) {
                        Eyz1(mm, nn, pp) = g->be_z1[pp] * Eyz1(mm, nn, pp) + g->ce_z1[pp] * (Hx(mm, nn, pp) - Hx(mm, nn, pp - 1)) / dz;
                        Ey(mm, nn, pp) = Ey(mm, nn, pp) + Ceyh(mm, nn, pp) * Eyz1(mm, nn, pp);
                    }
                    kk = zPML2 - 1;
                    for (pp = SizeZ - zPML2 - 1; pp < SizeZ - 1; pp++) {
                        Eyz2(mm, nn, kk) = g->be_z2[kk] * Eyz2(mm, nn, kk) + g->ce_z2[kk] * (Hx(mm, nn, pp) - Hx(mm, nn, pp - 1)) / dz;
                        Ey(mm, nn, pp) = Ey(mm, nn, pp) + Ceyh(mm, nn, pp) * Eyz2(mm, nn, kk);
                        kk -= 1;
                    }
                }
            }
            for (pp = 0; pp < SizeZ - 1; pp++) {
                for (nn = 1; nn < SizeY - 1; nn++) {
                    for (mm = 1; mm < xPML1; mm++) {
                        Ezx1(mm, nn, pp) = g->be_x1[mm] * Ezx1(mm, nn, pp) + g->ce_x1[mm] * (Hy(mm, nn, pp) - Hy(mm - 1, nn, pp)) / dx;
                        Ez(mm, nn, pp) = Ez(mm, nn, pp) + Cezh(mm, nn, pp) * Ezx1(mm, nn, pp);
                    }
                    ii = xPML2 - 1;
                    for (mm = SizeX - xPML2; mm < SizeX - 1; mm++) {
                        Ezx2(ii, nn, pp) = g->be_x2[ii] * Ezx2(ii, nn, pp) + g->ce_x2[ii] * (Hy(mm, nn, pp) - Hy(mm - 1, nn, pp)) / dx;
                        Ez(mm, nn, pp) = Ez(mm, nn, pp) + Cezh(mm, nn, pp) * Ezx2(ii, nn, pp);
                        ii -= 1;
                    }
                }
                for (mm = 1; mm < SizeX - 1; mm++) {
                    for (nn = 1; nn < yPML1; nn++) {
                        Ezy1(mm, nn, pp) = g->be_y1[nn] * Ezy1(mm, nn, pp) - g->ce_y1[nn] * (Hx(mm, nn, pp) - Hx(mm, nn - 1, pp)) / dy;
                        Ez(mm, nn, pp) = Ez(mm, nn, pp) + Cezh(mm, nn, pp) * Ezy1(mm, nn, pp);
                    }
                    jj = yPML2 - 1;
                    for (nn = SizeY - yPML2; nn < SizeY - 1; nn++) {
                        Ezy2(mm, jj, pp) = g->be_y2[jj] * Ezy2(mm, jj, pp) - g->ce_y2[jj] * (Hx(mm, nn, pp) - Hx(mm, nn - 1, pp)) / dy;
                        Ez(mm, nn, pp) = Ez(mm, nn, pp) + Cezh(mm, nn, pp) * Ezy2(mm, jj, pp);
                        jj -= 1;
                    }
                }
            }
        }
    }

    return;
}
