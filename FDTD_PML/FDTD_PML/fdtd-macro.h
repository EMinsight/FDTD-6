#ifndef _FDTD_MACRO_H
#define _FDTD_MACRO_H

#include "fdtd-grid.h"

#define HxG(G, M, N, P)     G->hx[((M) * (SizeYG(G) - 1) + N) * (SizeZG(G) - 1) + P]
#define ChxhG(G, M, N, P)   G->chxh[((M) * (SizeYG(G)) + N) * (SizeZG(G)) + P]
#define ChxeG(G, M, N, P)   G->chxe[((M) * (SizeYG(G)) + N) * (SizeZG(G)) + P]

#define HyG(G, M, N, P)     G->hy[((M) * SizeYG(G) + N) * (SizeZG(G) - 1) + P]
#define ChyhG(G, M, N, P)   G->chyh[((M) * SizeYG(G) + N) * (SizeZG(G)) + P]
#define ChyeG(G, M, N, P)   G->chye[((M) * SizeYG(G) + N) * (SizeZG(G)) + P]

#define HzG(G, M, N, P)     G->hz[((M) * (SizeYG(G) - 1) + N) * SizeZG(G) + P]
#define ChzhG(G, M, N, P)   G->chzh[((M) * (SizeYG(G)) + N) * SizeZG(G) + P]
#define ChzeG(G, M, N, P)   G->chze[((M) * (SizeYG(G)) + N) * SizeZG(G) + P]

#define ExG(G, M, N, P)     G->ex[((M) * SizeYG(G) + N) * SizeZG(G) + P]
#define CexeG(G, M, N, P)   G->cexe[((M) * SizeYG(G) + N) * SizeZG(G) + P]
#define CexhG(G, M, N, P)   G->cexh[((M) * SizeYG(G) + N) * SizeZG(G) + P]

#define EyG(G, M, N, P)     G->ey[((M) * (SizeYG(G) - 1) + N) * SizeZG(G) + P]
#define CeyeG(G, M, N, P)   G->ceye[((M) * (SizeYG(G)) + N) * SizeZG(G) + P]
#define CeyhG(G, M, N, P)   G->ceyh[((M) * (SizeYG(G)) + N) * SizeZG(G) + P]

#define EzG(G, M, N, P)     G->ez[((M) * SizeYG(G) + N) * (SizeZG(G) - 1) + P]
#define CezeG(G, M, N, P)   G->ceze[((M) * SizeYG(G) + N) * (SizeZG(G)) + P]
#define CezhG(G, M, N, P)   G->cezh[((M) * SizeYG(G) + N) * (SizeZG(G)) + P]

#define SigG(G, M, N, P)    G->sig[((M) * SizeYG(G) + N) * (SizeZG(G)) + P]
#define EpsG(G, M, N, P)    G->eps[((M) * SizeYG(G) + N) * (SizeZG(G)) + P]

#define Exy1G(G, M, N, P)   G->p_Exy1[((M) * (yPML1) + N) * (SizeZG(G)) + P]
#define Exy2G(G, M, N, P)   G->p_Exy2[((M) * (yPML2) + N) * (SizeZG(G)) + P]
#define Exz1G(G, M, N, P)   G->p_Exz1[((M) * (SizeYG(G)) + N) * (zPML1) + P]
#define Exz2G(G, M, N, P)	G->p_Exz2[((M) * (SizeYG(G)) + N) * (zPML2) + P]
#define Eyx1G(G, M, N, P)	G->p_Eyx1[((M) * (SizeYG(G) - 1) + N) * (SizeZG(G)) + P]
#define Eyx2G(G, M, N, P)	G->p_Eyx2[((M) * (SizeYG(G) - 1) + N) * (SizeZG(G)) + P]
#define Eyz1G(G, M, N, P)	G->p_Eyz1[((M) * (SizeYG(G) - 1) + N) * (zPML1) + P]
#define Eyz2G(G, M, N, P)	G->p_Eyz2[((M) * (SizeYG(G) - 1) + N) * (zPML2) + P]
#define Ezx1G(G, M, N, P)	G->p_Ezx1[((M) * (SizeYG(G)) + N) * (SizeZG(G) - 1) + P]
#define Ezx2G(G, M, N, P)	G->p_Ezx2[((M) * (SizeYG(G)) + N) * (SizeZG(G) - 1) + P]
#define Ezy1G(G, M, N, P)	G->p_Ezy1[((M) * (yPML1) + N) * (SizeZG(G) - 1) + P]
#define Ezy2G(G, M, N, P)	G->p_Ezy2[((M) * (yPML2) + N) * (SizeZG(G) - 1) + P]

#define Hxy1G(G, M, N, P)   G->p_Hxy1[((M) * (yPML1 - 1) + N) * (SizeZG(G) - 1) + P]
#define Hxy2G(G, M, N, P)	G->p_Hxy2[((M) * (yPML2 - 1) + N) * (SizeZG(G) - 1) + P]
#define Hxz1G(G, M, N, P)	G->p_Hxz1[((M) * (SizeYG(G) - 1) + N) * (zPML1 - 1) + P]
#define Hxz2G(G, M, N, P)	G->p_Hxz2[((M) * (SizeYG(G) - 1) + N) * (zPML2 - 1) + P]
#define Hyx1G(G, M, N, P)	G->p_Hyx1[((M) * (SizeYG(G)) + N) * (SizeZG(G) - 1) + P]
#define Hyx2G(G, M, N, P)	G->p_Hyx2[((M) * (SizeYG(G)) + N) * (SizeZG(G) - 1) + P]
#define Hyz1G(G, M, N, P)	G->p_Hyz1[((M) * (SizeYG(G)) + N) * (zPML1 - 1) + P]
#define Hyz2G(G, M, N, P)	G->p_Hyz2[((M) * (SizeYG(G)) + N) * (zPML2 - 1) + P]
#define Hzx1G(G, M, N, P)	G->p_Hzx1[((M) * (SizeYG(G) - 1) + N) * (SizeZG(G)) + P]
#define Hzx2G(G, M, N, P)	G->p_Hzx2[((M) * (SizeYG(G) - 1) + N) * (SizeZG(G)) + P]
#define Hzy1G(G, M, N, P)	G->p_Hzy1[((M) * (yPML1 - 1) + N) * (SizeZG(G)) + P]
#define Hzy2G(G, M, N, P)	G->p_Hzy2[((M) * (yPML2 - 1) + N) * (SizeZG(G)) + P]


#define SizeXG(G)     G->sizeX
#define SizeYG(G)     G->sizeY
#define SizeZG(G)     G->sizeZ
#define TimeG(G)      G->time
#define MaxTimeG(G)   G->maxTime
#define CdtdsG(G)     G->cdtds
#define TypeG(G)      G->type

#define Hx(M, N, P)     HxG(g, M, N, P)
#define Chxh(M, N, P)   ChxhG(g, M, N, P)
#define Chxe(M, N, P)   ChxeG(g, M, N, P)

#define Hy(M, N, P)     HyG(g, M, N, P)
#define Chyh(M, N, P)   ChyhG(g, M, N, P)
#define Chye(M, N, P)   ChyeG(g, M, N, P)

#define Hz(M, N, P)     HzG(g, M, N, P)
#define Chzh(M, N, P)   ChzhG(g, M, N, P)
#define Chze(M, N, P)   ChzeG(g, M, N, P)


#define Ex(M, N, P)     ExG(g, M, N, P)
#define Cexe(M, N, P)   CexeG(g, M, N, P)
#define Cexh(M, N, P)   CexhG(g, M, N, P)

#define Ey(M, N, P)     EyG(g, M, N, P)
#define Ceye(M, N, P)   CeyeG(g, M, N, P)
#define Ceyh(M, N, P)   CeyhG(g, M, N, P)

#define Ez(M, N, P)     EzG(g, M, N, P)
#define Ceze(M, N, P)   CezeG(g, M, N, P)
#define Cezh(M, N, P)   CezhG(g, M, N, P)

#define Sig(M, N, P)    SigG(g, M, N, P)
#define Eps(M, N, P)    EpsG(g, M, N, P)

#define Exy1(M, N, P)   Exy1G(g, M, N, P)
#define Exy2(M, N, P)   Exy2G(g, M, N, P)
#define Exz1(M, N, P)   Exz1G(g, M, N, P)
#define Exz2(M, N, P)   Exz2G(g, M, N, P)
#define Eyx1(M, N, P)   Eyx1G(g, M, N, P)
#define Eyx2(M, N, P)   Eyx2G(g, M, N, P)
#define Eyz1(M, N, P)   Eyz1G(g, M, N, P)
#define Eyz2(M, N, P)   Eyz2G(g, M, N, P)
#define Ezx1(M, N, P)   Ezx1G(g, M, N, P)
#define Ezx2(M, N, P)   Ezx2G(g, M, N, P)
#define Ezy1(M, N, P)   Ezy1G(g, M, N, P)
#define Ezy2(M, N, P)   Ezy2G(g, M, N, P)

#define Hxy1(M, N, P)   Hxy1G(g, M, N, P)
#define Hxy2(M, N, P)   Hxy2G(g, M, N, P)
#define Hxz1(M, N, P)   Hxz1G(g, M, N, P)
#define Hxz2(M, N, P)   Hxz2G(g, M, N, P)
#define Hyx1(M, N, P)   Hyx1G(g, M, N, P)
#define Hyx2(M, N, P)   Hyx2G(g, M, N, P)
#define Hyz1(M, N, P)   Hyz1G(g, M, N, P)
#define Hyz2(M, N, P)   Hyz2G(g, M, N, P)
#define Hzx1(M, N, P)   Hzx1G(g, M, N, P)
#define Hzx2(M, N, P)   Hzx2G(g, M, N, P)
#define Hzy1(M, N, P)   Hzy1G(g, M, N, P)
#define Hzy2(M, N, P)   Hzy2G(g, M, N, P)

#define SizeX     SizeXG(g)
#define SizeY     SizeYG(g)
#define SizeZ     SizeZG(g)
#define Time      TimeG(g)
#define MaxTime   MaxTimeG(g)
#define Cdtds     CdtdsG(g)
#define Type      TypeG(g)


#endif
