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

#define SizeX     SizeXG(g)
#define SizeY     SizeYG(g)
#define SizeZ     SizeZG(g)
#define Time      TimeG(g)
#define MaxTime   MaxTimeG(g)
#define Cdtds     CdtdsG(g)
#define Type      TypeG(g)


#endif
