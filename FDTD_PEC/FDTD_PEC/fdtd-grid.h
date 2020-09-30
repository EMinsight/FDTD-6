#ifndef _FDTD_GRID_H
#define _FDTD_GRID_H

//定义三维网格
struct Grid {
	double* hx, * chxh, * chxe;
	double* hy, * chyh, * chye;
	double* hz, * chzh, * chze;
	double* ex, * cexe, * cexh;
	double* ey, * ceye, * ceyh;
	double* ez, * ceze, * cezh;
	double* sig, * eps;
	int sizeX, sizeY, sizeZ;
	int time, maxTime;
};

typedef struct Grid Grid;

#endif
