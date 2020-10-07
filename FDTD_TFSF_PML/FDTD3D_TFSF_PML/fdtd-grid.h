#ifndef _FDTD_GRID_H
#define _FDTD_GRID_H

enum GRIDTYPE { oneDGrid, teZGrid, tmZGrid, threeGrid };

struct Grid {
	double* hx, * chxh, * chxe;
	double* hy, * chyh, * chye;
	double* hz, * chzh, * chze;
	double* ex, * cexe, * cexh;
	double* ey, * ceye, * ceyh;
	double* ez, * ceze, * cezh;
	double* sig, * eps;
	double* p_Exy1, * p_Exy2;
	double* p_Exz1, * p_Exz2;
	double* p_Eyz1, * p_Eyz2;
	double* p_Eyx1, * p_Eyx2;
	double* p_Ezx1, * p_Ezx2;
	double* p_Ezy1, * p_Ezy2;
	double* p_Hxy1, * p_Hxy2;
	double* p_Hxz1, * p_Hxz2;
	double* p_Hyz1, * p_Hyz2;
	double* p_Hyx1, * p_Hyx2;
	double* p_Hzx1, * p_Hzx2;
	double* p_Hzy1, * p_Hzy2;
	double* be_x1, * ce_x1, * alpe_x1, * sige_x1, * kape_x1;
	double* be_x2, * ce_x2, * alpe_x2, * sige_x2, * kape_x2;
	double* be_y1, * ce_y1, * alpe_y1, * sige_y1, * kape_y1;
	double* be_y2, * ce_y2, * alpe_y2, * sige_y2, * kape_y2;
	double* be_z1, * ce_z1, * alpe_z1, * sige_z1, * kape_z1;
	double* be_z2, * ce_z2, * alpe_z2, * sige_z2, * kape_z2;
	double* bh_x1, * ch_x1, * alph_x1, * sigh_x1, * kaph_x1;
	double* bh_x2, * ch_x2, * alph_x2, * sigh_x2, * kaph_x2;
	double* bh_y1, * ch_y1, * alph_y1, * sigh_y1, * kaph_y1;
	double* bh_y2, * ch_y2, * alph_y2, * sigh_y2, * kaph_y2;
	double* bh_z1, * ch_z1, * alph_z1, * sigh_z1, * kaph_z1;
	double* bh_z2, * ch_z2, * alph_z2, * sigh_z2, * kaph_z2;
	double* den_ex, * den_hx, * den_ey, * den_hy, * den_ez, * den_hz;
	int sizeX, sizeY, sizeZ;
	int time, maxTime;
	int type;
	double cdtds;
};

typedef struct Grid Grid;

#endif
