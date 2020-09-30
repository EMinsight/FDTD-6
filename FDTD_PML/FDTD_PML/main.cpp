/***********************************************************/
/*                FDTD3D with PML boundary                 */
/***********************************************************/
/*                                                         */
/*                  Author：Tianshaowen                    */
/*                    Date：2020/9/30                      */
/*               tianshaowenwen@buaa.edu.cn                */
/***********************************************************/
#include "fdtd-alloc.h"
#include "fdtd-macro.h"
#include "fdtd-proto.h"
#include <math.h>

//网格大小
int sizeX = 100;
int sizeY = 100;
int sizeZ = 100;

//CPML边界网格大小
int xPML1 = 20;
int xPML2 = 20;
int yPML1 = 20;
int yPML2 = 20;
int zPML1 = 20;
int zPML2 = 20;

//网格步长
double dx = 1.0e-3;
double dy = 1.0e-3;
double dz = 1.0e-3;

//时间步长
double dt;

//总运行时间步长及开始记录时间
int maxTime = 300;

//网格参数
int m = 3;
int ma = 1;

double mu0;
double eps0;
double eps;
double epsR;
double sig_x;
double sig_y;
double sig_z;
double alp_x;
double alp_y;
double alp_z;
double kap_x;
double kap_y;
double kap_z;

double c0 = 299792458;
double PI = 3.141592653589793238463;

int main()
{
	//参数设置
	dt = dx / (sqrt(3) * c0);
	epsR = 1.0;
	mu0 = 4.0 * PI * 1e-7;
	eps0 = 1.0 / (c0 * c0 * mu0);
	eps = eps0 * epsR;
	sig_x = 0.75 * (0.8 * (m + 1.0) / sqrt(mu0 / eps0 * epsR) / dx);
	sig_y = 0.75 * (0.8 * (m + 1.0) / sqrt(mu0 / eps0 * epsR) / dy);
	sig_z = 0.75 * (0.8 * (m + 1.0) / sqrt(mu0 / eps0 * epsR) / dz);
	alp_x = 0.24;
	alp_y = 0.24;
	alp_z = 0.24;
	kap_x = 15.0;
	kap_y = 15.0;
	kap_z = 15.0;

	Grid* g;
	ALLOC_1D(g, 1, Grid);
	gridInit(g);

	for (Time = 0; Time < MaxTime; Time++) {
		//更新电场
		updateE(g);
		//添加波源
		for (int pp = 0; pp < sizeZ - 1; pp++) {
			Ez(50, 50, pp) += ezInc(Time * 1.0, 0.0);
		}
		//更新磁场
		updateH(g);
		//记录数据
		snapshot3d(g);
	}
}
double ezInc(double time, double location)
{
	double arg;
	double ppw = 20.0;

	arg = PI * ((time - location) / ppw - 1.0);
	arg = arg * arg;
	double ans = (1.0 - 2.0 * arg) * exp(-arg);
	return ans;

}
