/***********************************************************/
/*                FDTD3D with PML boundary                 */
/***********************************************************/
/*                                                         */
/*                  Author��Tianshaowen                    */
/*                    Date��2020/9/23                      */
/*               tianshaowenwen@buaa.edu.cn                */
/***********************************************************/
#include "fdtd-alloc.h"
#include "fdtd-macro.h"
#include "fdtd-proto.h"
#include <math.h>

//�����С
int sizeX = 100;
int sizeY = 100;
int sizeZ = 100;

//���񲽳�
double dx = 1.0e-3;
double dy = 1.0e-3;
double dz = 1.0e-3;

//ʱ�䲽��
double dt;

//������ʱ�䲽������ʼ��¼ʱ��
int maxTime = 300;

//�������
double mu0;
double eps0;
double eps;
double epsR = 1.0;

double c0 = 299792458;
double PI = 3.141592653589793238463;

int main()
{
	//��������
	dt = dx / (sqrt(3) * c0);
	mu0 = 4.0 * PI * 1e-7;
	eps0 = 1.0 / (c0 * c0 * mu0);
	eps = eps0 * epsR;

	Grid* g;
	ALLOC_1D(g, 1, Grid);
	gridInit(g);

	for (Time = 0; Time < MaxTime; Time++) {
		//���µ糡
		updateE(g);
		//��Ӳ�Դ
		for (int pp = 0; pp < sizeZ - 1; pp++) {
			Ez(50, 50, pp) += ezInc(Time * 1.0, 0.0);
		}
		//���´ų�
		updateH(g);
		//��¼����
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
