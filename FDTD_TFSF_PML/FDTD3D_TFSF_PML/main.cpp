/***********************************************************/
/*           ����TFSF�������ⷽ������ƽ�沨+PML�߽�        */
/***********************************************************/
/*                                                         */
/*                  Author��Tianshaowen                    */
/*                    Date��2020/10/7                      */
/*               tianshaowenwen@buaa.edu.cn                */
/***********************************************************/
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "fdtd-alloc.h"
#include "fdtd-macro.h"
#include "fdtd-proto.h"

using namespace std;

//���������Ƿ���ý�����    1->��   0->��
int    isSpherePresent;
double radius = 15.0;
//���䲨����
double phi;
double thi;
int pol = 0;              //���䲨��������  0->thet����   1->phi����
//�����С
int sizeX = 120;
int sizeY = 120;
int sizeZ = 120;
//���񲽳�
double dx = 1.0e-3;
double dy = 1.0e-3;
double dz = 1.0e-3;
//ʱ�䲽��
double dt;
//CPML�߽������С
int xPML1 = 15;
int xPML2 = 15;
int yPML1 = 15;
int yPML2 = 15;
int zPML1 = 15;
int zPML2 = 15;
//������ʱ�䲽������ʼ��¼ʱ��
int maxTime = 300;
int startTime;
//TFSF�߽�
int firstX = 35;
int firstY = 35;
int firstZ = 35;
int lastX = 85;
int lastY = 85;
int lastZ = 85;
//��¼���䷽��   ��Ӧ����
int inDirection = 0;

double er[3];
double ths;
double phs;

//�������
int m = 3;
int ma = 1;
double mu0;
double eps0;
double eps;
double sig_x;
double sig_y;
double sig_z;
double alp_x;
double alp_y;
double alp_z;
double kap_x;
double kap_y;
double kap_z;
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
	sig_x = 0.75 * (0.8 * (m + 1.0) / sqrt(mu0 / eps0 * epsR) / dx);
	sig_y = 0.75 * (0.8 * (m + 1.0) / sqrt(mu0 / eps0 * epsR) / dy);
	sig_z = 0.75 * (0.8 * (m + 1.0) / sqrt(mu0 / eps0 * epsR) / dz);
	alp_x = 0.24;
	alp_y = 0.24;
	alp_z = 0.24;
	kap_x = 15.0;
	kap_y = 15.0;
	kap_z = 15.0;
	//����������
	std::fstream fp("./InputPara.txt");
	char temp[100];
	fp >> temp >> thi;
	fp >> temp >> phi;
	fp >> temp >> pol;
	fp >> temp >> startTime;
	fp >> temp >> isSpherePresent;
	fp.close();
	if (thi >= 0 && thi < 90) {
		if (phi >= 0 && phi < 90)
			inDirection = 1;
		else if (phi >= 90 && phi < 180)
			inDirection = 2;
		else if (phi >= 180 && phi < 270)
			inDirection = 3;
		else if(phi >= 270 && phi <= 360)
			inDirection = 4;
	}
	else if (thi >= 90 && thi <= 180) {
		if (phi >= 0 && phi < 90)
			inDirection = 5;
		else if (phi >= 90 && phi < 180)
			inDirection = 6;
		else if (phi >= 180 && phi < 270)
			inDirection = 7;
		else if (phi >= 270 && phi <= 360)
			inDirection = 8;
	}

	ths = thi * PI / 180.;
	phs = phi * PI / 180.;
	if (pol == 0) {
		er[0] = -cos(ths) * cos(phs);
		er[1] = -cos(ths) * sin(phs);
		er[2] = sin(ths);
	}
	else if (pol == 1) {
		er[0] = sin(phs);
		er[1] = -cos(phs);
		er[2] = 0.;
	}

	Grid* g;
	ALLOC_1D(g, 1, Grid);
	gridInit(g, isSpherePresent, false);
	tfsfInit(g);
	
	for (Time = 0; Time < MaxTime; Time++) {

		updateH(g, false);
		tfsf(g);
		updateE(g, false);
		
		snapshot3d(g);
	}
	return 0;
}
