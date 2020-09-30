#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "fdtd-macro.h"

static int  frameX = 0, frameY = 0, frameZ = 0;
static char basename[80] = "plane";

void snapshot3d(Grid* g)
{
    int mm, nn, pp;
    float dim1, dim2, temp;
    char filename[100];
    FILE* out;
    if (Time >= 0) {
        int isxplane = 0;              //记录x平面 1->记录  0->不记录
        int isyplane = 1;
        int iszplane = 1;

        if (isxplane == 1) {

            sprintf(filename, "./data/%s-x.%d", basename, frameX++);
            out = fopen(filename, "wb");
            dim1 = SizeY;
            dim2 = SizeZ;
            fwrite(&dim1, sizeof(float), 1, out);
            fwrite(&dim2, sizeof(float), 1, out);

            mm = (SizeX) / 2;
            for (nn = SizeY - 1; nn >= 0; nn--) {
                for (pp = 0; pp < SizeZ; pp++) {
                    temp = (float)Ez(mm, nn, pp);
                    fwrite(&temp, sizeof(float), 1, out);
                }
            }
            fclose(out);
        }

        if (isyplane == 1) {
            sprintf(filename, "./data/%s-y.%d", basename, frameY++);
            out = fopen(filename, "wb");
            dim1 = SizeX - 1;
            dim2 = SizeZ;
            fwrite(&dim1, sizeof(float), 1, out);
            fwrite(&dim2, sizeof(float), 1, out);

            nn = SizeY / 2;
            for (mm = SizeX - 1; mm >= 0; mm--) {
                for (pp = 0; pp < SizeZ - 1; pp++) {
                    temp = (float)Ez(mm, nn, pp);
                    fwrite(&temp, sizeof(float), 1, out);
                }
            }
            fclose(out);
        }

        if (iszplane == 1) {
            sprintf(filename, "./data/%s-z.%d", basename, frameZ++);
            out = fopen(filename, "wb");
            dim1 = SizeX - 1;
            dim2 = SizeY;
            fwrite(&dim1, sizeof(float), 1, out);
            fwrite(&dim2, sizeof(float), 1, out);

            pp = 1;
            for (mm = SizeX - 1; mm >= 0; mm--) {
                for (nn = 0; nn < SizeY - 1; nn++) {
                    temp = (float)(Ez(mm, nn, pp));
                    fwrite(&temp, sizeof(float), 1, out);
                }
            }
            fclose(out);
        }
    }
}
