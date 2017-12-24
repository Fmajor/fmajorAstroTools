#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void correlate(double* v1, int n1, double* v2, int n2, int* lag, int nlag, double* value) {
    int i, i1, i2, iStart;
    for (i=0; i<nlag; i++) {
        value[i]=0;
        if (lag[i]<0) {
            iStart = -lag[i];
            for (i1=0, i2=iStart; i2<n2 && i1<n1; i1++, i2++)
                value[i] += v1[i1] * v2[i2];
        } else {
            iStart = lag[i];
            for (i1=iStart, i2=0; i2<n2 && i1<n1; i1++, i2++)
                value[i] += v1[i1] * v2[i2];
        }
    }
}

void buildLamp(double* lam, int nLam, double* modLam, int nx, int npad, double* intensity, double* model) {
    int iline, iloc;
    double dpix, dx;
    if (modLam[1] > modLam[0]) {
        for (iline=0; iline<nLam; iline++) {
            iloc = 0;
            while (lam[iline]>modLam[iloc]) {
                if (iloc>=nx) {
                    iloc=0;
                    break;
                } else {
                    iloc++;
                }}
            if (iloc>0) {
                iloc = iloc-1;
                dx = lam[iline] - modLam[iloc];
                if (iloc>0 && iloc<nx-2) {
                    dpix = modLam[iloc+1] - modLam[iloc];
                    model[npad+iloc] += intensity[iline] * (1-dx/dpix);
                    model[npad+iloc+1] += intensity[iline] * dx/dpix;
                }
            }
            //printf("iLine: %5d  loc: %5d  dx:%8g  dpix:%8g\n", iline, iloc, dx, dpix);
        }
    } else {
        for (iline=0; iline<nLam; iline++) {
            iloc = 0;
            while (lam[iline]<modLam[iloc]) {
                if (iloc>=nx) {
                    iloc=0;
                    break;
                } else {
                    iloc++;
                }}
            if (iloc>0) {
                iloc = iloc-1;
                dx = modLam[iloc] - lam[iline];
                if (iloc>0 && iloc<nx-2) {
                    dpix = modLam[iloc] - modLam[iloc+1];
                    model[npad+iloc] += intensity[iline] * (1-dx/dpix);
                    model[npad+iloc+1] += intensity[iline] * dx/dpix;
                }
            }
            //printf("iLine: %5d  loc: %5d  dx:%8g  dpix:%8g\n", iline, iloc, dx, dpix);
        }
    }
}

void getSubIndex(long ii, long* nsteps, long n, long* result) {
    long *toMod, *toDiv;
    long i;
    toMod = (long *)malloc(sizeof(long)*n);
    toDiv = (long *)malloc(sizeof(long)*n);
    toMod[0] = nsteps[0];
    toDiv[0] = 1;
    for (i=0; i<n-1; i++) {
        toMod[i+1] = toMod[i] * nsteps[i+1];
        toDiv[i+1] = toDiv[i] * nsteps[i];
    }
    for (i=n-1; i>=0; i--) {
        result[i] = (ii%toMod[i])/toDiv[i];
    }
}