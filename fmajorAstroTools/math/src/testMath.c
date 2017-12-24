#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "astroMath.h"

#define dPrintArray(a, n){ for(int tempIndex=0;tempIndex<n;tempIndex++) printf("%3f ", a[tempIndex]); printf("\n"); }
#define fPrintArray(a, n){ for(int tempIndex=0;tempIndex<n;tempIndex++) printf("%3f ", a[tempIndex]); printf("\n"); }
#define iPrintArray(a, n){ for(int tempIndex=0;tempIndex<n;tempIndex++) printf("%5i ", a[tempIndex]); printf("\n");}
#define dFormatPrintArray(a, n, format){ for(int tempIndex=0;tempIndex<n;tempIndex++) printf(format, a[tempIndex]); printf("\n"); }
/* example of usage
    dPrintArray(testData, 10);
    dFormatPrintArray(testData, 10, "%3f ");
end of example of usage*/




int main(void)
{
    printf("isfinite(NAN)         = %d\n", isfinite(NAN));
    printf("isfinite(INFINITY)    = %d\n", isfinite(INFINITY));
    printf("isfinite(0.0)         = %d\n", isfinite(0.0));
    printf("isfinite(DBL_MIN/2.0) = %d\n", isfinite(DBL_MIN/2.0));
    printf("isfinite(1.0)         = %d\n", isfinite(1.0));
    printf("isfinite(exp(800))    = %d\n", isfinite(exp(800)));
    double testData[]={1,9,2,8,3,7,4,6,5,-43};
    int indexArray[]={-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
    /*printf("test data:\n\t");*/
    /*dPrintArray(testData, 10);*/
    /*printf("median: %f\n", percentile(testData, 10, 0.1));*/
    printf("!!!!!\n");
    dPrintArray(testData, 10);
    doubleIndexSort(testData, 10, indexArray);
    dPrintArray(testData, 10);
    iPrintArray(indexArray, 10);
}
