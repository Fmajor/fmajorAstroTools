#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "htm.h"

void
tests(void){
    printf("Test for htmN\n"); {
        for (int i=0; i<=20; i++)
            printf("%-2d: %ld\n", i, htmN(i));
    }

    printf("Test for idRange\n"); {
        long int *temp;
        for (int i=0; i<=20; i++) {
            temp = idRange(i);
            printf("%-2d: %ld, %ld\n", i, temp[0], temp[1]);
        }
    }

    printf("Test for intLog2\n");
        //for (int i=0; i<=1000; i++)
            //printf("%-2d: %d\n", i, intLog2(i));

    printf("Test for printBin\n"); {
        for (int i=0; i<=200; i++) {
            printf("%-4d: ", i);
            printBin(i);
            printf("\n");
        }
    }

    //printf("Test for htmDepth\n"); {
        //long int *temp;
        //for (int i=0; i<=1024; i++) {
            //temp = idRange(i);
            //printf("%-4d: %d\n", i, htmDepth(i));
        //}
    //}

//    printf("Test for subdivide\n");{
//        long int *temp;
//        temp = subdivide(1023, 1); printf("%ld %ld\n", temp[0], temp[1]);
//        temp = subdivide(1023, 2); printf("%ld %ld\n", temp[0], temp[1]);
//        temp = subdivide(1023, 3); printf("%ld %ld\n", temp[0], temp[1]);
//        temp = subdivide(1023, 4); printf("%ld %ld\n", temp[0], temp[1]);
//        temp = subdivide(1023,17); printf("%ld %ld\n", temp[0], temp[1]);
//    }
//    printf("Test for fathers\n");{
//        long int *temp;
//        long int tempi = 17179869183;
//        /*long int tempi = 17592186044415;*/
//        temp = fathers(tempi); printf("%ld (init)\n", tempi);
//        for (int i=1; i<temp[0]; i++) {
//            printf("%ld\n", temp[i]);
//        }
//    }

//    printf("Test for upperTriangle"); {
//        for (int i=0; i<=2000; i++) {
//            printf("%-4d: ", i);
//            printBin(i);
//            printf(" ");
//            printf("%2d\n", upperTriangle(i));
//        }
//    }
    printf("Test for inTriangle\n");{
        struct Vector3 A={1,0,0};
        struct Vector3 Am; neg(&Am, A);
        struct Vector3 B={0,1,0};
        struct Vector3 Bm; neg(&Bm, B);
        struct Vector3 C={0,0,1};
        struct Vector3 Cm; neg(&Cm, C);
        struct threePoint a8[8]={
            {A, B, C}, {B, Am, C} ,{Am, Bm, C} ,{Bm, A, C},
            {A, B, Cm},{B, Am, Cm},{Am, Bm, Cm},{Bm, A, Cm}};
        struct Vector3 P={143,-113,-1321};
        /*productDouble(&P, P, 600);*/
        int temp;
        int dirSigns[8];
        int results[8];
        for (int i=0; i<8; i++) {
            dirSigns[i] = upperTriangle(i+8);
            /*temp = inTriangle(a8[i].p1, a8[i].p2, a8[i].p3, P, dirSigns[i]);*/
            temp = inTriangle(a8[i].p1, a8[i].p2, a8[i].p3, P);
            printf("%2d: %d\n", i, temp);
        }
        /*inMultiTriangle(results, P, 8, a8, dirSigns);*/
        inMultiTriangle(results, P, 8, a8);
        for (int i=0; i<8; i++) {
            printf("%d: %d\n", i, results[i]);
        }
        printf("last line\n");
        double * dtemp;
        dtemp = (double *)a8;
        for (int i=0; i<8; i++) {
            printf("%d\n",i);
            for (int j=0; j<3; j++) {
                printf("\t");
                for (int k=0; k<3; k++) {
                    printf("%f ", dtemp[i*9 + j*3 + k]);
                }
                printf("\n");
            }
            
        }
    }
}

int main(int argc,char *argv[]) {
    tests();
    return 0;
}

