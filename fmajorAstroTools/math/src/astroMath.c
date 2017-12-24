#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>
// isfinite()

int randomized_partition(double A[], int p, int r);
int randomInt(int p, int r);
int partition(double A[], int p, int r);

double randomized_select(double A[], int p, int r, int i) {
    int q, k;
    if (p==r)
        return A[p];
    q = randomized_partition(A, p, r);
    k = q-p+1;
    if (i==k)
        return A[q];
    else if (i<k)
        return randomized_select(A, p, q-1, i);
    else
        return randomized_select(A, q+1, r, i-k);
}

double median(double A[], int n) {
    A--; // make it work for A[0..n-1]
    int i;
    double result;
    if (n%2==0) {
        result = (randomized_select(A, 1, n, n/2) + randomized_select(A, 1, n, n/2+1))/2.0;
    } else {
        i = (int) (1 + (n-1)/2);
        result = randomized_select(A, 1, n, i);
    }
    A++;
    return result;
}

//!!  sub of A ~ [1,...n]
double percentile(double A[], int n, double p) {
    A--; // make it work for A[0..n-1]
    int i;
    double result;
    if (p<0 || p>1) {
        printf("in percentile: p = %10.5g should be in between 0 and 1\n", p);
        exit(1);
    }
    i = (int) (1 + (n-1) * p);
    result = randomized_select(A, 1, n, i);
    A++;
    return result;
}

#undef SWAP
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
int randomized_partition(double A[], int p, int r) {
    int i;
    double temp;
    i=randomInt(p, r);
    SWAP(A[r], A[i])
    return partition(A, p, r);
}

int partition(double A[], int p, int r){
    double x=A[r];
    double temp;
    int i=p-1;
    for (int j=p; j<r; j++) {
        if (A[j]<=x) {
            i++;
            SWAP(A[i], A[j])
        }
    }
    SWAP(A[i+1], A[r]);
    return i+1;
}
#undef SWAP

int randomInt(int p, int r) {
    srand(0);
    return (int)(((double)rand()/RAND_MAX)*(r-p) + p);
}

// the sub of A is [0, n-1]
int doubleMin(const double A[], int n, double *result) {
    *result = INFINITY;
    int index;
    for (int i=0; i<n; i++) {
        if (A[i]<*result) {
            *result = A[i];
            index = i;
        }
    }
    return index;
}

// the sub of A is [0, n-1]
int doubleMax(const double A[], int n, double *result) {
    *result = -INFINITY;
    int index;
    for (int i=0; i<n; i++) {
        if (A[i]>*result) {
            *result = A[i];
            index = i;
        }
    }
    return index;
}

double doubleSum(const double A[], int n) {
    double result=0;
    for (int i=0; i<n; i++)
        result += A[i];
    return result;
}
double doubleSub(const double A[], const double B[], double C[], int n) {
    double result=0;
    for (int i=0; i<n; i++)
        result += A[i];
    return result;
}
int shortSum(const short A[], int n) {
    int result=0;
    for (int i=0; i<n; i++)
        result += A[i];
    return result;
}
void shortNot(const short A[], short B[], int n){
    for (int i=0; i<n; i++)
        B[i] = !A[i];
}

struct doubleIndexArray {
    int index;
    double value;
};

int compareDoubleIndexArray(const void *a, const void *b) {
    return ((*(struct doubleIndexArray*)a).value > (*(struct doubleIndexArray*)b).value) ? 1 : -1;
}

// only have index array, origin array not sorted
void doubleIndexSort(const double a[], int n, int indexArray[]) {
    struct doubleIndexArray *obj;
    obj = (struct doubleIndexArray *)malloc(sizeof(struct doubleIndexArray)*n);
    for (int i=0;i<n;i++) {
        obj[i].index = i;
        obj[i].value = a[i];
    }
    qsort(obj, n, sizeof(struct doubleIndexArray), compareDoubleIndexArray);
    for (int i=0;i<n;i++) {
        indexArray[i] = obj[i].index;
    }
    free(obj);
}

// masktype==0 badmask, masktype==1 goodmask
double maskMedian(const double *data, short *mask, double *tempdata, int N, int masktype) {
    int i, goodN=0;
    if (masktype) {// goodmask
        for (i=0;i<N;i++) {
            if (mask[i]) {
                tempdata[goodN++] = data[i];
            }
        }
    } else {// badmask
        for (i=0;i<N;i++) {
            if (!mask[i]) {
                tempdata[goodN++] = data[i];
            }
        }
    }
    return median(tempdata, goodN);
}


//!! data and invvar must be clear from nan and inf
// invvar[0] == -1 =====> no invvar set
//!! mask is the init data badmask
//!!    it must have been updated by invvar
//!! new badMask is the new add badmask after sigma rejection

int iterstat1d(const double data[], int M,
               const double invvar[], double lsigrej, double hsigrej,
               int maxiter, short const initmask[], short badmask[], short newBadMask[],
               short stdRej, short useMedian, int maxRej, double maxDev,
               double result[]) {
    //#region: variables
    int niter, goodN, i, j;
    int newBadMaskNum;
    double numerator, denominator, mean, std, sigrej;
    double *newBadValue, *diff, aux, *badness;
    short *goodmask, allBad;
    int *newBadIndex;
    int *auxIndex, *newBadValueIndexTodo, newBadValueIndexTodoNum;
    int aux0, aux1;
    double *medianTemp, thisdiffsub;
    //#endregion: variables
    //#region: malloc
    diff = (double *)malloc(sizeof(double)*M);
    badness = (double *)malloc(sizeof(double)*M);
    goodmask = (short *)malloc(sizeof(short)*M);
    newBadIndex = (int *)malloc(sizeof(int)*M);
    newBadValue = (double *)malloc(sizeof(double)*M);
    auxIndex = (int *)malloc(sizeof(int)*M);
    newBadValueIndexTodo = (int *)malloc(sizeof(int)*M);
    medianTemp = (double *)malloc(sizeof(double)*M);
    //#endregion: malloc
    //#region: !! first cycle
    //!! here data number must > 1
    allBad = 0;
    shortNot(badmask, goodmask, M);
    goodN = shortSum(goodmask, M);
    //!! all bad
    if (goodN==0) {
        mean = NAN;
        std = NAN;
        allBad = 1;
    //!! goodN > 0
    } else {
        if (invvar[0] != -1) {// invvar set, weight mean
            numerator = 0;
            denominator = 0;
            for (int i=0; i<M; i++) {
                numerator += data[i] * invvar[i];
                denominator += invvar[i];
            }
            mean = numerator / denominator;
        } else {// normal mean
            numerator = 0;
            for (int i=0; i<M; i++)
                numerator += data[i] * goodmask[i];
            mean = numerator / goodN;
        }

        // goodN > 1, can cal std
        if (goodN>1) {
            numerator = 0;
            for (int i=0; i<M; i++)
                numerator += goodmask[i] * (data[i] - mean) * (data[i] - mean);
            std = sqrt(numerator / (goodN - 1));
        // goodN == 1, can not cal std
        } else {
            std = NAN;
        }
    }
    if isfinite(std) {
        if (useMedian) {
            thisdiffsub = maskMedian(data, badmask, medianTemp, M, 0);
        } else {
            thisdiffsub = mean;
        }
        if (stdRej)
            for (int i=0; i<M; i++) {
                diff[i] = fabs(data[i] - thisdiffsub);
                sigrej = data[i] - thisdiffsub >0 ? hsigrej : lsigrej;
                badness[i] = diff[i] / std > sigrej ? diff[i] / std - sigrej : 0;
            }
        else
            for (int i=0; i<M; i++) {
                diff[i] = fabs(data[i] - thisdiffsub);
                sigrej = data[i] - thisdiffsub >0 ? hsigrej : lsigrej;
                badness[i] = diff[i] * sqrt(invvar[i]) > sigrej ? diff[i] * sqrt(invvar[i]) - sigrej : 0;
            }
        if (maxDev!=-1) {
            for (int i=0; i<M; i++) {
                diff[i] = diff[i] > maxDev ? (diff[i]-maxDev)/maxDev : 0;
                badness[i] += diff[i];
            }
        }
        newBadMaskNum = 0;
        for (int i=0; i<M; i++) { // update new bad mask
            newBadMask[i] = 0;
            if (badness[i] && !badmask[i]) {
                newBadMask[i] = 1;
                newBadMaskNum++;
            }
        }
    } else {
        newBadMaskNum = 0;
    }
    #ifdef DEBUGALL
    printf("\t c:      iter: %d mean: %.10f std: %.10f badnumber: %d newBadNumber: %d\n", niter, mean, std, shortSum(badmask, M), newBadMaskNum);
    #endif
    //#endregion: first cycle
    //!! main loop
    niter = 0;
    while (niter<maxiter && newBadMaskNum) {
        niter++;
        //!! do the rejection
        if (maxRej != -1) {
            for (i=0, j=0; i<M; i++) {
                if (newBadMask[i]) {
                    newBadIndex[j] = i;
                    newBadValue[j++] = badness[i];
                }
            }
            if (maxRej==1) {
                newBadValueIndexTodo[0] = newBadIndex[doubleMax(newBadValue, newBadMaskNum, &aux)];
                newBadValueIndexTodoNum = 1;
            } else {
                newBadValueIndexTodoNum = maxRej > newBadMaskNum ? newBadMaskNum : maxRej;
                doubleIndexSort(newBadValue, newBadMaskNum, auxIndex);
                for (i=newBadMaskNum - newBadValueIndexTodoNum, j=0; i<newBadMaskNum; i++, j++) {
                    newBadValueIndexTodo[j] = newBadIndex[auxIndex[i]];
                }
            }
            // update badmask accroding to maxRej
            #ifdef DEBUGALL
            printf("\t\t");
            #endif
            for (int i=0; i<newBadValueIndexTodoNum; i++) {
                aux0 = newBadValueIndexTodo[i];
                #ifdef DEBUGALL
                printf("%d:%.3f ", aux0, badness[aux0]);
                #endif
                badmask[aux0] = 1;
            }
            #ifdef DEBUGALL
            printf("\n");
            #endif

        } else {// directly update badmask
            #ifdef DEBUGALL
            printf("\t\t");
            #endif
            for (i=0; i<M; i++)
                if (newBadMask[i]) {
                    #ifdef DEBUGALL
                    printf("%d:%.3f ", aux0, badness[aux0]);
                    #endif
                    badmask[i] = 1;
                }
            #ifdef DEBUGALL
            printf("\n");
            #endif
        }

        shortNot(badmask, goodmask, M);
        goodN = shortSum(goodmask, M);
        //!! all bad
        if (goodN==0) {
            mean = NAN;
            std = NAN;
            allBad = 1;
            break;
        //!! goodN > 0
        } else {
            if (invvar[0] != -1) {// invvar set
                numerator = 0;
                denominator = 0;
                for (int i=0; i<M; i++) {
                    numerator += data[i] * invvar[i];
                    denominator += invvar[i];
                }
                mean = numerator / denominator;
            } else {
                numerator = 0;
                for (int i=0; i<M; i++)
                    numerator += data[i] * goodmask[i];
                mean = numerator / goodN;
            }

            // goodN > 1, can cal std
            if (goodN>1) {
                numerator = 0;
                for (int i=0; i<M; i++)
                    numerator += goodmask[i] * (data[i] - mean) * (data[i] - mean);
                std = sqrt(numerator / (goodN - 1));
            // goodN == 1, can not cal std
            } else {
                std = NAN;
                break;
            }
        }
        if (useMedian) {
            thisdiffsub = maskMedian(data, badmask, medianTemp, M, 0);
        } else {
            thisdiffsub = mean;
        }

        if (stdRej)
            for (int i=0; i<M; i++) {
                diff[i] = fabs(data[i] - thisdiffsub);
                sigrej = data[i] - thisdiffsub >0 ? hsigrej : lsigrej;
                badness[i] = diff[i] / std > sigrej ? diff[i] / std - sigrej : 0;
            }
        else
            for (int i=0; i<M; i++) {
                diff[i] = fabs(data[i] - thisdiffsub);
                sigrej = data[i] - thisdiffsub >0 ? hsigrej : lsigrej;
                badness[i] = diff[i] * sqrt(invvar[i]) > sigrej ? diff[i] * sqrt(invvar[i]) - sigrej : 0;
            }
        if (maxDev!=-1) {
            for (int i=0; i<M; i++) {
                diff[i] = diff[i] > maxDev ? (diff[i]-maxDev)/maxDev : 0;
                badness[i] += diff[i];
            }
        }
        newBadMaskNum = 0;
        for (int i=0; i<M; i++) { // update new bad mask
            newBadMask[i] = 0;
            if (badness[i] && !badmask[i]) {
                newBadMask[i] = 1;
                newBadMaskNum++;
            }
        }

        #ifdef DEBUGALL
        printf("\t c:      iter: %d mean: %.10f std: %.10f badnumber: %d newBadNumber: %d\n", niter, mean, std, shortSum(badmask, M), newBadMaskNum);
        #endif
    }
    //#region: gen results
    result[0] = mean;
    result[1] = std;
    for (i=0; i<M; i++) {
        newBadMask[i] = badmask[i] - initmask[i];
    }
    if (allBad) {
        result[2] = NAN;
    } else {
        result[2] = maskMedian(data, badmask, medianTemp, M, 0);
    }
    //#endregion: gen results
    //#region: free
    free(diff);
    free(badness);
    free(goodmask);
    free(newBadIndex);
    free(newBadValue);
    free(auxIndex);
    free(newBadValueIndexTodo);
    free(medianTemp);
    //#endregion: free
    return niter;
}

int iterstat1dCycle(const double data[], int M,
               const double invvar[], double lsigrej, double hsigrej,
               int maxiter, short const initmask[], short badmask[], short newBadMask[],
               short stdRej, short useMedian, int maxRej, double maxDev,
               double result[], double* badness,
               double *diff, short *goodmask, int *newBadIndex, double *newBadValue,
               int *auxIndex, int *newBadValueIndexTodo) {
    //#region: variables
    int niter, goodN, i, j;
    int newBadMaskNum;
    double numerator, denominator, mean, std, sigrej;
    double aux;
    short allBad;
    int newBadValueIndexTodoNum;
    int aux0, aux1;
    double *medianTemp, thisdiffsub;
    medianTemp = (double *)malloc(sizeof(double)*M);
    //#endregion: variables
    //#region: !! first cycle
    //!! here data number must > 1
    allBad = 0;
    shortNot(badmask, goodmask, M);
    goodN = shortSum(goodmask, M);
    //!! all bad
    if (goodN==0) {
        mean = NAN;
        std = NAN;
        allBad = 1;
    //!! goodN > 0
    } else {
        if (invvar[0] != -1) {// invvar set
            numerator = 0;
            denominator = 0;
            for (int i=0; i<M; i++) {
                numerator += data[i] * invvar[i];
                denominator += invvar[i];
            }
            mean = numerator / denominator;
        } else {
            numerator = 0;
            for (int i=0; i<M; i++)
                numerator += data[i] * goodmask[i];
            mean = numerator / goodN;
        }

        // goodN > 1, can cal std
        if (goodN>1) {
            numerator = 0;
            for (int i=0; i<M; i++)
                numerator += goodmask[i] * (data[i] - mean) * (data[i] - mean);
            std = sqrt(numerator / (goodN - 1));
        // goodN == 1, can not cal std
        } else {
            std = NAN;
        }
    }
    if isfinite(std) {
        if (useMedian) {
            thisdiffsub = maskMedian(data, badmask, medianTemp, M, 0);
        } else {
            thisdiffsub = mean;
        }
        if (stdRej)
            for (int i=0; i<M; i++) {
                diff[i] = fabs(data[i] - thisdiffsub);
                sigrej = data[i] - thisdiffsub >0 ? hsigrej : lsigrej;
                badness[i] = diff[i] / std > sigrej ? diff[i] / std - sigrej : 0;
            }
        else
            for (int i=0; i<M; i++) {
                diff[i] = fabs(data[i] - thisdiffsub);
                sigrej = data[i] - thisdiffsub >0 ? hsigrej : lsigrej;
                badness[i] = diff[i] * sqrt(invvar[i]) > sigrej ? diff[i] * sqrt(invvar[i]) - sigrej : 0;
            }
        if (maxDev!=-1) {
            for (int i=0; i<M; i++) {
                diff[i] = diff[i] > maxDev ? diff[i]-maxDev : 0;
                badness[i] += diff[i];
            }
        }
        newBadMaskNum = 0;
        for (int i=0; i<M; i++) { // update new bad mask
            newBadMask[i] = 0;
            if (badness[i] && !badmask[i]) {
                newBadMask[i] = 1;
                newBadMaskNum++;
            }
        }
    } else {
        newBadMaskNum = 0;
    }
    #ifdef DEBUGALL
    printf("\t c:      iter: %d mean: %.10f std: %.10f badnumber: %d newBadNumber: %d\n", niter, mean, std, shortSum(badmask, M), newBadMaskNum);
    #endif
    //#endregion: first cycle
    //!! main loop
    niter = 0;
    while (niter<maxiter && newBadMaskNum) {
        niter++;
        //!! do the rejection
        if (maxRej != -1) {
            for (i=0, j=0; i<M; i++) {
                if (newBadMask[i]) {
                    newBadIndex[j] = i;
                    newBadValue[j++] = badness[i];
                }
            }
            if (maxRej==1) {
                newBadValueIndexTodo[0] = newBadIndex[doubleMax(newBadValue, newBadMaskNum, &aux)];
                newBadValueIndexTodoNum = 1;
            } else {
                newBadValueIndexTodoNum = maxRej > newBadMaskNum ? newBadMaskNum : maxRej;
                doubleIndexSort(newBadValue, newBadMaskNum, auxIndex);
                for (i=newBadMaskNum - newBadValueIndexTodoNum, j=0; i<newBadMaskNum; i++, j++) {
                    newBadValueIndexTodo[j] = newBadIndex[auxIndex[i]];
                }
            }
            // update badmask accroding to maxRej
            #ifdef DEBUGALL
            printf("\t\t");
            #endif
            for (int i=0; i<newBadValueIndexTodoNum; i++) {
                aux0 = newBadValueIndexTodo[i];
                #ifdef DEBUGALL
                printf("%d ", aux0);
                #endif
                badmask[aux0] = 1;
            }
            #ifdef DEBUGALL
            printf("\n");
            #endif

        } else {// directly update badmask
            for (i=0; i<M; i++)
                if (newBadMask[i])
                    badmask[i] = 1;
        }

        shortNot(badmask, goodmask, M);
        goodN = shortSum(goodmask, M);
        //!! all bad
        if (goodN==0) {
            mean = NAN;
            std = NAN;
            allBad = 1;
            break;
        //!! goodN > 0
        } else {
            if (invvar[0] != -1) {// invvar set
                numerator = 0;
                denominator = 0;
                for (int i=0; i<M; i++) {
                    numerator += data[i] * invvar[i];
                    denominator += invvar[i];
                }
                mean = numerator / denominator;
            } else {
                numerator = 0;
                for (int i=0; i<M; i++)
                    numerator += data[i] * goodmask[i];
                mean = numerator / goodN;
            }

            // goodN > 1, can cal std
            if (goodN>1) {
                numerator = 0;
                for (int i=0; i<M; i++)
                    numerator += goodmask[i] * (data[i] - mean) * (data[i] - mean);
                std = sqrt(numerator / (goodN - 1));
            // goodN == 1, can not cal std
            } else {
                std = NAN;
                break;
            }
        }
        if (useMedian) {
            thisdiffsub = maskMedian(data, badmask, medianTemp, M, 0);
        } else {
            thisdiffsub = mean;
        }

        if (stdRej)
            for (int i=0; i<M; i++) {
                diff[i] = fabs(data[i] - thisdiffsub);
                sigrej = data[i] - thisdiffsub >0 ? hsigrej : lsigrej;
                badness[i] = diff[i] / std > sigrej ? diff[i] / std - sigrej : 0;
            }
        else
            for (int i=0; i<M; i++) {
                diff[i] = fabs(data[i] - thisdiffsub);
                sigrej = data[i] - thisdiffsub >0 ? hsigrej : lsigrej;
                badness[i] = diff[i] * sqrt(invvar[i]) > sigrej ? diff[i] * sqrt(invvar[i]) - sigrej : 0;
            }
        if (maxDev!=-1) {
            for (int i=0; i<M; i++) {
                diff[i] = diff[i] > maxDev ? diff[i]-maxDev : 0;
                badness[i] += diff[i];
            }
        }
        newBadMaskNum = 0;
        for (int i=0; i<M; i++) { // update new bad mask
            newBadMask[i] = 0;
            if (badness[i] && !badmask[i]) {
                newBadMask[i] = 1;
                newBadMaskNum++;
            }
        }

        #ifdef DEBUGALL
        printf("\t c:      iter: %d mean: %.10f std: %.10f badnumber: %d newBadNumber: %d\n", niter, mean, std, shortSum(badmask, M), newBadMaskNum);
        #endif
    }
    //#region: gen results
    result[0] = mean;
    result[1] = std;
    for (i=0; i<M; i++) {
        newBadMask[i] = badmask[i] - initmask[i];
    }
    if (allBad) {
        result[2] = NAN;
    } else {
        result[2] = maskMedian(data, badmask, medianTemp, M, 0);
    }
    free(medianTemp);
    //#endregion: gen results
    return niter;
}

void iterstat3d(double data[], int c, int b, int M,
                double invvar[], double lsigrej, double hsigrej,
                int maxiter, short initmask[], short badmask[], short newBadMask[],
                short stdRej, short useMedian, int maxRej, double maxDev,
                double result[]) {
    //#region: variables
    int niter, goodN, i, j, eachb, eachc;
    int resultStep0, resultStep1;
    int newBadMaskNum;
    double numerator, denominator, mean, std;
    double *newBadValue, *diff, aux, *badness;
    short *goodmask, allBad;
    int *newBadIndex;
    int *auxIndex, *newBadValueIndexTodo, newBadValueIndexTodoNum;
    int aux0, aux1;

    double *datap, *invvarp, *resultp;
    short *initmaskp, *badmaskp, *newBadMaskp;
    int dataStep0, dataStep1;
    resultStep0 = 3*b;
    resultStep1 = 3;
    dataStep0 = M*b;
    dataStep1 = M;
    //#endregion: variables
    //#region: malloc
    diff = (double *)malloc(sizeof(double)*M);
    badness = (double *)malloc(sizeof(double)*M);
    goodmask = (short *)malloc(sizeof(short)*M);
    newBadIndex = (int *)malloc(sizeof(int)*M);
    newBadValue = (double *)malloc(sizeof(double)*M);
    auxIndex = (int *)malloc(sizeof(int)*M);
    newBadValueIndexTodo = (int *)malloc(sizeof(int)*M);
    //#endregion: malloc

    datap = data; invvarp = invvar; resultp = result;
    initmaskp = initmask; badmaskp = badmask; newBadMaskp = newBadMask;

    for (eachc=0; eachc<c; eachc++) {
        for (eachb=0; eachb<b; eachb++,
                datap += dataStep1,
                invvarp += dataStep1,
                resultp += resultStep1,
                initmaskp += dataStep1,
                badmaskp += dataStep1,
                newBadMaskp += dataStep1) {
                #ifdef DEBUGALL
                    printf("%d/%d, %d/%d\n", eachc, c, eachb, b);
                #endif
            iterstat1dCycle(datap, M,
                   invvarp, lsigrej, hsigrej,
                   maxiter, initmaskp, badmaskp, newBadMaskp,
                   stdRej, useMedian, maxRej, maxDev,
                   resultp, badness,
                   diff, goodmask, newBadIndex, newBadValue,
                   auxIndex, newBadValueIndexTodo);
        }
    }
    //#region: free
    free(diff);
    free(badness);
    free(goodmask);
    free(newBadIndex);
    free(newBadValue);
    free(auxIndex);
    free(newBadValueIndexTodo);
    //#endregion: free
}

// 0,1,2,3,4,5,6,7,8,9
// <0, return 0
// >8, return 8
int sortedFindLeft(double thisx, int start, double x[], int N) {
    int tryRight = 2;
    int i, mid;
    int end=N-1;
    if (thisx<x[0]) return 0;
    if (thisx>x[N-2]) return N-2;
    for (i=0; i<tryRight; i++) {
        if (thisx<x[start+1])
            return start;
        start++;
    }
    while (end>start+1) {
        mid = (start+end)/2;
        if (thisx>x[mid])
            start = mid;
        else
            end = mid;
    }
    return start;

}
void getLeftIndex(double x[], int N, double xx[], int indexs[], int M) {
    int i,j;
    int start = 0;
    for (i=0;i<M;i++){
        start = sortedFindLeft(xx[i], start, x, N);
        indexs[i] = start;
    }
}


// x and xx should be sorted
void linearInterpWithError(double x[], double y[], double sigma[], int N,
                     double xx[], double yy[], double sigmasigma[], int M){
    //#region: variables
    int *indexs, i, ii;
    double hi;
    //#endregion: variables
    //#region: malloc
    indexs = (int *)malloc(sizeof(int)*M);
    //#endregion: malloc

    getLeftIndex(x, N, xx, indexs, M);
    for (i=0;i<M;i++) {
        ii = indexs[i];
        hi = (x[ii+1] - xx[i]) / (x[ii+1] - x[ii]);
        yy[i] = hi * y[ii] + (1-hi) * y[ii+1];
        sigmasigma[i] = sqrt(hi * sigma[ii] * sigma[ii] + (1-hi) * sigma[ii+1] * sigma[ii+1]);
    }
    //#region: free
    free(indexs);
    //#endregion: free
}

void yAx(const double y[], const double A[], const double x[],
         int N, int xOrder, int yOrder,
         double result[]) {

    int i,j,k,n,m,coeffOffset,xOffset,yOffset;
    //printf("in c:\n\t%d %d %d %.15lf %.15lf %.15lf %.15lf\n",
    //                  N, xOrder, yOrder, y[3], A[3], x[3], result[3]);
    //for (i=0;i<3;i++){
    //    yOffset = i*yOrder;
    //    for (j=0;j<yOrder;j++)
    //        printf("%lf ",y[yOffset+j]);
    //    printf("\n");
    //}

    for (i=0; i<N; i++) {
        xOffset = i*xOrder;
        yOffset = i*yOrder;
        for (n=0; n<yOrder; n++) {
            coeffOffset = n*xOrder;
            for (m=0; m<xOrder; m++) {
                result[i] += A[coeffOffset + m] * x[xOffset + m] * y[yOffset + n];
            }
        }
    }
}

//!! interpolate functions
struct Vertex {
    double x;
    double y;
    struct Vertex *prev;
    struct Vertex *next;
    int flag;
};
struct Polygon {
    struct Vertex *P;
    int num;
};

struct Vector2 {
    double x;
    double y;
};
struct Vector3 {
    double x;
    double y;
    double z;
};

struct Vector2* newVector2(double x, double y) {
    struct Vector2 *newV = (struct Vector2 *)malloc(sizeof(struct Vector2));
    newV->x = x;
    newV->y = y;
    return newV;
}
double dot2(struct Vector2* a, struct Vector2* b) {
    return a->x*b->x + a->y*b->y;
}
struct Vector2* plus2(struct Vector2* a, struct Vector2* b) {
    struct Vector2 *newV = (struct Vector2 *)malloc(sizeof(struct Vector2));
    newV->x = a->x + b->x;
    newV->y = a->y + b->y;
    return newV;
}
struct Vector2* minus2(struct Vector2* a, struct Vector2* b) {
    struct Vector2 *newV = (struct Vector2 *)malloc(sizeof(struct Vector2));
    newV->x = a->x - b->x;
    newV->y = a->y - b->y;
    return newV;
}
struct Vector2* neg2(struct Vector2* v) {
    struct Vector2 *newV = (struct Vector2 *)malloc(sizeof(struct Vector2));
    newV->x = -v->x;
    newV->y = -v->y;
    return newV;
}



double norm2(struct Vector2* v) {
    return sqrt(dot2(v,v));
}

struct Vector3* newVector3(double x, double y, double z) {
    struct Vector3 *newV = (struct Vector3 *)malloc(sizeof(struct Vector3));
    newV->x = x;
    newV->y = y;
    newV->z = z;
    return newV;
}

double dot3(struct Vector3* a, struct Vector3* b) {
    return a->x*b->x + a->y*b->y + a->z*b->z;
}
double norm3(struct Vector3* v) {
    return sqrt(dot3(v,v));
}
struct Vector3* plus3(struct Vector3* a, struct Vector3* b) {
    struct Vector3 *newV = (struct Vector3 *)malloc(sizeof(struct Vector3));
    newV->x = a->x + b->x;
    newV->y = a->y + b->y;
    newV->z = a->z + b->z;
    return newV;
}
struct Vector3* minus3(struct Vector3* a, struct Vector3* b) {
    struct Vector3 *newV = (struct Vector3 *)malloc(sizeof(struct Vector3));
    newV->x = a->x - b->x;
    newV->y = a->y - b->y;
    newV->z = a->z - b->z;
    return newV;
}
struct Vector3* neg3(struct Vector3* v) {
    struct Vector3 *newV = (struct Vector3 *)malloc(sizeof(struct Vector3));
    newV->x = -v->x;
    newV->y = -v->y;
    newV->z = -v->z;
    return newV;
}

struct Vector3* cross3(struct Vector3* a, struct Vector3* b) {
    struct Vector3 *newV = (struct Vector3 *)malloc(sizeof(struct Vector3));
    newV->x = (a->y * b->z - a->z * b->y);
    newV->y = (a->z * b->x - a->x * b->z);
    newV->z = (a->x * b->y - a->y * b->x);
    return newV;
}

void setValue3(struct Vector3* v, double x, double y, double z) {
    v->x = x;
    v->y = y;
    v->z = z;
}
void setValue2(struct Vector3* v, double x, double y) {
    v->x = x;
    v->y = y;
}

struct Polygon* newPolygon(){
    // return pointer of the polygon
    struct Polygon *newObj = (struct Polygon *)malloc(sizeof(struct Polygon));
    newObj->P = NULL;
    newObj->num = 0;
    return newObj;
}

struct Vertex* newVertex(double x, double y, int flag) {
    struct Vertex *newV = (struct Vertex *)malloc(sizeof(struct Vertex));
    newV->x = x;
    newV->y = y;
    newV->flag = flag;
    return newV;
}

struct Vertex* addVertex(struct Polygon *obj, double x, double y, int flag, int index) {
    // return pointer of the inserted point
    struct Vertex *thisV = obj->P;
    struct Vertex *thisVprev, *thisVnext;
    if ((index<0)||(index>obj->num)) {
        printf("insert out of range %d/%d\n", index, obj->num);
        return NULL;
    }
    if (index==0) {// insert into first pos
        struct Vertex *newV = (struct Vertex *)malloc(sizeof(struct Vertex));
        newV->x = x;
        newV->y = y;
        newV->flag = flag;
        obj->P = newV;
        obj->num++;
        if (!thisV) {// insert into empty Polygon
            newV->prev = newV;
            newV->next = newV;
            return newV;
        } else {// insert into nonempty Polygon
            thisVprev = thisV->prev;
            newV->next = thisV;
            newV->prev = thisV->prev;
            thisV->prev = newV;
            thisVprev->next = newV;
            return newV;
        }
    } else {// insert into latter pos
        if (!thisV) {
            printf("Can not insert into >0 pos in a empty Vertex");
            return NULL;
        } else {
            struct Vertex *newV = (struct Vertex *)malloc(sizeof(struct Vertex));
            newV->x = x;
            newV->y = y;
            newV->flag = flag;
            obj->num++;
            for (int i=0;i<index;i++) {// thisV point to insert pos
                thisV = thisV->next;
            }
            thisVprev = thisV->prev;
            newV->next = thisV;
            newV->prev = thisV->prev;
            thisV->prev = newV;
            thisVprev->next = newV;
            return newV;
        }
    }
}

/*void addVertexV(struct Polygon *obj, struct Vertex *newV, int index) {*/
    /*// return pointer of the inserted point*/
    /*struct Vertex *thisV = obj->P;*/
    /*struct Vertex *thisVprev, *thisVnext;*/
    /*if ((index<0)||(index>obj->num)) {*/
        /*printf("insert out of range %d/%d\n", index, obj->num);*/
        /*return NULL;*/
    /*}*/
    /*if (index==0) {// insert into first pos*/
        /*obj->P = newV;*/
        /*obj->num++;*/
        /*if (!thisV) {// insert into empty Polygon*/
            /*newV->prev = newV;*/
            /*newV->next = newV;*/
            /*return newV;*/
        /*} else {// insert into nonempty Polygon*/
            /*thisVprev = thisV->prev;*/
            /*newV->next = thisV;*/
            /*newV->prev = thisV->prev;*/
            /*thisV->prev = newV;*/
            /*thisVprev->next = newV;*/
            /*return newV;*/
        /*}*/
    /*} else {// insert into latter pos*/
        /*if (!thisV) {*/
            /*printf("Can not insert into >0 pos in a empty Vertex");*/
            /*return NULL;*/
        /*} else {*/
            /*obj->num++;*/
            /*for (int i=0;i<index;i++) {// thisV point to insert pos*/
                /*thisV = thisV->next;*/
            /*}*/
            /*thisVprev = thisV->prev;*/
            /*newV->next = thisV;*/
            /*newV->prev = thisV->prev;*/
            /*thisV->prev = newV;*/
            /*thisVprev->next = newV;*/
            /*return newV;*/
        /*}*/
    /*}*/
/*}*/

struct Vertex* delVertex(struct Polygon *obj, int index) {
    // return pointer of the next Vertex (return first if delete the last)
    struct Vertex *thisV = obj->P;
    struct Vertex *thisVprev, *thisVnext;
    if ((index<0)||(index>obj->num-1)) {
        printf("delete out of range %d/%d\n", index, obj->num);
        return NULL;
    }
    if (!thisV) {
        printf("delete Vertex form empty Polygon");
        return NULL;
    } else {
        for (int i=0;i<index;i++) {// thisV point to insert pos
            thisV = thisV->next;
        }
        thisVprev = thisV->prev;
        thisVnext = thisV->next;
        thisVprev->next = thisVnext;
        thisVnext->prev = thisVprev;
        free(thisV);
    }
    if (index==0) obj->P = thisVnext;
    obj->num--;
    if (!obj->num) obj->P = NULL;
    return thisVnext;
}

void freePolygon(struct Polygon *obj) {
    struct Vertex *thisV, *nextV, *initV;
    if (obj->P) {
       thisV = obj->P;
       initV = obj->P;
       do {
           nextV = thisV->next;
           free(thisV);
           thisV = nextV;
       } while(thisV!=initV);
    }
    free(obj);
}

int showPolygon(struct Polygon *obj) {
    struct Vertex *thisV = obj->P;
    struct Vertex *firstV = obj->P;
    if ((!obj->P)&&(!obj->num)) {// empty polygon
        printf("empty polygon\n\n");
        return 0;
    } else if (obj->P&&obj->num) {// normal non empty polygon
        printf("number:%4d\n", obj->num);
        for (int i=0;i<obj->num;i++) {
            printf("\taddr:%10x prev:%10x next:%10x x:%10.5lf y:%10.5lf flag:%5d\n",
                   (unsigned int)(thisV),
                   (unsigned int)(thisV->prev),
                   (unsigned int)(thisV->next),
                   thisV->x, thisV->y, thisV->flag);
            thisV = thisV->next;
        }
        printf("check:%4ld is 0\n\n", ((long)thisV) - ((long)firstV) );
    }
    return obj->num;
}

struct Vertex* intersect(struct Vertex *A, struct Vertex *B,
                         struct Vertex *C, struct Vertex *D) {
    double y_CD, x_AC, y_AC, x_CD, x_AB, y_AB;
    y_CD = D->y - C->y;
    y_AC = C->y - A->y;
    y_AB = B->y - A->y;
    x_AC = C->x - A->x;
    x_CD = D->x - C->x;
    x_AB = B->x - A->x;
    double denominator = y_CD*x_AB - x_CD*y_AB;
    if (fabs(denominator)<100*DBL_EPSILON) {
        printf("!!!!singular intersection!!!!");
        exit(1);
    }
    double m = (y_CD*x_AC - y_AC*x_CD)/denominator;
    double x = A->x + m * x_AB;
    double y = A->y + m * y_AB;
    struct Vertex *newV = (struct Vertex *)malloc(sizeof(struct Vertex));
    newV->x = x;
    newV->y = y;
    newV->flag = -1;
    return newV;
}

void CutPolygon(struct Polygon *poly,
                struct Vertex *D, struct Vertex *E, struct Vertex *F) {
    if (!poly->num) return;
    struct Vector3 *DE = newVector3(E->x-D->x, E->y-D->y, 0);
    struct Vector3 *DF = newVector3(F->x-D->x, F->y-D->y, 0);
    struct Vector3 *zD = cross3(DE, DF);
    struct Vector3 *vD = cross3(zD, DE);
    struct Vector3 *tempV = newVector3(0,0,0);
    struct Vertex  *thisV, *nextV, *newV, *prevV;
    struct Vertex  *in  = NULL;
    struct Vertex  *out = NULL;
    int countIn = 0;
    int countOut = 0;
    // update flag
    thisV = poly->P;
    for (int i=0; i<poly->num; i++) {
        setValue3(tempV, thisV->x - D->x, thisV->y - D->y, 0);
        if (dot3(tempV, vD)>0) {
            thisV->flag = 1;// >0 same side
            countIn++;
        } else {
            thisV->flag = 0;// <=0 not same side
            countOut++;
        }
        thisV = thisV->next;
    }
    //printf("after update flag, countIn:%d, countOut:%d\n", countIn, countOut);
    //showPolygon(poly);
    if ((countIn==0)||(countOut==0)) {// have all or reject all
        if (countIn==0) {// reject all
            thisV = poly->P;
            thisV->prev->next=NULL;
            while (thisV) {
                nextV = thisV->next;
                free(thisV);
                thisV = nextV;
            }
            poly->P   = 0;
            poly->num = 0;
        }
        // if countOut==0, do nothing, wait for clean and return
    } else { // have in and out
        // calculate intersect, update polygon
        thisV = poly->P;
        for (int i=0; i<poly->num; i++) {
            nextV = thisV->next;
            if (thisV->flag==0&&nextV->flag==1) { // ... to same  0 3 1
                newV = intersect(D, E, thisV, nextV);
                //printf("find new intersect:\n");
                //printf("\t(%7.3lf, %7.3lf) ===> (%7.3lf, %7.3lf)\n",
                //        thisV->x, thisV->y, nextV->x, nextV->y);
                //printf("\t\t(%7.3lf, %7.3lf)\n", newV->x, newV->y);
                newV->flag = 3;
                newV->prev = thisV;
                newV->next = nextV;
                nextV->prev = newV;
                thisV->next = newV;
                poly->num++;
                //showPolygon(poly);
                if (in) {
                    printf("find more than one in, can not be that!!");
                    exit(1);
                }
                in = newV;
                thisV = nextV;
            } else if (thisV->flag==1&&nextV->flag==0) { // same to ... 1 2 0
                newV = intersect(D, E, thisV, nextV);
                //printf("find new intersect:\n");
                //printf("\t(%7.3lf, %7.3lf) ===> (%7.3lf, %7.3lf)\n",
                //        thisV->x, thisV->y, nextV->x, nextV->y);
                //printf("\t\t(%7.3lf, %7.3lf)\n", newV->x, newV->y);
                newV->flag = 2;
                newV->prev = thisV;
                newV->next = nextV;
                nextV->prev = newV;
                thisV->next = newV;
                poly->num++;
                //showPolygon(poly);
                if (out) {
                    printf("find more than one out, can not be that!!");
                    exit(1);
                }
                out = newV;
                thisV = nextV;
            } else {
                thisV = thisV->next;
            }
        }
        if (((!in)&&(out))||((in)&&(!out))) {
            printf("only in or only out, can not be that!!");
            exit(1);
        }
        // free out boundary points
        //printf("do frees\n");
        thisV = out->next;
        nextV = thisV->next;
        while(thisV != in) {
            poly->num--;
            free(thisV);
            thisV = nextV;
            nextV = thisV->next;
        }
        out->next = in;
        in->prev = out;
        poly->P = in;
    }
    free(DE);free(DF);free(zD);free(vD);free(tempV);
}

// triangle cut by triangle
struct Polygon* TriangleTrianglePolygon(
        double Ax, double Ay, double Bx, double By, double Cx, double Cy,
        double Dx, double Dy, double Ex, double Ey, double Fx, double Fy) {
    struct Vertex *D = newVertex(Dx, Dy, 0);
    struct Vertex *E = newVertex(Ex, Ey, 0);
    struct Vertex *F = newVertex(Fx, Fy, 0);
    struct Polygon *initPolygon = newPolygon();
    addVertex(initPolygon, Ax, Ay, 0, 0);
    addVertex(initPolygon, Bx, By, 0, 1);
    addVertex(initPolygon, Cx, Cy, 0, 2);
    //printf("init triangle:\n");
    //showPolygon(initPolygon);
    //printf("first  cut:\n");
    CutPolygon(initPolygon, D, E, F);
    //showPolygon(initPolygon);
    //printf("second cut:\n");
    CutPolygon(initPolygon, E, F, D);
    //showPolygon(initPolygon);
    //printf("third  cut:\n");
    CutPolygon(initPolygon, F, D, E);
    //showPolygon(initPolygon);
    free(D);free(E);free(F);
    return initPolygon;
}

// triangle cut by square
struct Polygon* TriangleSquarePolygon(
        double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4,
        double Dx, double Dy, double Ex, double Ey, double Fx, double Fy) {
    struct Vertex *P1 = newVertex(x1, y1, 0);
    struct Vertex *P2 = newVertex(x2, y2, 0);
    struct Vertex *P3 = newVertex(x3, y3, 0);
    struct Vertex *P4 = newVertex(x4, y4, 0);
    struct Polygon *initPolygon = newPolygon();
    addVertex(initPolygon, Dx, Dy, 0, 0);
    addVertex(initPolygon, Ex, Ey, 0, 1);
    addVertex(initPolygon, Fx, Fy, 0, 2);
    CutPolygon(initPolygon, P1, P2, P3);
    CutPolygon(initPolygon, P2, P3, P4);
    CutPolygon(initPolygon, P3, P4, P1);
    CutPolygon(initPolygon, P4, P1, P2);
    free(P1);free(P2);free(P3);free(P4);
    return initPolygon;
}

double areaPolygon(struct Polygon *poly) {
    int num = poly->num;
    int i;
    struct Vertex *thisV = poly->P;
    if (!thisV) return 0.0;// None polygon
    double * xs = (double *)malloc(sizeof(double)*(num+1));
    double * ys = (double *)malloc(sizeof(double)*(num+1));
    double area=0;
    for (i=0;i<num;i++) {
        xs[i] = thisV->x;
        ys[i] = thisV->y;
        thisV = thisV->next;
    }
    xs[i] = thisV->x;
    ys[i] = thisV->y;
    for (i=0;i<num;i++) {
        area += xs[i]*ys[i+1] - xs[i+1]*ys[i];
    }
    area = fabs(area/2);
    free(xs);free(ys);
    return area;
}

// triangle cut by square
void areaTriangleSquare(
        double x1,double y1, double x2,double y2, double x3,double y3, double x4,double y4,
        double Dx, double Dy, double Ex, double Ey, double Fx, double Fy,
        double *fullArea, double *overlapArea, double *squareArea) {
    struct Vector3 *DE = newVector3(Ex-Dx, Ey-Dy, 0);
    struct Vector3 *DF = newVector3(Fx-Dx, Fy-Dy, 0);
    struct Vector3 *zD = cross3(DE, DF);
    *fullArea = 0.5*norm3(zD);
    struct Polygon* poly = TriangleSquarePolygon(
            x1, y1, x2, y2, x3, y3, x4, y4,
            Dx, Dy, Ex, Ey, Fx, Fy);
    *overlapArea = areaPolygon(poly);
    struct Polygon *square = newPolygon();
    addVertex(square, x1, y1, 0, 0);
    addVertex(square, x2, y2, 0, 1);
    addVertex(square, x3, y3, 0, 2);
    addVertex(square, x4, y4, 0, 3);
    *squareArea = areaPolygon(square);
    free(DE);free(DF);free(zD);
    freePolygon(square);
    freePolygon(poly);
}

void copyArray(int nx, int ny, double *from, double *to) {
    int i,j;
    for (j=0;j<ny;j++) {
        for (i=0;i<nx;i++) {
            to[nx*j + i] = from[nx*j + i];
        }
    }
}

void _copy(int nx, int ny, int i, int j, int direction,
          double *ULx, double *ULy,
          double *URx, double *URy,
          double *LLx, double *LLy,
          double *LRx, double *LRy) {
    int xgood, ygood;
    if (direction==0x11) { // UR
        xgood = i!=(nx-1);
        ygood = j!=(ny-1);
        if (xgood) {
            ULx[nx*j + i+1] = URx[nx*j + i];
            ULy[nx*j + i+1] = URy[nx*j + i];
        }
        if (ygood) {
            LRx[nx*(j+1) + i] = URx[nx*j + i];
            LRy[nx*(j+1) + i] = URy[nx*j + i];
        }
        if (xgood&&ygood) {
            LLx[nx*(j+1) + (i+1)] = URx[nx*j + i];
            LLy[nx*(j+1) + (i+1)] = URy[nx*j + i];
        }
        return;
    }
    if (direction==0x10) { // UL
        xgood = i!=(0);
        ygood = j!=(ny-1);
        if (xgood) {
            URx[nx*j + i-1] = ULx[nx*j + i];
            URy[nx*j + i-1] = ULy[nx*j + i];
        }
        if (ygood) {
            LLx[nx*(j+1) + i] = ULx[nx*j + i];
            LLy[nx*(j+1) + i] = ULy[nx*j + i];
        }
        if (xgood&&ygood) {
            LRx[nx*(j+1) + (i-1)] = ULx[nx*j + i];
            LRy[nx*(j+1) + (i-1)] = ULy[nx*j + i];
        }
        return;
    }
    if (direction==0x01) { // LR
        xgood = i!=(nx-1);
        ygood = j!=(0);
        if (xgood) {
            LLx[nx*j + i+1] = LRx[nx*j + i];
            LLy[nx*j + i+1] = LRy[nx*j + i];
        }
        if (ygood) {
            URx[nx*(j-1) + i] = LRx[nx*j + i];
            URy[nx*(j-1) + i] = LRy[nx*j + i];
        }
        if (xgood&&ygood) {
            ULx[nx*(j-1) + (i+1)] = LRx[nx*j + i];
            ULy[nx*(j-1) + (i+1)] = LRy[nx*j + i];
        }
        return;
    }
    if (direction==0x00) { // LL
        xgood = i!=(0);
        ygood = j!=(0);
        if (xgood) {
            LRx[nx*j + i-1] = LLx[nx*j + i];
            LRy[nx*j + i-1] = LLy[nx*j + i];
        }
        if (ygood) {
            ULx[nx*(j-1) + i] = LLx[nx*j + i];
            ULy[nx*(j-1) + i] = LLy[nx*j + i];
        }
        if (xgood&&ygood) {
            URx[nx*(j-1) + (i-1)] = LLx[nx*j + i];
            URy[nx*(j-1) + (i-1)] = LLy[nx*j + i];
        }
        return;
    }
}

int _find3(int nx, int ny, int i, int j, double *this, int direction) {
    if (direction==0x11) { // UR
        if (i==0||j==0) return 0;
        if (isfinite(this[nx*(j-1) + (i-1)])&&
            isfinite(this[nx*(j-1) + (i  )])&&
            isfinite(this[nx*(j  ) + (i-1)]))
            return 1;
        else
            return 0;
    }
    if (direction==0x10) { // UL
        if (i==nx-1||j==0) return 0;
        if (isfinite(this[nx*(j-1) + (i+1)])&&
            isfinite(this[nx*(j-1) + (i  )])&&
            isfinite(this[nx*(j  ) + (i+1)]))
            return 1;
        else
            return 0;
    }
    if (direction==0x01) { // LR
        if (i==0||j==ny-1) return 0;
        if (isfinite(this[nx*(j+1) + (i-1)])&&
            isfinite(this[nx*(j+1) + (i  )])&&
            isfinite(this[nx*(j  ) + (i-1)]))
            return 1;
        else
            return 0;
    }
    if (direction==0x00) { // LL
        if (i==nx-1||j==ny-1) return 0;
        if (isfinite(this[nx*(j+1) + (i+1)])&&
            isfinite(this[nx*(j+1) + (i  )])&&
            isfinite(this[nx*(j  ) + (i+1)]))
            return 1;
        else
            return 0;
    }
    return 0;
}

void _use3(int nx, int ny, int i, int j,
           double *thisx, double *thisy, int direction,
           double *ULx, double *ULy,
           double *URx, double *URy,
           double *LLx, double *LLy,
           double *LRx, double *LRy) {
    if (direction==0x11) { // UR
        thisx[nx*j + i] = -thisx[nx*(j-1) + (i-1)]
                       +thisx[nx*(j-1) + (i  )]
                       +thisx[nx*(j  ) + (i-1)];
        thisy[nx*j + i] = -thisy[nx*(j-1) + (i-1)]
                       +thisy[nx*(j-1) + (i  )]
                       +thisy[nx*(j  ) + (i-1)];
        _copy(nx, ny, i, j, direction, ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
        return;
    }
    if (direction==0x10) { // UL
        thisx[nx*j + i] = -thisx[nx*(j-1) + (i+1)]
                       +thisx[nx*(j-1) + (i  )]
                       +thisx[nx*(j  ) + (i+1)];
        thisy[nx*j + i] = -thisy[nx*(j-1) + (i+1)]
                       +thisy[nx*(j-1) + (i  )]
                       +thisy[nx*(j  ) + (i+1)];
        _copy(nx, ny, i, j, direction, ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
        return;
    }
    if (direction==0x01) { // LR
        thisx[nx*j +i] = -thisx[nx*(j+1) + (i-1)]
                      +thisx[nx*(j+1) + (i  )]
                      +thisx[nx*(j  ) + (i-1)];
        thisy[nx*j +i] = -thisy[nx*(j+1) + (i-1)]
                      +thisy[nx*(j+1) + (i  )]
                      +thisy[nx*(j  ) + (i-1)];
        _copy(nx, ny, i, j, direction, ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
        return;
    }
    if (direction==0x00) { // LL
        thisx[nx*j + i] = -thisx[nx*(j+1) + (i+1)]
                       +thisx[nx*(j+1) + (i  )]
                       +thisx[nx*(j  ) + (i+1)];
        thisy[nx*j + i] = -thisy[nx*(j+1) + (i+1)]
                       +thisy[nx*(j+1) + (i  )]
                       +thisy[nx*(j  ) + (i+1)];
        _copy(nx, ny, i, j, direction, ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
        return;
    }
}

int _find2x(int nx, int ny, int i, int j, double *this, int direction) {
    if (direction==0x11 || direction==0x01) { // UR or LR
        if (i==0||i==1) return 0;
        if (isfinite(this[nx*(j) + (i-1)])&&
            isfinite(this[nx*(j) + (i-2)]))
            return 1;
        else
            return 0;
    }
    if (direction==0x10 || direction==0x00) { // UL or LL
        if (i==nx-1||i==nx-2) return 0;
        if (isfinite(this[nx*(j) + (i+1)])&&
            isfinite(this[nx*(j) + (i+2)]))
            return 1;
        else
            return 0;
    }
    return 0;
}

void _use2x(int nx, int ny, int i, int j,
           double *thisx, double *thisy, int direction,
           double *ULx, double *ULy,
           double *URx, double *URy,
           double *LLx, double *LLy,
           double *LRx, double *LRy) {
    if (direction==0x11 || direction==0x01) { // UR or LR
        thisx[nx*j + i] = 2*thisx[nx*(j) + (i-1)] - thisx[nx*(j) + (i-2)];
        thisy[nx*j + i] = 2*thisy[nx*(j) + (i-1)] - thisy[nx*(j) + (i-2)];
        _copy(nx, ny, i, j, direction, ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
        return;
    }
    if (direction==0x10 || direction==0x00) { // UL or LL
        thisx[nx*j + i] = 2*thisx[nx*(j) + (i+1)] - thisx[nx*(j) + (i+2)];
        thisy[nx*j + i] = 2*thisy[nx*(j) + (i+1)] - thisy[nx*(j) + (i+2)];
        _copy(nx, ny, i, j, direction, ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
        return;
    }
}

int _find2y(int nx, int ny, int i, int j, double *this, int direction) {
    if (direction==0x11 || direction==0x10) { // UR or UL
        if (j==0||j==1) return 0;
        if (isfinite(this[nx*(j-1) + (i)])&&
            isfinite(this[nx*(j-2) + (i)]))
            return 1;
        else
            return 0;
    }
    if (direction==0x01 || direction==0x00) { // LR or LL
        if (j==ny-1||j==ny-2) return 0;
        if (isfinite(this[nx*(j+1) + (i)])&&
            isfinite(this[nx*(j+2) + (i)]))
            return 1;
        else
            return 0;
    }
    return 0;
}

void _use2y(int nx, int ny, int i, int j,
           double *thisx, double *thisy, int direction,
           double *ULx, double *ULy,
           double *URx, double *URy,
           double *LLx, double *LLy,
           double *LRx, double *LRy) {
    if (direction==0x11 || direction==0x10) { // UR or UL
        thisx[nx*j + i] = 2*thisx[nx*(j-1) + (i)] - thisx[nx*(j-2) + (i)];
        thisy[nx*j + i] = 2*thisy[nx*(j-1) + (i)] - thisy[nx*(j-2) + (i)];
        _copy(nx, ny, i, j, direction, ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
        return;
    }
    if (direction==0x01 || direction==0x00) { // LR or LL
        thisx[nx*j + i] = 2*thisx[nx*(j+1) + (i)] - thisx[nx*(j+2) + (i)];
        thisy[nx*j + i] = 2*thisy[nx*(j+1) + (i)] - thisy[nx*(j+2) + (i)];
        _copy(nx, ny, i, j, direction, ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
        return;
    }
}

// old and good
void get4map_backup(int nx, int ny,
             const long* mask,
             const double* xx, const double* yy,
             double *ULx, double *ULy,
             double *URx, double *URy,
             double *LLx, double *LLy,
             double *LRx, double *LRy) {
    int i,j,k;
    int totalN = nx*ny;
    double *Tx = (double *)malloc(sizeof(double)*totalN);
    double *Ty = (double *)malloc(sizeof(double)*totalN);

    // in latter comments, [1:3] = [1,2,3], CONTAINS the end point!!
    // calculate T
    //     j      i
    //Tx[0:ny-2,0:nx-2]=
    //              j       i        j+1      i        j        i+1       j+1     i+1
    //         (xx[0:ny-2,0:nx-2]+xx[1:ny-1,0:nx-2]+xx[0:ny-2,1:nx-1]+xx[1:ny-1,1:nx-1])
    //Ty is the same
    for (j=0;j<=ny-2;j++) {
        for (i=0;i<=nx-2;i++) {
            Tx[nx*j + i] = (xx[nx*j + i] + xx[nx*(j+1) + i] +
                            xx[nx*(j+1) + i+1] + xx[nx*(j) + i+1])/4;
            Ty[nx*j + i] = (yy[nx*j + i] + yy[nx*(j+1) + i] +
                            yy[nx*(j+1) + i+1] + yy[nx*(j) + i+1])/4;
        }
    }
    /*T[0:ny-2, 0:ny-2]
          ===> UR[0:ny-2, 0:nx-2]
          ===> UL[0:ny-2, 1:nx-1]
          ===> LR[1:ny-1, 0:nx-2]
          ===> LL[1:ny-1, 1:nx-1]*/

    /*
    leave UR[ny-1, 0:nx-2], UR[0:ny-2, nx-1] and UR[ny-1, nx-1]
          UR[ny-1, 0:nx-2] = 2*UR[ny-2, 0:nx-2] - UR[ny-3, 0:nx-2]
          UR[0:ny-2, nx-1] = 2*UR[0:ny-2, nx-2] - UR[0:ny-2, nx-3]
          UR[ny-1, nx-1]   = 3*UR[ny-2, nx-2] - UR[ny-3, nx-2] - UR[ny-2, nx-3]
    leave UL[ny-1, 1:nx-1], UL[0:ny-2, 0] and UL[ny-1, 0]
          UL[ny-1, 1:nx-1] = 2*UL[ny-2, 1:nx-1] - UL[ny-3, 1:nx-1]
          UL[0:ny-2, 0]    = 2*UL[0:ny-2, 1] - UL[0:ny-2, 2]
          UL[ny-1, 0]      = 3*UL[ny-2, 1] - UL[ny-3, 1] - UL[ny-2, 2]
    leave LR[0, 0:nx-2], LR[1:ny-1, nx-1] and LR[0, nx-1]
          LR[0, 0:nx-2]    = 2*LR[1, 0:nx-2] - LR[2, 0:nx-2]
          LR[1:ny-1, nx-1] = 2*LR[1:ny-1, nx-2] - LR[1:ny-1, nx-3]
          LR[0, nx-1]      = 3*LR[1, nx-2] - LR[2, nx-2] - LR[1, nx-3]
    leave LL[0, 1:nx-1], LL[1:ny-1, 0] and LL[0, 0]
          LL[0, 1:nx-1] = 2*LL[1, 1:nx-1] - LL[2, 1:nx-1]
          LL[1:ny-1, 0] = 2*LL[1:ny-1, 1] - LL[1:ny-1, 2]
          LL[0, 0]      = 3*LL[1,1] - LL[1,2] - LL[2,1]*/
    // doing
    /* doing main part
            j       i
          T[0:ny-2, 0:ny-2]
                            j       i
                    ===> UR[0:ny-2, 0:nx-2]
                            j       i+1
                    ===> UL[0:ny-2, 1:nx-1]
                            j+1     i
                    ===> LR[1:ny-1, 0:nx-2]
                            j+1     i+1
                    ===> LL[1:ny-1, 1:nx-1]*/
    for (j=0;j<=ny-2;j++) {
        for (i=0;i<=nx-2;i++) {
            URx[nx*j + i]       = Tx[nx*j + i];
            URy[nx*j + i]       = Ty[nx*j + i];
            ULx[nx*j + i+1]     = Tx[nx*j + i];
            ULy[nx*j + i+1]     = Ty[nx*j + i];
            LRx[nx*(j+1) + i] = Tx[nx*j + i];
            LRy[nx*(j+1) + i] = Ty[nx*j + i];
            LLx[nx*(j+1) + i+1] = Tx[nx*j + i];
            LLy[nx*(j+1) + i+1] = Ty[nx*j + i];
        }
    }
    // doing [0:ny-2, *]
    //    UR[0:ny-2, nx-1] = 2*UR[0:ny-2, nx-2] - UR[0:ny-2, nx-3]
    //    UL[0:ny-2, 0]    = 2*UL[0:ny-2, 1] - UL[0:ny-2, 2]
    for (j=0;j<=ny-2;j++) {
        URx[nx*j + nx-1] = 2*URx[nx*j + nx-2] - URx[nx*j + nx-3];
        URy[nx*j + nx-1] = 2*URy[nx*j + nx-2] - URy[nx*j + nx-3];
        ULx[nx*j + 0]    = 2*ULx[nx*j + 1] - ULx[nx*j + 2];
        ULy[nx*j + 0]    = 2*ULy[nx*j + 1] - ULy[nx*j + 2];
    }
    // doing [1:ny-1, *]
    //    LR[1:ny-1, nx-1] = 2*LR[1:ny-1, nx-2] - LR[1:ny-1, nx-3]
    //    LL[1:ny-1, 0] = 2*LL[1:ny-1, 1] - LL[1:ny-1, 2]
    for (j=1;j<=ny-1;j++) {
        LRx[nx*j + nx-1] = 2*LRx[nx*j + nx-2] - LRx[nx*j + nx-3];
        LRy[nx*j + nx-1] = 2*LRy[nx*j + nx-2] - LRy[nx*j + nx-3];
        LLx[nx*j + 0]    = 2*LLx[nx*j + 1] - LLx[nx*j + 2];
        LLy[nx*j + 0]    = 2*LLy[nx*j + 1] - LLy[nx*j + 2];
    }
    // doing [*, 0:nx-2]
    //    UR[ny-1, 0:nx-2] = 2*UR[ny-2, 0:nx-2] - UR[ny-3, 0:nx-2]
    //    LR[0, 0:nx-2]    = 2*LR[1, 0:nx-2] - LR[2, 0:nx-2]
    for (i=0;i<=nx-2;i++) {
        URx[nx*(ny-1) + i] = 2*URx[nx*(ny-2) + i] - URx[nx*(ny-3) + i];
        URy[nx*(ny-1) + i] = 2*URy[nx*(ny-2) + i] - URy[nx*(ny-3) + i];
        LRx[nx*(0) + i]    = 2*LRx[nx*(1) + i] - LRx[nx*(2) + i];
        LRy[nx*(0) + i]    = 2*LRy[nx*(1) + i] - LRy[nx*(2) + i];
    }
    // doing [*, 1:nx-1]
    //    UL[ny-1, 1:nx-1] = 2*UL[ny-2, 1:nx-1] - UL[ny-3, 1:nx-1]
    //    LL[0, 1:nx-1] = 2*LL[1, 1:nx-1] - LL[2, 1:nx-1]
    for (i=1;i<=nx-1;i++) {
        ULx[nx*(ny-1) + i] = 2*ULx[nx*(ny-2) + i] - ULx[nx*(ny-3) + i];
        ULy[nx*(ny-1) + i] = 2*ULy[nx*(ny-2) + i] - ULy[nx*(ny-3) + i];
        LLx[nx*0 + i] = 2*LLx[nx*1 + i] - LLx[nx*2 + i];
        LLy[nx*0 + i] = 2*LLy[nx*1 + i] - LLy[nx*2 + i];
    }
    // doing singles
    URx[nx*(ny-1)+nx-1]=3*URx[nx*(ny-2)+nx-2] - URx[nx*(ny-3)+nx-2] - URx[nx*(ny-2)+nx-3];
    URy[nx*(ny-1)+nx-1]=3*URy[nx*(ny-2)+nx-2] - URy[nx*(ny-3)+nx-2] - URy[nx*(ny-2)+nx-3];
    ULx[nx*(ny-1)+0]   =3*ULx[nx*(ny-2)+1]    - ULx[nx*(ny-3)+1]    - ULx[nx*(ny-2)+2];
    ULy[nx*(ny-1)+0]   =3*ULy[nx*(ny-2)+1]    - ULy[nx*(ny-3)+1]    - ULy[nx*(ny-2)+2];
    LRx[nx*0+nx-1]     = 3*LRx[nx*1+nx-2]     - LRx[nx*2+nx-2]      - LRx[nx*1+nx-3];
    LRy[nx*0+nx-1]     = 3*LRy[nx*1+nx-2]     - LRy[nx*2+nx-2]      - LRy[nx*1+nx-3];
    LLx[nx*0+0]        = 3*LLx[nx*1+1]        - LLx[nx*1+2]         - LLx[nx*2+1];
    LLy[nx*0+0]        = 3*LLy[nx*1+1]        - LLy[nx*1+2]         - LLy[nx*2+1];
    free(Tx); free(Ty);

    // fix all points with nan value but good mask
    double xbad, ybad;
    int modify, thisdirection;
    double *thisx, *thisy;
    while(1) {
        xbad = NAN; ybad = NAN; modify=0;
        // correct UR
        thisx = URx;
        thisy = URy;
        thisdirection = 0x11;
        for (j=0;j<=ny-1;j++) {
            for (i=0;i<=nx-1;i++) {
                if (!mask[nx*j + i]) continue;
                if (!isfinite(thisx[nx*j + i])) {
                    // try to use 3 point to fix it
                    if (_find3(nx, ny, i, j, thisx, thisdirection)) {
                        _use3(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 x points to fixt it
                    } else if (_find2x(nx, ny, i, j, thisx, thisdirection)) {
                        _use2x(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 y points to fixt it
                    } else if (_find2y(nx, ny, i, j, thisx, thisdirection)) {
                        _use2y(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // record this point as bad point
                    } else {
                        xbad = (float)i;
                        ybad = (float)j;
                    }
                }
            }
        }
        // correct UL
        thisx = ULx;
        thisy = ULy;
        thisdirection = 0x10;
        for (j=0;j<=ny-1;j++) {
            for (i=0;i<=nx-1;i++) {
                if (!mask[nx*j + i]) continue;
                if (!isfinite(thisx[nx*j + i])) {
                    // try to use 3 point to fix it
                    if (_find3(nx, ny, i, j, thisx, thisdirection)) {
                        _use3(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 x points to fixt it
                    } else if (_find2x(nx, ny, i, j, thisx, thisdirection)) {
                        _use2x(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 y points to fixt it
                    } else if (_find2y(nx, ny, i, j, thisx, thisdirection)) {
                        _use2y(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // record this point as bad point
                    } else {
                        xbad = (float)i;
                        ybad = (float)j;
                    }
                }
            }
        }
        // correct LR
        thisx = LRx;
        thisy = LRy;
        thisdirection = 0x01;
        for (j=0;j<=ny-1;j++) {
            for (i=0;i<=nx-1;i++) {
                if (!mask[nx*j + i]) continue;
                if (!isfinite(thisx[nx*j + i])) {
                    // try to use 3 point to fix it
                    if (_find3(nx, ny, i, j, thisx, thisdirection)) {
                        _use3(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 x points to fixt it
                    } else if (_find2x(nx, ny, i, j, thisx, thisdirection)) {
                        _use2x(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 y points to fixt it
                    } else if (_find2y(nx, ny, i, j, thisx, thisdirection)) {
                        _use2y(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // record this point as bad point
                    } else {
                        xbad = (float)i;
                        ybad = (float)j;
                    }
                }
            }
        }
        //correct LL
        thisx = LLx;
        thisy = LLy;
        thisdirection = 0x00;
        for (j=0;j<=ny-1;j++) {
            for (i=0;i<=nx-1;i++) {
                if (!mask[nx*j + i]) continue;
                if (!isfinite(thisx[nx*j + i])) {
                    // try to use 3 point to fix it
                    if (_find3(nx, ny, i, j, thisx, thisdirection)) {
                        _use3(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 x points to fixt it
                    } else if (_find2x(nx, ny, i, j, thisx, thisdirection)) {
                        _use2x(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 y points to fixt it
                    } else if (_find2y(nx, ny, i, j, thisx, thisdirection)) {
                        _use2y(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // record this point as bad point
                    } else {
                        xbad = (float)i;
                        ybad = (float)j;
                    }
                }
            }
        }
        if ((!modify)&&(isfinite(xbad)||isfinite(ybad))) {
            printf("have pixel that are not able to be auto corrected, exit\n");
            printf("\t(%5.0lf, %5.0lf)", xbad, ybad);
        } else if (!modify) {
            break;
        }
    }
    // mask value for bad pixels
    for (j=0;j<=ny-1;j++) {
        for (i=0;i<=nx-1;i++) {
            if (!mask[nx*j + i]) {
                ULx[nx*j + i] = NAN;
                ULy[nx*j + i] = NAN;
                URx[nx*j + i] = NAN;
                URy[nx*j + i] = NAN;
                LLx[nx*j + i] = NAN;
                LLy[nx*j + i] = NAN;
                LRx[nx*j + i] = NAN;
                LRy[nx*j + i] = NAN;
            }
        }
    }
    // confirm value for good pixels
    for (j=0;j<=ny-1;j++) {
        for (i=0;i<=nx-1;i++) {
            if (mask[nx*j + i]) {
                if (!isfinite(ULx[nx*j + i]))
                    printf("(%d, %d) miss ULx", i, j);
                if (!isfinite(ULy[nx*j + i]))
                    printf("(%d, %d) miss ULy", i, j);
                if (!isfinite(URx[nx*j + i]))
                    printf("(%d, %d) miss URx", i, j);
                if (!isfinite(URy[nx*j + i]))
                    printf("(%d, %d) miss URy", i, j);
                if (!isfinite(LLx[nx*j + i]))
                    printf("(%d, %d) miss LLx", i, j);
                if (!isfinite(LLy[nx*j + i]))
                    printf("(%d, %d) miss LLy", i, j);
                if (!isfinite(LRx[nx*j + i]))
                    printf("(%d, %d) miss LRx", i, j);
                if (!isfinite(LRy[nx*j + i]))
                    printf("(%d, %d) miss LRy", i, j);
            }
        }
    }
}

int get4map(int nx, int ny,
             const long* mask,
             const double* xx, const double* yy,
             double *ULx, double *ULy,
             double *URx, double *URy,
             double *LLx, double *LLy,
             double *LRx, double *LRy) {
    int i,j,k;
    int totalN = nx*ny;
    double *Tx = (double *)malloc(sizeof(double)*totalN);
    double *Ty = (double *)malloc(sizeof(double)*totalN);

    // in latter comments, [1:3] = [1,2,3], CONTAINS the end point!!
    // calculate T
    //     j      i
    //Tx[0:ny-2,0:nx-2]=
    //              j       i        j+1      i        j        i+1       j+1     i+1
    //         (xx[0:ny-2,0:nx-2]+xx[1:ny-1,0:nx-2]+xx[0:ny-2,1:nx-1]+xx[1:ny-1,1:nx-1])
    //Ty is the same
    for (j=0;j<=ny-2;j++) {
        for (i=0;i<=nx-2;i++) {
            Tx[nx*j + i] = (xx[nx*j + i] + xx[nx*(j+1) + i] +
                            xx[nx*(j+1) + i+1] + xx[nx*(j) + i+1])/4;
            Ty[nx*j + i] = (yy[nx*j + i] + yy[nx*(j+1) + i] +
                            yy[nx*(j+1) + i+1] + yy[nx*(j) + i+1])/4;
        }
    }
    /*T[0:ny-2, 0:ny-2]
          ===> UR[0:ny-2, 0:nx-2]
          ===> UL[0:ny-2, 1:nx-1]
          ===> LR[1:ny-1, 0:nx-2]
          ===> LL[1:ny-1, 1:nx-1] */

    /*
    leave UR[ny-1, 0:nx-2], UR[0:ny-2, nx-1] and UR[ny-1, nx-1]
          UR[ny-1, 0:nx-2] = 2*UR[ny-2, 0:nx-2] - UR[ny-3, 0:nx-2]
          UR[0:ny-2, nx-1] = 2*UR[0:ny-2, nx-2] - UR[0:ny-2, nx-3]
          UR[ny-1, nx-1]   = 3*UR[ny-2, nx-2] - UR[ny-3, nx-2] - UR[ny-2, nx-3]
    leave UL[ny-1, 1:nx-1], UL[0:ny-2, 0] and UL[ny-1, 0]
          UL[ny-1, 1:nx-1] = 2*UL[ny-2, 1:nx-1] - UL[ny-3, 1:nx-1]
          UL[0:ny-2, 0]    = 2*UL[0:ny-2, 1] - UL[0:ny-2, 2]
          UL[ny-1, 0]      = 3*UL[ny-2, 1] - UL[ny-3, 1] - UL[ny-2, 2]
    leave LR[0, 0:nx-2], LR[1:ny-1, nx-1] and LR[0, nx-1]
          LR[0, 0:nx-2]    = 2*LR[1, 0:nx-2] - LR[2, 0:nx-2]
          LR[1:ny-1, nx-1] = 2*LR[1:ny-1, nx-2] - LR[1:ny-1, nx-3]
          LR[0, nx-1]      = 3*LR[1, nx-2] - LR[2, nx-2] - LR[1, nx-3]
    leave LL[0, 1:nx-1], LL[1:ny-1, 0] and LL[0, 0]
          LL[0, 1:nx-1] = 2*LL[1, 1:nx-1] - LL[2, 1:nx-1]
          LL[1:ny-1, 0] = 2*LL[1:ny-1, 1] - LL[1:ny-1, 2]
          LL[0, 0]      = 3*LL[1,1] - LL[1,2] - LL[2,1] */
    // doing
    /* doing main part
            j       i
          T[0:ny-2, 0:ny-2]
                            j       i
                    ===> UR[0:ny-2, 0:nx-2]
                            j       i+1
                    ===> UL[0:ny-2, 1:nx-1]
                            j+1     i
                    ===> LR[1:ny-1, 0:nx-2]
                            j+1     i+1
                    ===> LL[1:ny-1, 1:nx-1] */
    for (j=0;j<=ny-2;j++) {
        for (i=0;i<=nx-2;i++) {
            URx[nx*j + i]       = Tx[nx*j + i];
            URy[nx*j + i]       = Ty[nx*j + i];
            ULx[nx*j + i+1]     = Tx[nx*j + i];
            ULy[nx*j + i+1]     = Ty[nx*j + i];
            LRx[nx*(j+1) + i] = Tx[nx*j + i];
            LRy[nx*(j+1) + i] = Ty[nx*j + i];
            LLx[nx*(j+1) + i+1] = Tx[nx*j + i];
            LLy[nx*(j+1) + i+1] = Ty[nx*j + i];
        }
    }
    // doing [0:ny-2, *]
    //    UR[0:ny-2, nx-1] = 2*UR[0:ny-2, nx-2] - UR[0:ny-2, nx-3]
    //    UL[0:ny-2, 0]    = 2*UL[0:ny-2, 1] - UL[0:ny-2, 2]
    for (j=0;j<=ny-2;j++) {
        URx[nx*j + nx-1] = 2*URx[nx*j + nx-2] - URx[nx*j + nx-3];
        URy[nx*j + nx-1] = 2*URy[nx*j + nx-2] - URy[nx*j + nx-3];
        ULx[nx*j + 0]    = 2*ULx[nx*j + 1] - ULx[nx*j + 2];
        ULy[nx*j + 0]    = 2*ULy[nx*j + 1] - ULy[nx*j + 2];
    }
    // doing [1:ny-1, *]
    //    LR[1:ny-1, nx-1] = 2*LR[1:ny-1, nx-2] - LR[1:ny-1, nx-3]
    //    LL[1:ny-1, 0] = 2*LL[1:ny-1, 1] - LL[1:ny-1, 2]
    for (j=1;j<=ny-1;j++) {
        LRx[nx*j + nx-1] = 2*LRx[nx*j + nx-2] - LRx[nx*j + nx-3];
        LRy[nx*j + nx-1] = 2*LRy[nx*j + nx-2] - LRy[nx*j + nx-3];
        LLx[nx*j + 0]    = 2*LLx[nx*j + 1] - LLx[nx*j + 2];
        LLy[nx*j + 0]    = 2*LLy[nx*j + 1] - LLy[nx*j + 2];
    }
    // doing [*, 0:nx-2]
    //    UR[ny-1, 0:nx-2] = 2*UR[ny-2, 0:nx-2] - UR[ny-3, 0:nx-2]
    //    LR[0, 0:nx-2]    = 2*LR[1, 0:nx-2] - LR[2, 0:nx-2]
    for (i=0;i<=nx-2;i++) {
        URx[nx*(ny-1) + i] = 2*URx[nx*(ny-2) + i] - URx[nx*(ny-3) + i];
        URy[nx*(ny-1) + i] = 2*URy[nx*(ny-2) + i] - URy[nx*(ny-3) + i];
        LRx[nx*(0) + i]    = 2*LRx[nx*(1) + i] - LRx[nx*(2) + i];
        LRy[nx*(0) + i]    = 2*LRy[nx*(1) + i] - LRy[nx*(2) + i];
    }
    // doing [*, 1:nx-1]
    //    UL[ny-1, 1:nx-1] = 2*UL[ny-2, 1:nx-1] - UL[ny-3, 1:nx-1]
    //    LL[0, 1:nx-1] = 2*LL[1, 1:nx-1] - LL[2, 1:nx-1]
    for (i=1;i<=nx-1;i++) {
        ULx[nx*(ny-1) + i] = 2*ULx[nx*(ny-2) + i] - ULx[nx*(ny-3) + i];
        ULy[nx*(ny-1) + i] = 2*ULy[nx*(ny-2) + i] - ULy[nx*(ny-3) + i];
        LLx[nx*0 + i] = 2*LLx[nx*1 + i] - LLx[nx*2 + i];
        LLy[nx*0 + i] = 2*LLy[nx*1 + i] - LLy[nx*2 + i];
    }
    // doing singles
    URx[nx*(ny-1)+nx-1]=3*URx[nx*(ny-2)+nx-2] - URx[nx*(ny-3)+nx-2] - URx[nx*(ny-2)+nx-3];
    URy[nx*(ny-1)+nx-1]=3*URy[nx*(ny-2)+nx-2] - URy[nx*(ny-3)+nx-2] - URy[nx*(ny-2)+nx-3];
    ULx[nx*(ny-1)+0]   =3*ULx[nx*(ny-2)+1]    - ULx[nx*(ny-3)+1]    - ULx[nx*(ny-2)+2];
    ULy[nx*(ny-1)+0]   =3*ULy[nx*(ny-2)+1]    - ULy[nx*(ny-3)+1]    - ULy[nx*(ny-2)+2];
    LRx[nx*0+nx-1]     = 3*LRx[nx*1+nx-2]     - LRx[nx*2+nx-2]      - LRx[nx*1+nx-3];
    LRy[nx*0+nx-1]     = 3*LRy[nx*1+nx-2]     - LRy[nx*2+nx-2]      - LRy[nx*1+nx-3];
    LLx[nx*0+0]        = 3*LLx[nx*1+1]        - LLx[nx*1+2]         - LLx[nx*2+1];
    LLy[nx*0+0]        = 3*LLy[nx*1+1]        - LLy[nx*1+2]         - LLy[nx*2+1];
    free(Tx); free(Ty);

    // fix all points with nan value but good mask
    double xbad, ybad;
    int modify, thisdirection;
    double *thisx, *thisy, *thisdx, *thisdy;
    // first use diagonal points
    while(1) {
        modify=0;
        // correct UR
        thisx = URx;
        thisy = URy;
        thisdx = LLx;
        thisdy = LLy;
        thisdirection = 0x11;
        for (j=0;j<=ny-1;j++) {
            for (i=0;i<=nx-1;i++) {
                if (!mask[nx*j + i]) continue;
                if (!isfinite(thisx[nx*j + i])) {
                    if (isfinite(thisdx[nx*j + i])) {
                        thisx[nx*j + i] = 2*xx[nx*j + i] - thisdx[nx*j + i];
                        thisy[nx*j + i] = 2*yy[nx*j + i] - thisdy[nx*j + i];
                        _copy(nx, ny, i, j, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify = 1;
                    }
                }
            }
        }
        // correct UL
        thisx = ULx;
        thisy = ULy;
        thisdx = LRx;
        thisdy = LRy;
        thisdirection = 0x10;
        for (j=0;j<=ny-1;j++) {
            for (i=0;i<=nx-1;i++) {
                if (!mask[nx*j + i]) continue;
                if (!isfinite(thisx[nx*j + i])) {
                    if (isfinite(thisdx[nx*j + i])) {
                        thisx[nx*j + i] = 2*xx[nx*j + i] - thisdx[nx*j + i];
                        thisy[nx*j + i] = 2*yy[nx*j + i] - thisdy[nx*j + i];
                        _copy(nx, ny, i, j, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify = 1;
                    }
                }
            }
        }
        // correct LL
        thisx = LLx;
        thisy = LLy;
        thisdx = URx;
        thisdy = URy;
        thisdirection = 0x00;
        for (j=0;j<=ny-1;j++) {
            for (i=0;i<=nx-1;i++) {
                if (!mask[nx*j + i]) continue;
                if (!isfinite(thisx[nx*j + i])) {
                    if (isfinite(thisdx[nx*j + i])) {
                        thisx[nx*j + i] = 2*xx[nx*j + i] - thisdx[nx*j + i];
                        thisy[nx*j + i] = 2*yy[nx*j + i] - thisdy[nx*j + i];
                        _copy(nx, ny, i, j, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify = 1;
                    }
                }
            }
        }
        // correct LR
        thisx = LRx;
        thisy = LRy;
        thisdx = ULx;
        thisdy = ULy;
        thisdirection = 0x01;
        for (j=0;j<=ny-1;j++) {
            for (i=0;i<=nx-1;i++) {
                if (!mask[nx*j + i]) continue;
                if (!isfinite(thisx[nx*j + i])) {
                    if (isfinite(thisdx[nx*j + i])) {
                        thisx[nx*j + i] = 2*xx[nx*j + i] - thisdx[nx*j + i];
                        thisy[nx*j + i] = 2*yy[nx*j + i] - thisdy[nx*j + i];
                        _copy(nx, ny, i, j, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify = 1;
                    }
                }
            }
        }
        if (!modify) {
            break;
        }
    }
    // then use 3 points or 2 points
    while(1) {
        xbad = NAN; ybad = NAN; modify=0;
        // correct UR
        thisx = URx;
        thisy = URy;
        thisdirection = 0x11;
        for (j=0;j<=ny-1;j++) {
            for (i=0;i<=nx-1;i++) {
                if (!mask[nx*j + i]) continue;
                if (!isfinite(thisx[nx*j + i])) {
                    // try to use 3 point to fix it
                    if (_find3(nx, ny, i, j, thisx, thisdirection)) {
                        _use3(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 x points to fixt it
                    } else if (_find2x(nx, ny, i, j, thisx, thisdirection)) {
                        _use2x(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 y points to fixt it
                    } else if (_find2y(nx, ny, i, j, thisx, thisdirection)) {
                        _use2y(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // record this point as bad point
                    } else {
                        xbad = (float)i;
                        ybad = (float)j;
                    }
                }
            }
        }
        // correct UL
        thisx = ULx;
        thisy = ULy;
        thisdirection = 0x10;
        for (j=0;j<=ny-1;j++) {
            for (i=0;i<=nx-1;i++) {
                if (!mask[nx*j + i]) continue;
                if (!isfinite(thisx[nx*j + i])) {
                    // try to use 3 point to fix it
                    if (_find3(nx, ny, i, j, thisx, thisdirection)) {
                        _use3(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 x points to fixt it
                    } else if (_find2x(nx, ny, i, j, thisx, thisdirection)) {
                        _use2x(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 y points to fixt it
                    } else if (_find2y(nx, ny, i, j, thisx, thisdirection)) {
                        _use2y(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // record this point as bad point
                    } else {
                        xbad = (float)i;
                        ybad = (float)j;
                    }
                }
            }
        }
        // correct LR
        thisx = LRx;
        thisy = LRy;
        thisdirection = 0x01;
        for (j=0;j<=ny-1;j++) {
            for (i=0;i<=nx-1;i++) {
                if (!mask[nx*j + i]) continue;
                if (!isfinite(thisx[nx*j + i])) {
                    // try to use 3 point to fix it
                    if (_find3(nx, ny, i, j, thisx, thisdirection)) {
                        _use3(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 x points to fixt it
                    } else if (_find2x(nx, ny, i, j, thisx, thisdirection)) {
                        _use2x(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 y points to fixt it
                    } else if (_find2y(nx, ny, i, j, thisx, thisdirection)) {
                        _use2y(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // record this point as bad point
                    } else {
                        xbad = (float)i;
                        ybad = (float)j;
                    }
                }
            }
        }
        //correct LL
        thisx = LLx;
        thisy = LLy;
        thisdirection = 0x00;
        for (j=0;j<=ny-1;j++) {
            for (i=0;i<=nx-1;i++) {
                if (!mask[nx*j + i]) continue;
                if (!isfinite(thisx[nx*j + i])) {
                    // try to use 3 point to fix it
                    if (_find3(nx, ny, i, j, thisx, thisdirection)) {
                        _use3(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 x points to fixt it
                    } else if (_find2x(nx, ny, i, j, thisx, thisdirection)) {
                        _use2x(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // try to use 2 y points to fixt it
                    } else if (_find2y(nx, ny, i, j, thisx, thisdirection)) {
                        _use2y(nx, ny, i, j, thisx, thisy, thisdirection,
                              ULx, ULy, URx, URy, LLx, LLy, LRx, LRy);
                        modify=1;
                    // record this point as bad point
                    } else {
                        xbad = (float)i;
                        ybad = (float)j;
                    }
                }
            }
        }
        if ((!modify)&&(isfinite(xbad)||isfinite(ybad))) {
            printf("have pixel that are not able to be auto corrected, exit\n");
            printf("\t(%5.0lf, %5.0lf)", xbad, ybad);
            return -1;
            break;
        } else if (!modify) {
            break;
        }
    }
    // mask value for bad pixels
    for (j=0;j<=ny-1;j++) {
        for (i=0;i<=nx-1;i++) {
            if (!mask[nx*j + i]) {
                ULx[nx*j + i] = NAN;
                ULy[nx*j + i] = NAN;
                URx[nx*j + i] = NAN;
                URy[nx*j + i] = NAN;
                LLx[nx*j + i] = NAN;
                LLy[nx*j + i] = NAN;
                LRx[nx*j + i] = NAN;
                LRy[nx*j + i] = NAN;
            }
        }
    }
    // confirm value for good pixels
    for (j=0;j<=ny-1;j++) {
        for (i=0;i<=nx-1;i++) {
            if (mask[nx*j + i]) {
                if (!isfinite(ULx[nx*j + i]))
                    printf("(%d, %d) miss ULx", i, j);
                if (!isfinite(ULy[nx*j + i]))
                    printf("(%d, %d) miss ULy", i, j);
                if (!isfinite(URx[nx*j + i]))
                    printf("(%d, %d) miss URx", i, j);
                if (!isfinite(URy[nx*j + i]))
                    printf("(%d, %d) miss URy", i, j);
                if (!isfinite(LLx[nx*j + i]))
                    printf("(%d, %d) miss LLx", i, j);
                if (!isfinite(LLy[nx*j + i]))
                    printf("(%d, %d) miss LLy", i, j);
                if (!isfinite(LRx[nx*j + i]))
                    printf("(%d, %d) miss LRx", i, j);
                if (!isfinite(LRy[nx*j + i]))
                    printf("(%d, %d) miss LRy", i, j);
            }
        }
    }
    return 1;
}

// search for first larger than  or equal to value
int searchsorted(int n, double* array, double value) {
    double this;
    int thissmallindex, thislargeindex, thisindex;
    if (value>array[n-1]) return n; // no one
    if (value<array[0]) return 0; // no one
    thissmallindex = 0;
    thislargeindex = n-1;
    while (1) {
        thisindex = (thissmallindex + thislargeindex)/2; this = array[thisindex];
        if (thisindex==thissmallindex) return thisindex+1;
        if (this>value) {
            thislargeindex = thisindex;
        } else if (this<value) {
            thissmallindex = thisindex;
        } else {
            return thisindex;
        }
    }
    return -1;
}

// assume that outxx, outyy is square grid!
int interpolate2D(int innx, int inny,
              const double* inxx,   const double* inyy,
              const double* influx, const double* invar,
              long   *inmask, double *inoutmask,
              double *inULx, double *inULy,
              double *inURx, double *inURy,
              double *inLLx, double *inLLy,
              double *inLRx, double *inLRy,
              int outnx, int outny,
              const double *outxx,  const double *outyy,
              double *outflux, double *outvar,
              double *outULx,  double *outULy,
              double *outURx,  double *outURy,
              double *outLLx,  double *outLLy,
              double *outLRx,  double *outLRy,
              double *outmask
              ) {
    int inTotal = innx*inny;
    int outTotal = outnx*outny;
    int i,j,oi,oj;
    int success;
    long *tempout = (long*)malloc(outTotal*sizeof(long));// outmask is all good..

    for (i=0;i<outTotal;i++) {// mask to get outmaps
        tempout[i] = 1;
    }
    // get map for in grid
    success = get4map(innx, inny, inmask, inxx, inyy,
            inULx, inULy, inURx, inURy,
            inLLx, inLLy, inLRx, inLRy);
    if (success<0) return success;
    // get map for out grid
    success = get4map(outnx, outny, tempout, outxx, outyy,
            outULx, outULy, outURx, outURy,
            outLLx, outLLy, outLRx, outLRy);
    if (success<0) return success;
    double xmin, xmax, ymin, ymax;
    int xMinIndex, xMaxIndex, yMinIndex, yMaxIndex;
    double *array = (double*)malloc(sizeof(double)*4);
    double *thisxUR, *thisxUL, *thisxLR, *thisxLL;
    double *thisyUR = (double*)malloc(sizeof(double)*outny);
    double *thisyUL = (double*)malloc(sizeof(double)*outny);
    double *thisyLR = (double*)malloc(sizeof(double)*outny);
    double *thisyLL = (double*)malloc(sizeof(double)*outny);
    double fullArea0, overlapArea0;
    double fullArea1, overlapArea1;
    double fullArea, overlapArea, factor3, factor4, squareArea;
    double dtemp;
    // fixed 1D outx and outy
    thisxUL = outULx;
    thisxUR = outURx;
    thisxLR = outLRx;
    thisxLL = outLLx;
    for (oj=0;oj<outny;oj++) {
        thisyUL[oj] = outULy[outnx*oj];
        thisyUR[oj] = outURy[outnx*oj];
        thisyLR[oj] = outLRy[outnx*oj];
        thisyLL[oj] = outLLy[outnx*oj];
    }
    // project each output pixel
    for (j=0;j<inny;j++) {
        for (i=0;i<innx;i++) {
            if (inmask[innx*j + i]) {
                //printf("do output (%d, %d), ", i, j);
                array[0] = inURx[innx*j + i];
                array[1] = inLRx[innx*j + i];
                array[2] = inULx[innx*j + i];
                array[3] = inLLx[innx*j + i];
                doubleMin(array, 4, &xmin);
                doubleMax(array, 4, &xmax);
                array[0] = inURy[innx*j + i];
                array[1] = inLRy[innx*j + i];
                array[2] = inULy[innx*j + i];
                array[3] = inLLy[innx*j + i];
                doubleMin(array, 4, &ymin);
                doubleMax(array, 4, &ymax);
                //printf("[%7.3lf, %7.3lf], [%7.3lf, %7.3lf] ==> ", xmin, xmax, ymin, ymax);

                // xmin
                array[0] = (double)searchsorted(outnx, thisxUL, xmin)-1;
                array[1] = (double)searchsorted(outnx, thisxUR, xmin)-1;
                array[2] = (double)searchsorted(outnx, thisxLR, xmin)-1;
                array[3] = (double)searchsorted(outnx, thisxLL, xmin)-1;
                doubleMin(array, 4, &dtemp);
                xMinIndex = (int)dtemp;
                // xmax
                array[0] = (double)searchsorted(outnx, thisxUL, xmax);
                array[1] = (double)searchsorted(outnx, thisxUR, xmax);
                array[2] = (double)searchsorted(outnx, thisxLR, xmax);
                array[3] = (double)searchsorted(outnx, thisxLL, xmax);
                doubleMax(array, 4, &dtemp);
                xMaxIndex = (int)dtemp;
                // ymin
                array[0] = (double)searchsorted(outny, thisyUL, ymin)-1;
                array[1] = (double)searchsorted(outny, thisyUR, ymin)-1;
                array[2] = (double)searchsorted(outny, thisyLR, ymin)-1;
                array[3] = (double)searchsorted(outny, thisyLL, ymin)-1;
                doubleMin(array, 4, &dtemp);
                yMinIndex = (int)dtemp;
                // ymax
                array[0] = (double)searchsorted(outny, thisyUL, ymax);
                array[1] = (double)searchsorted(outny, thisyUR, ymax);
                array[2] = (double)searchsorted(outny, thisyLR, ymax);
                array[3] = (double)searchsorted(outny, thisyLL, ymax);
                doubleMax(array, 4, &dtemp);
                yMaxIndex = (int)dtemp;
                //printf("[%d, %d], [%d, %d] ? %d\n",
                //        xMinIndex, xMaxIndex, yMinIndex, yMaxIndex,
                //        xMaxIndex>=0&&xMinIndex<=outnx-1&&yMaxIndex>=0&&yMinIndex<=outny
                //        );
                if (xMaxIndex>=0&&xMinIndex<=outnx-1&&yMaxIndex>=0&&yMinIndex<=outny) {
                    if (xMinIndex<0) xMinIndex=0;
                    if (yMinIndex<0) yMinIndex=0;
                    if (xMaxIndex>outnx-1) xMaxIndex=outnx-1;
                    if (yMaxIndex>outny-1) yMaxIndex=outny-1;
                    for (oj=yMinIndex;oj<=yMaxIndex;oj++){
                        for (oi=xMinIndex;oi<xMaxIndex;oi++){
                            areaTriangleSquare(
                                outULx[outnx*oj+oi], outULy[outnx*oj+oi],
                                outLLx[outnx*oj+oi], outLLy[outnx*oj+oi],
                                outLRx[outnx*oj+oi], outLRy[outnx*oj+oi],
                                outURx[outnx*oj+oi], outURy[outnx*oj+oi],
                                inULx[innx*j+i], inULy[innx*j+i],
                                inLLx[innx*j+i], inLLy[innx*j+i],
                                inURx[innx*j+i], inURy[innx*j+i],
                                &fullArea0,  &overlapArea0, &squareArea);
                            areaTriangleSquare(
                                outULx[outnx*oj+oi], outULy[outnx*oj+oi],
                                outLLx[outnx*oj+oi], outLLy[outnx*oj+oi],
                                outLRx[outnx*oj+oi], outLRy[outnx*oj+oi],
                                outURx[outnx*oj+oi], outURy[outnx*oj+oi],
                                inURx[innx*j+i], inURy[innx*j+i],
                                inLLx[innx*j+i], inLLy[innx*j+i],
                                inLRx[innx*j+i], inLRy[innx*j+i],
                                &fullArea1,  &overlapArea1, &squareArea);
                            fullArea    = fullArea0+fullArea1;
                            overlapArea = overlapArea0+overlapArea1;
                            factor3 = overlapArea/fullArea;
                            factor4 = overlapArea/squareArea;
                            /* // debug print
                            printf("\t(%10.5lf, %10.5lf), (%10.5lf, %10.5lf), (%10.5lf, %10.5lf), (%10.5lf, %10.5lf)\n",
                                    outULx[outnx*oj+oi], outULy[outnx*oj+oi],
                                    outLLx[outnx*oj+oi], outLLy[outnx*oj+oi],
                                    outLRx[outnx*oj+oi], outLRy[outnx*oj+oi],
                                    outURx[outnx*oj+oi], outURy[outnx*oj+oi]
                                    );
                            printf("\t(%10.5lf, %10.5lf), (%10.5lf, %10.5lf), (%10.5lf, %10.5lf)\n",
                                inULx[innx*j+i], inULy[innx*j+i],
                                inLLx[innx*j+i], inLLy[innx*j+i],
                                inURx[innx*j+i], inURy[innx*j+i]
                                );
                            printf("\t(%10.5lf, %10.5lf), (%10.5lf, %10.5lf), (%10.5lf, %10.5lf)\n",
                                inURx[innx*j+i], inURy[innx*j+i],
                                inLLx[innx*j+i], inLLy[innx*j+i],
                                inLRx[innx*j+i], inLRy[innx*j+i]
                                );
                            printf("\t\t%10.5lf/%10.5lf, %10.5lf/%10.5lf, %10.5lf/%10.5lf, %10.5lf/%10.5lf\n",
                                   overlapArea0, fullArea0,
                                   overlapArea1, fullArea1,
                                   overlapArea, fullArea,
                                   overlapArea, squareArea);*/
                            if (factor3>0) {
                                outflux[outnx*oj + oi] += factor3 * influx[innx*j + i];
                                // 0*inf = nan!!!!
                                outvar[outnx*oj + oi]  += factor3 * invar[innx*j + i];
                            }
                            outmask[outnx*oj + oi] += factor4;
                            inoutmask[innx*j + i] += factor3;
                        }
                    }
                }
            }
        }
    }
    free(array);free(tempout);
    free(thisyUR);free(thisyUL);free(thisyLR);free(thisyLL);
    return 0;
}


int get2map(int nx, const double* xx, const long* mask, double *Lx, double *Rx) {
    int i;
    for (i=0;i<=nx-2;i++) {
        Rx[i] = (xx[i]+xx[i+1])/2;
        Lx[i+1] = (xx[i]+xx[i+1])/2;
    }
    Rx[nx-1] = 2*Rx[nx-2] - Rx[nx-3];
    Lx[0] = 2*Lx[1] - Lx[2];
    // correct Lx
    for (i=0;i<=nx-1;i++) {
        if (!mask[i]) continue;
        if (!isfinite(Lx[i])) {
            if (i==nx-2||i==nx-1||(!isfinite(Lx[i+1])||!isfinite(Lx[i+2]))) {
                printf("Lx have pixel that are not able to be auto corrected(index = %d), exit\n",i);
                return -1;
            } else {
                Lx[i] = 2*Lx[i+1] - Lx[i+2];
            }
        }
    }
    // correct Rx
    for (i=0;i<=nx-1;i++) {
        if (!mask[i]) continue;
        if (!isfinite(Rx[i])) {
            if (i==0||i==1||(!isfinite(Lx[i-1])||!isfinite(Lx[i-2]))) {
                printf("Rx have pixel that are not able to be auto corrected(index = %d), exit\n",i);
                return -1;
            } else {
                Rx[i] = 2*Rx[i-1] - Rx[i-2];
            }
        }
    }
    // mask value for bad pixels
    for (i=0;i<=nx-1;i++) {
        if (!mask[i]) {
            Rx[i] = NAN;
            Lx[i] = NAN;
        }
    }
    // confirm value for good pixels
    for (i=0;i<=nx-1;i++) {
        if (mask[i]) {
            if (!isfinite(Rx[i]))
                printf("%d miss Rx", i);
            if (!isfinite(Lx[i]))
                printf("%d miss Lx", i);
        }
    }
    return 0;
}

int interpolate1D(int innx, const double* inxx,
                   const double* influx, const double* invar,
                   long   *inmask, double *inoutmask,
                   double *inLx, double *inRx,
                   int outnx, const double *outxx,
                   double *outflux, double *outvar,
                   double *outLx,  double *outRx,
                   double *outmask
                   ) {
    int i,j,oi,oj;
    int success;
    long *tempout = (long*)malloc(sizeof(long)*outnx);
    for (i=0;i<outnx;i++) {
        tempout[i] = 1;
    }
    success = get2map(innx, inxx, inmask, inLx, inRx);
    if (success<0) return success;
    success = get2map(outnx, outxx, tempout, outLx, outRx);
    if (success<0) return success;
    // do project
    double xmin, xmax;
    double overlapArea, inArea, outArea, infactor, outfactor;
    double outxmin, outxmax;
    double overlapmax, overlapmin;
    int a0, a1;
    int xMinIndex, xMaxIndex;
    for (i=0;i<innx;i++) {
        xmin = inLx[i]<inRx[i]?inLx[i]:inRx[i];
        xmax = inLx[i]>inRx[i]?inLx[i]:inRx[i];
        a0 = searchsorted(outnx, outLx, xmin)-1;
        a1 = searchsorted(outnx, outRx, xmin)-1;
        xMinIndex = a0<a1?a0:a1;
        a0 = searchsorted(outnx, outLx, xmax);
        a1 = searchsorted(outnx, outRx, xmax);
        xMaxIndex = a0>a1?a0:a1;
        //printf("do input %5d, (%7.3lf, %7.3lf) [%5d, %5d]\n",
        //                  i, xmin, xmax, xMinIndex, xMaxIndex);
        if (xMaxIndex>=0&&xMinIndex<=outnx-1) {
            if (xMinIndex<0) xMinIndex=0;
            if (xMaxIndex>outnx-1) xMaxIndex=outnx-1;
                for (oi=xMinIndex;oi<=xMaxIndex;oi++){
                    outxmin = outRx[oi]<outLx[oi]?outRx[oi]:outLx[oi];
                    outxmax = outRx[oi]>outLx[oi]?outRx[oi]:outLx[oi];
                    inArea = xmax-xmin;
                    outArea = xmax-xmin;
                    overlapmin = xmin>outxmin?xmin:outxmin;
                    overlapmax = xmax<outxmax?xmax:outxmax;
                    overlapArea = overlapmax-overlapmin>0?overlapmax-overlapmin:0;

                    infactor  = overlapArea/inArea;
                    outfactor = overlapArea/outArea;
                    outflux[oi] += infactor * influx[i];
                    outvar[oi]  += infactor * invar[i];
                    outmask[oi] += outfactor;
                    inoutmask[i] += infactor;
                }
        }
    }
    free(tempout);
    return 0;
}

// mask is bad mask
int fillHollowness(long* mask, int ny, int nx) {
    int i,j, modify=0;
    do {
        modify=0;
        for (j=1;j<ny-1;j++) {
            for (i=1;i<nx-1;i++) {
                if ((!mask[nx*j+i])&&
                    ((mask[nx*(j-1)+(i  )]&&mask[nx*(j  )+(i+1)]&&mask[nx*(j  )+(i-1)])||
                     (mask[nx*(j+1)+(i  )]&&mask[nx*(j  )+(i+1)]&&mask[nx*(j  )+(i-1)])||
                     (mask[nx*(j+1)+(i  )]&&mask[nx*(j-1)+(i  )]&&mask[nx*(j  )+(i-1)])||
                     (mask[nx*(j+1)+(i  )]&&mask[nx*(j-1)+(i  )]&&mask[nx*(j  )+(i+1)]))) {
                    mask[nx*j+i] = 1;
                    modify=1;
                }
            }
        }
    } while(modify==1);
    j=0;
        for (i=1;i<nx-1;i++) {
            if (!mask[nx*j+i]&&
                mask[nx*(j+1)+(i  )]&&
                mask[nx*(j  )+(i+1)]&&
                mask[nx*(j  )+(i-1)]) {
                mask[nx*j+i] = 1;
            }
        }
    j=ny-1;
        for (i=1;i<nx-1;i++) {
            if (!mask[nx*j+i]&&
                mask[nx*(j-1)+(i  )]&&
                mask[nx*(j  )+(i+1)]&&
                mask[nx*(j  )+(i-1)]) {
                mask[nx*j+i] = 1;
            }
        }
    for (j=1;j<ny-1;j++) {
        i=0;
            if (!mask[nx*j+i]&&
                mask[nx*(j+1)+(i  )]&&
                mask[nx*(j-1)+(i  )]&&
                mask[nx*(j  )+(i+1)]) {
                mask[nx*j+i] = 1;
            }
    }
    for (j=1;j<ny-1;j++) {
        i=nx-1;
            if (!mask[nx*j+i]&&
                mask[nx*(j+1)+(i  )]&&
                mask[nx*(j-1)+(i  )]&&
                mask[nx*(j  )+(i-1)]) {
                mask[nx*j+i] = 1;
            }
    }
    return 0;
}
