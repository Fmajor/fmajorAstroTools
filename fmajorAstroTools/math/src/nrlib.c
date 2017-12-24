#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "./nrlib.h"
#define NR_END 1
#define FREE_ARG char*

static char fileName[255];
static int blockNum, myDebug, finishGoodCount;
static double mindchi2=0.001, dlamda=10.0;
static int gIgnoreBadFit=0;
static double minA = -100, badAverageForReject=0.7;


// Utilities
/* Numerical Recipes standard error handler */
void nrerror(char error_text[]) {
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}
/* allocate a float vector with subscript range v[nl..nh] */
float *vector(long nl, long nh) {
    float *v;
    v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    if (!v) nrerror("allocation failure in vector()");
    return v-nl+NR_END;
}
/* allocate an int vector with subscript range v[nl..nh] */
int *ivector(long nl, long nh) {
    int *v;
    v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
    if (!v) nrerror("allocation failure in ivector()");
    return v-nl+NR_END;
}
/* allocate an unsigned char vector with subscript range v[nl..nh] */
unsigned char *cvector(long nl, long nh) {
    unsigned char *v;
    v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
    if (!v) nrerror("allocation failure in cvector()");
    return v-nl+NR_END;
}
/* allocate an unsigned long vector with subscript range v[nl..nh] */
unsigned long *lvector(long nl, long nh) {
    unsigned long *v;
    v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
    if (!v) nrerror("allocation failure in lvector()");
    return v-nl+NR_END;
}
/* allocate a double vector with subscript range v[nl..nh] */
double *dvector(long nl, long nh) {
    double *v;
    v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
    if (!v) nrerror("allocation failure in dvector()");
    return v-nl+NR_END;
}
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
float **matrix(long nrl, long nrh, long ncl, long nch) {
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    float **m;
    /* allocate pointers to rows */
    m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;
    /* allocate rows and set pointers to them */
    m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(long nrl, long nrh, long ncl, long nch) {
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    double **m;
    /* allocate pointers to rows */
    m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;
    /* allocate rows and set pointers to them */
    m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
int **imatrix(long nrl, long nrh, long ncl, long nch) {
    long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
    int **m;
    /* allocate pointers to rows */
    m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m += NR_END;
    m -= nrl;
    /* allocate rows and set pointers to them */
    m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
    if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}

/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
    long newrl, long newcl) {
    long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
    float **m;
    /* allocate array of pointers to rows */
    m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
    if (!m) nrerror("allocation failure in submatrix()");
    m += NR_END;
    m -= newrl;
    /* set pointers to rows */
    for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch) {
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
    float **m;
    /* allocate pointers to rows */
    m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
    if (!m) nrerror("allocation failure in convert_matrix()");
    m += NR_END;
    m -= nrl;
    /* set pointers to rows */
    m[nrl]=a-ncl;
    for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
    /* return pointer to array of pointers to rows */
    return m;
}
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh) {
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
    float ***t;
    /* allocate pointers to pointers to rows */
    t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
    if (!t) nrerror("allocation failure 1 in f3tensor()");
    t += NR_END;
    t -= nrl;
    /* allocate pointers to rows and set pointers to them */
    t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
    if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;
    /* allocate rows and set pointers to them */
    t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
    if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;
    for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
    for(i=nrl+1;i<=nrh;i++) {
        t[i]=t[i-1]+ncol;
        t[i][ncl]=t[i-1][ncl]+ncol*ndep;
        for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
    }
    /* return pointer to array of pointers to rows */
    return t;
}
/* free a float vector allocated with vector() */
void free_vector(float *v, long nl, long nh) {
    free((FREE_ARG) (v+nl-NR_END));
}
/* free an int vector allocated with ivector() */
void free_ivector(int *v, long nl, long nh) {
    free((FREE_ARG) (v+nl-NR_END));
}
/* free an unsigned char vector allocated with cvector() */
void free_cvector(unsigned char *v, long nl, long nh) {
    free((FREE_ARG) (v+nl-NR_END));
}
/* free an unsigned long vector allocated with lvector() */
void free_lvector(unsigned long *v, long nl, long nh) {
    free((FREE_ARG) (v+nl-NR_END));
}
/* free a double vector allocated with dvector() */
void free_dvector(double *v, long nl, long nh) {
    free((FREE_ARG) (v+nl-NR_END));
}
/* free a float matrix allocated by matrix() */
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch) {
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}
/* free a double matrix allocated by dmatrix() */
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch) {
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}
/* free an int matrix allocated by imatrix() */
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch) {
    free((FREE_ARG) (m[nrl]+ncl-NR_END));
    free((FREE_ARG) (m+nrl-NR_END));
}
/* free a submatrix allocated by submatrix() */
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch) {
    free((FREE_ARG) (b+nrl-NR_END));
}
/* free a matrix allocated by convert_matrix() */
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch) {
    free((FREE_ARG) (b+nrl-NR_END));
}
/* free a float f3tensor allocated by f3tensor() */
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
    long ndl, long ndh) {
    free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
    free((FREE_ARG) (t[nrl]+ncl-NR_END));
    free((FREE_ARG) (t+nrl-NR_END));
}
// recipes
/* Levenberg-Marquardt method */
//mrqminGB(x, y, sig,
//         xcen, xcenMin, xcenMax,
//         sigma, sigmaMi, sigmaMax,
//         a, ia, goodSum, nTrace, nPoly,
//         covar, alpha, gmod, bmod, &chisq, &alamda);
int mrqminGB(int ndata, const double x[], const double y[], const double sig[],
        const double xcenMin[], const double xcenMax[],
        const double sigmaMin[], const double sigmaMax[],
        double a[], int ia[], double goodSum[], int nTrace, int nPoly,
        double **covar, double **alpha, double gmod[], double bmod[], double reject,
        double *chisq, double *alamda) {
    int i,j,k,l,m, bad;
    FILE *fp;
    int ma=3*nTrace+nPoly;
    int ixcen, isigma, haveBadLine;
    double big, flag, temp, tempDa, tempMax, sign, tempMin;
    static int mfit, goodCount=0, badCount=0, finishCount=0, restart=0;
    static double ochisq,*atry,*beta,*da,**oneda, lastBestChi2;
    double polyStep = 5;
    double rejectFactor=2;
    double rejectFactor2=3;

   // first call
    if (*alamda < 0.0) {
        atry=dvector(1,ma);
        beta=dvector(1,ma);
        da=dvector(1,ma);
        for (mfit=0,j=1;j<=ma;j++)
            if (ia[j]) mfit++;
        oneda=dmatrix(1,mfit,1,1);
        *alamda=0.001;
        //a[nTrace*3+1] = 1;
        //for (j=nTrace*3+2; j<=ma; j++) a[j]=0;
        mrqcofGB(x,y,sig,ndata,a,ia,goodSum,nTrace,nPoly,alpha,beta,gmod,bmod,chisq);
        ochisq=(*chisq);
        lastBestChi2 = ochisq;
        for (j=1;j<=ma;j++) atry[j]=a[j];
            if (myDebug>2){
            sprintf(fileName,"%dx.dat", blockNum);
            fp=fopen(fileName,"wt");
            for(i=1; i<=ndata; i++) {
                fprintf(fp, "%20.15e ", x[i]);
            }
            fclose(fp);
            sprintf(fileName,"%dy.dat", blockNum);
            fp=fopen(fileName,"wt");
            for(i=1; i<=ndata; i++) {
                fprintf(fp, "%20.15e ", y[i]);
            }
            fclose(fp);
            sprintf(fileName,"%dsig.dat", blockNum);
            fp=fopen(fileName,"wt");
            for(i=1; i<=ndata; i++) {
                fprintf(fp, "%20.15e ", sig[i]);
            }
            fclose(fp);
            sprintf(fileName,"%dgmod.dat", blockNum);
            fp=fopen(fileName,"wt");
            for(i=1; i<=ndata; i++) {
                fprintf(fp, "%20.15e ", gmod[i]);
            }
            fclose(fp);
            sprintf(fileName,"%dbmod.dat", blockNum);
            fp=fopen(fileName,"wt");
            for(i=1; i<=ndata; i++) {
                fprintf(fp, "%20.15e ", bmod[i]);
            }
            fclose(fp);
            sprintf(fileName,"%da.dat", blockNum);
            fp=fopen(fileName,"wt");
            for(i=1; i<=nTrace*3; i++) {
                fprintf(fp, "%20.15e ", atry[i]);
            }
            fclose(fp);
            printf("goodSum:\n");
            for(i=1; i<=nTrace; i++)
                printf("%10g ", goodSum[i]);
            printf("\n");
            printf("alamda= %30.10g, chisq= %30.10g\n", *alamda, *chisq);
            }
    }
   // normal call
        if (myDebug>2) {
            printf("\tia:\n");
            for(i=1; i<=3*nTrace; i++) {
                printf("%10d ", ia[i]);
            }
            printf("\n");
        }
    for (j=0,l=1;l<=ma;l++) {
         if (ia[l]) {
             for (j++,k=0,m=1;m<=ma;m++) {
                 if (ia[m]) {
                    k++;
                    covar[j][k]=alpha[j][k];
                }
            }
            covar[j][j]=alpha[j][j]*(1.0+(*alamda));
            oneda[j][1]=beta[j];
        }
    }
    if (myDebug>3){
        printf("\n");
        printf("alpha :\n");
        sprintf(fileName,"%dalpha.dat", blockNum);
        fp=fopen(fileName,"wt");
        for(i=1; i<=mfit; i++) {
            for (j=1; j<=mfit; j++) {
                printf("%+6.4e ", alpha[i][j]);
                fprintf(fp, "%20.15e ", alpha[i][j]);
            }
            printf("\n");
            fprintf(fp,"\n");
        }
        fclose(fp);
        
        printf("covar :\n");
        sprintf(fileName,"%dcovar.dat", blockNum);
        fp=fopen(fileName,"wt");
        for(i=1; i<=mfit; i++) {
            for (j=1; j<=mfit; j++) {
                printf("%+6.4e ", covar[i][j]);
                fprintf(fp, "%20.15e ", covar[i][j]);
            }
            printf("\n");
            fprintf(fp,"\n");
        }
        printf("\n");
        fclose(fp);
    }
    bad = gaussj(covar,mfit,oneda,1);
    if (bad) {
        haveBadLine = 0;
        for (i=1; i<=nTrace; i++) {
            if (goodSum[i]<reject/rejectFactor && ia[(i-1)*3+2]) {// the ia[(i-1)*3+2] is important
                haveBadLine ++;
                a[(i-1)*3 + 1]  = -99;
                ia[(i-1)*3 + 2] = 0;//fix location, only
            }
        }
        for (int ind=1; ind<=ndata; ind++) {
            gmod[ind]=0;
            bmod[ind]=0;
        }
        if (myDebug>-1)
            printf("\tbad fit at (%10.5g, %10.5g) and drop %4d from %4d fibers, mfit=%4d\n", x[1], x[ndata], haveBadLine, nTrace, mfit);
        if (haveBadLine&&haveBadLine*3<mfit-nPoly) {
            free_dmatrix(oneda,1,mfit,1,1);
            free_dvector(da,1,ma);
            free_dvector(beta,1,ma);
            free_dvector(atry,1,ma);
            *alamda = -1;
            for (l=nTrace*3+1; l<=ma;l++) {
                a[l]=0;
            }
            if (myDebug>-1)
                printf("\tRemove bad lines and restart for this block\n");
            return 0;// restart!
        } else {
            for (int ind=1; ind<=nTrace; ind++) {
                goodSum[ind]=0;
            }
            for (l=nTrace*3+1; l<=ma;l++) {
                a[l]=0;
            }
            if (gIgnoreBadFit||haveBadLine*3==mfit-nPoly) {// all line masked.. normal bad fit
                return 3;// no last call
            } else {
                printf("not ignoreBadFit...\n");
                printf("\tnew a:\n");
                for(i=1; i<=mfit; i++) {
                    printf("%10.5g ", atry[i]);
                }
                printf("\n");
                printf("alamda= %30.10g, chisq= %30.10g\n", *alamda, *chisq);
                exit(0);
            }
        }
    }
    for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
   // last call
    if (*alamda == 0.0) {
         covsrt(covar,ma,ia,mfit);
         if (myDebug>3){
            printf("Last covar^-1:\n");
            for(i=1; i<=ma; i++) {
                for (j=1; j<=ma; j++) {
                    printf("%+6.4e ", covar[i][j]);
                }
                printf("\n");
            }
            printf("\n");
         }
         free_dmatrix(oneda,1,mfit,1,1);
         free_dvector(da,1,ma);
         free_dvector(beta,1,ma);
         free_dvector(atry,1,ma);
         return 2;
    }
    if (myDebug>2){
    printf("da:\n");
    for(i=1; i<=ma; i++) {
        printf("%10.5g ", da[i]);
    }
    printf("\n");
    }
    if (myDebug>3){
        printf("covar^-1:\n");
        for(i=1; i<=nTrace*3; i++) {
            for (j=1; j<=nTrace*3; j++) {
                printf("%+6.4e ", covar[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }
   // update atry based on constrains
    for (j=0, l=1, ixcen=1, isigma=1, flag=0; l<=nTrace*3;) {
        if (ia[l]) {
            tempDa = da[++j];
            tempDa = tempDa>5?5:tempDa;
            tempDa = tempDa<-5?-5:tempDa;
            atry[l]=a[l]+tempDa;
            if (atry[l]<minA)
                atry[l]=minA;// may cause singular error when a is too small
            l++;
        }// constrain for A
        if (ia[l]) { // constrain for xcen
                atry[l]=a[l]+da[++j];
            if (atry[l]>xcenMax[ixcen])
                atry[l]=xcenMax[ixcen];
            else
                if (atry[l]<xcenMin[ixcen])
                    atry[l]=xcenMin[ixcen];
        }
        if (atry[l]<atry[l-3]&&l>3) { // x2 before x1, swap and make flag
            flag=1;
            big = atry[l-3];
            atry[l-3] = atry[l]<xcenMin[ixcen-1]?xcenMin[ixcen-1]:atry[l];
            atry[l] = big>xcenMax[ixcen]?xcenMax[ixcen]:big;
        }
        ixcen++;l++;
        if (ia[l]) { // constrain for sigma
            atry[l]=a[l]+da[++j];
            if (atry[l]>sigmaMax[isigma])
                atry[l]=sigmaMax[isigma];
            else
                if (atry[l]<sigmaMin[isigma])
                    atry[l]=sigmaMin[isigma];
        } isigma++;l++;
    }
    l=nTrace*3+1;
    if (ia[l]){
        tempMax = dmax(y, ndata);
        tempMin = dmin(y, ndata);
        atry[l] = a[l] + da[++j];
        if (atry[l]>tempMax)
            atry[l] = tempMax;
        if (atry[l]<tempMin)
            atry[l] = tempMin;
    }
    for (l=nTrace*3+2; l<=ma;l++) { // change for poly
        if (ia[l]){
            atry[l]=a[l]+da[++j];
            sign = ((a[l]>0 && da[j]>0) || (a[l]<0 && da[j]<0))?1:-1;
            if (fabs(da[j])/(fabs(a[l])+1)>polyStep)
                atry[l]=a[l] * sign * polyStep;
        }
    }
    if (myDebug>2){
    printf("new a:\n");
    for(i=1; i<=ma; i++) {
        printf("%10.5g ", atry[i]);
    }
    printf("\n");
    }
   // build the next model based on atry
    mrqcofGB(x,y,sig,ndata,atry,ia,goodSum,nTrace,nPoly,covar,da,gmod,bmod,chisq);
    if (myDebug>1){
    printf("alamda= %30.10g, chisq= %30.10g\n", *alamda, *chisq);
    }
    // reject too bad ragions
    for (i=1, temp=0; i<=nTrace; i++) {
        temp += goodSum[i];
    }
    temp /= nTrace;
    
    if (temp<badAverageForReject) {
        badCount=99; // directly reject below..
        *chisq = 1e99;
        restart = 1;
    }
   // update ochisq
    if (*chisq < ochisq) {
        goodCount++;
        badCount=0;
        restart=0;
        if (goodCount>10) {
            *alamda *= 0.1;
            goodCount = 0;
        }
        lastBestChi2 = ochisq;
        ochisq=(*chisq);
        //printf("\n\n\n\nmindchi2: %20.10g\n", mindchi2);
        if (((lastBestChi2-ochisq)/lastBestChi2)<mindchi2) {
            finishCount++;
        } else {
            finishCount=0;
            //printf("\n\n\n\nmindchi2: %20.10g\n", mindchi2);
        }
        if (finishCount>finishGoodCount) {
            dlamda=10.0;
            goodCount = 0;
            badCount = 0;
            finishCount = 0;
            restart = 0;
            return 1;
        }
        for (j=0,l=1;l<=ma;l++) {
            if (ia[l]) {
                for (j++,k=0,m=1;m<=ma;m++) {
                    if (ia[m]) {
                        k++;
                        alpha[j][k]=covar[j][k];
                    }
                }
                beta[j]=da[j];
                a[l]=atry[l];
            }
        }
    } else {
        if (*chisq > ochisq) {
            *alamda *= dlamda;//10.0
            goodCount=0;
            finishCount=0;
            badCount++;
            if (badCount>20) {
                if (restart) {
                    haveBadLine = 0;
                    for (i=1; i<=nTrace; i++) {
                        if (goodSum[i]<reject/rejectFactor && ia[(i-1)*3+2]) {
                            haveBadLine += 1;
                            a[(i-1)*3 + 1]  = -99;
                            ia[(i-1)*3 + 2] = 0;
                        }
                    }
                    if (myDebug>-1)
                        printf("\tbad fit at (%10.5g, %10.5g) and drop %4d from %4d fibers, mfit=%4d\n", x[1], x[ndata], haveBadLine, nTrace, mfit);

                    if (haveBadLine&&haveBadLine*3<mfit-nPoly) {
                        free_dmatrix(oneda,1,mfit,1,1);
                        free_dvector(da,1,ma);
                        free_dvector(beta,1,ma);
                        free_dvector(atry,1,ma);
                        *alamda = -1;
                        for (l=nTrace*3+1; l<=ma;l++) {
                            a[l]=0;
                        }
                        if (myDebug>-1)
                            printf("\tRemove bad lines and restart for this block\n");
                        return 0;// restart!
                    } else {
                        for (int ind=1; ind<=nTrace; ind++) {
                            goodSum[ind]=0;
                        }
                        for (l=nTrace*3+1; l<=ma;l++) {
                            a[l]=0;
                        }
                        if (gIgnoreBadFit||haveBadLine*3==mfit-nPoly) {
                            return 3;// no last call
                        } else {
                            printf("not ignoreBadFit...\n");
                            printf("\tnew a:\n");
                            for(i=1; i<=nTrace*3; i++) {
                                printf("%10.5g ", atry[i]);
                            }
                            printf("\n");
                            printf("alamda= %30.10g, chisq= %30.10g\n", *alamda, *chisq);
                            exit(0);
                        }
                    }
                } else {
                    restart = 1;
                    badCount = 0;
                    *alamda = 0.001;
                    dlamda=2.33;
                }}
        *chisq=ochisq;
        } else {
            dlamda=10.0;
            goodCount = 0;
            badCount = 0;
            restart = 0;
            finishCount = 0;
            return 1;
        }
    }
        if (myDebug>1){
        sprintf(fileName,"%dx.dat", blockNum);
        fp=fopen(fileName,"wt");
        for(i=1; i<=ndata; i++) {
            fprintf(fp, "%20.15e ", x[i]);
        }
        fclose(fp);
        sprintf(fileName,"%dy.dat", blockNum);
        fp=fopen(fileName,"wt");
        for(i=1; i<=ndata; i++) {
            fprintf(fp, "%20.15e ", y[i]);
        }
        fclose(fp);
        sprintf(fileName,"%dsig.dat", blockNum);
        fp=fopen(fileName,"wt");
        for(i=1; i<=ndata; i++) {
            fprintf(fp, "%20.15e ", sig[i]);
        }
        fclose(fp);
        sprintf(fileName,"%dgmod.dat", blockNum);
        fp=fopen(fileName,"wt");
        for(i=1; i<=ndata; i++) {
            fprintf(fp, "%20.15e ", gmod[i]);
        }
        fclose(fp);
        sprintf(fileName,"%dbmod.dat", blockNum);
        fp=fopen(fileName,"wt");
        for(i=1; i<=ndata; i++) {
            fprintf(fp, "%20.15e ", bmod[i]);
        }
        fclose(fp);
        sprintf(fileName,"%da.dat", blockNum);
        fp=fopen(fileName,"wt");
        for(i=1; i<=nTrace*3; i++) {
            fprintf(fp, "%20.15e ", atry[i]);
        }
        fclose(fp);
        printf("goodSum:\n");
        for(i=1; i<=nTrace; i++)
            printf("%10g ", goodSum[i]);
        printf("goodCount: %3d, badCount: %3d, finishCount: %3d", goodCount, badCount, finishCount);
        printf("\n");
        printf("alamda= %30.10g, chisq= %30.10g\n", *alamda, *chisq);
        }
    
    
    /*for (i=1; i<=nTrace; i++) {
        if (goodSum[i]<reject/rejectFactor2 && ia[(i-1)*3+2]) {// the ia[(i-1)*3+2] is important
            a[(i-1)*3 + 1]  = -99;
            ia[(i-1)*3 + 2] = 0;//fix location, only
        }
    }*/
    
    return 0;
}
/* utilities function for Levenberg-Marquardt method */
// mrqcofGB(x,y,sig,ndata,a,ia,ma,alpha,beta,gmod,bmod,chisq);
void mrqcofGB(const double x[], const double y[], const double sig[], int ndata, double a[], int ia[], double goodSum[],
    int nTrace, int nPoly, double **alpha, double beta[], double gmodArray[], double bmodArray[], double *chisq) {
    int i,j,k,l,m,mfit,n;
    int ma=3*nTrace+nPoly;
    double ymod, gmod, bmod, wt,sig2i,dy,*dyda;
    double *eachGoodSum;

    dyda=dvector(1,ma);
    eachGoodSum = dvector(1,nTrace);
    for(i=1;i<=nTrace;i++) goodSum[i]=0;

    for (j=1, mfit=0;j<=ma;j++)
        if (ia[j]) mfit++;
    for (j=1;j<=mfit;j++) {
        for (k=1;k<=j;k++) alpha[j][k]=0.0;
        beta[j]=0.0;
    }
    *chisq=0.0;
    for (i=1;i<=ndata;i++) {
        //dGaussian(x[i], a, &gmod, dyda, eachGoodSum, nTrace*3, 7, 1);
        dGaussian(x[i], a, &gmod, dyda, eachGoodSum, nTrace*3, 1, 1);
        dPoly(x[i], x[1], x[ndata], a+3*nTrace+1, &bmod, dyda+3*nTrace+1, nPoly);
        ymod = gmod + bmod;
        gmodArray[i] = gmod;
        bmodArray[i] = bmod;
        sig2i=1.0/(sig[i]*sig[i]);
        if (sig2i>0)
            for(n=1; n<=nTrace; n++)
                goodSum[n] += eachGoodSum[n];
        dy=y[i]-ymod;
        for (j=0,l=1;l<=ma;l++) {
            if (ia[l]) {
                wt=dyda[l]*sig2i;
                for (j++,k=0,m=1;m<=l;m++)
                    if (ia[m]) alpha[j][++k] += wt*dyda[m];
                beta[j] += dy*wt;
            }
        }
        *chisq += dy*dy*sig2i;
    }
    for (j=2;j<=mfit;j++)
        for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
    free_dvector(dyda,1,ma);
    free_dvector(eachGoodSum,1,nTrace);
}

void lastmrqcofGB(const double x[], const double y[], const double sig[], int ndata, double a[], int ia[], double goodSum[],
              int nTrace, int nPoly, double gmodArray[], double bmodArray[], double *chisq) {
    int i,j,k,l,m,mfit,n;
    int ma=3*nTrace+nPoly;
    double ymod, gmod, bmod, wt,sig2i,dy,*dyda;
    double *eachGoodSum;
    
    dyda=dvector(1,ma);
    eachGoodSum = dvector(1,nTrace);
    for(i=1;i<=nTrace;i++) goodSum[i]=0;
    
    for (j=1, mfit=0;j<=ma;j++)
        if (ia[j]) mfit++;

    *chisq=0.0;
    for (i=1;i<=ndata;i++) {
        //dGaussian(x[i], a, &gmod, dyda, eachGoodSum, nTrace*3, 7, 1);
        dGaussian(x[i], a, &gmod, dyda, eachGoodSum, nTrace*3, 1, 1);
        dPoly(x[i], x[1], x[ndata], a+3*nTrace+1, &bmod, dyda+3*nTrace+1, nPoly);
        ymod = gmod + bmod;
        gmodArray[i] = gmod;
        bmodArray[i] = bmod;
        sig2i=1.0/(sig[i]*sig[i]);
        if (sig2i>0)
            for(n=1; n<=nTrace; n++)
                goodSum[n] += eachGoodSum[n];
        dy=y[i]-ymod;
        *chisq += dy*dy*sig2i;
    }

    free_dvector(dyda,1,ma);
    free_dvector(eachGoodSum,1,nTrace);
}


#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
/* expand in storage the covariance matrix covar, so as to take into account
 * parameters that are being held fixed. */
void covsrt(double **covar, int ma, int ia[], int mfit) {
    int i,j,k;
    double swap;
    for (i=mfit+1;i<=ma;i++)
        for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
    k=mfit;
    for (j=ma;j>=1;j--) {
        if (ia[j]) {
            for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
            for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
            k--;
        }
    }
}
#undef SWAP
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
/* Linear equation solution by Gauss-Jordan elimination */
int gaussj(double **a, int n, double **b, int m) {
    int *indxc,*indxr,*ipiv;
    int i,icol,irow,j,k,l,ll;
    double big,dum,pivinv,temp;
    indxc=ivector(1,n);
    indxr=ivector(1,n);
    ipiv=ivector(1,n);
    for (j=1;j<=n;j++) ipiv[j]=0;
    for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
            if (ipiv[j] != 1)
                for (k=1;k<=n;k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j][k]) >= big) {
                            big=fabs(a[j][k]);
                            irow=j;
                            icol=k;
                        }
                    } else if (ipiv[k] > 1) {
                        //nrerror("gaussj: Singular Matrix-1");
                        if (!gIgnoreBadFit)
                            printf("gaussj: Singular Matrix-1!\n");
                        return 1;
                    }
                }
        ++(ipiv[icol]);
        if (irow != icol) {
            for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
            for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
        }
        indxr[i]=irow;
        indxc[i]=icol;
        if (a[icol][icol] == 0.0) {
            //nrerror("gaussj: Singular Matrix-2");
            if (!gIgnoreBadFit)
                printf("gaussj: Singular Matrix-2\n");
            return 1;
        }
        pivinv=1.0/a[icol][icol];
        a[icol][icol]=1.0;
        for (l=1;l<=n;l++) a[icol][l] *= pivinv;
        for (l=1;l<=m;l++) b[icol][l] *= pivinv;
        for (ll=1;ll<=n;ll++)
            if (ll != icol) {
                dum=a[ll][icol];
                a[ll][icol]=0.0;
                for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
                for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
            }
    }
    for (l=n;l>=1;l--) {
        if (indxr[l] != indxc[l])
            for (k=1;k<=n;k++)
                SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }
    free_ivector(ipiv,1,n);
    free_ivector(indxr,1,n);
    free_ivector(indxc,1,n);
    return 0;
}
#undef SWAP
// function written by Jin
void findBlocks(int *xmin, int*xmax,
                // blocks[0] is nblocks, blocks[1:] is the end of each block
                int* blocks, int* bMin, int* bMax,
                int** corMat, int nTrace, int* nBands) {
    int nBlock, blockMin, blockMax, iFib, actFibCount;
    int theMax, finding, found, i, j;
    nBlock = 0;
    blockMin = xmin[0]; blockMax=xmax[0];
    // init blockMin with min(xmin)
    for (iFib=0; iFib<nTrace; iFib++)
        if (xmin[iFib]<blockMin)
            blockMin=xmin[iFib];
    for (iFib=0; iFib<nTrace; iFib++) {
        actFibCount = 0;
        theMax = xmax[iFib];
        if (theMax > blockMax) blockMax=theMax;
        for (finding=iFib+1, found=iFib+1;finding<nTrace;finding++){
            if (xmin[finding]<theMax) {
                actFibCount += 1;
                corMat[iFib][found  ] = finding;
                corMat[found++][iFib] = finding;
                if (xmax[finding]>blockMax)
                    blockMax = xmax[finding];
            }
        }
        if (actFibCount>0) { // in the same block
            corMat[iFib][iFib] = actFibCount;
        } else { // new block
            bMin[nBlock] = blockMin;
            bMax[nBlock] = blockMax;
            nBlock += 1;
            blocks[nBlock] = iFib;
            if (iFib+1<nTrace)
                blockMin = xmin[iFib+1];
            blockMax=xmax[iFib];
        }
    }
    blocks[0] = nBlock;
    for (i=0; i<nTrace; i++) {
        j = 0;
        while (corMat[i][j]==0&&j<i) j++;
        nBands[i] = i-j;
    }
}
void dGaussian(double x, double a[], double *y, double dyda[], double goodSum[], int tTrace, int subpixs, double subpixWidth) {
    int i,j,k;
    double base;
    double diff, denom, epow;
    double delta = subpixWidth/(subpixs-1);
    double xdelta, sigma, mu, expA, A, tempy;
    if (subpixs >= 3) {
        *y = 0.0;
        for (i=0; i<tTrace; i+=3){
            sigma = a[i+3];
            mu    = a[i+2];
            expA  = exp(a[i+1]);
            /*A     = a[i+1];*/
            denom = 1.0/sqrt(6.283185307) / sigma;
            tempy=0;
            for(xdelta = -subpixWidth/2, j=0; j<subpixs; j++, xdelta+=delta)  {
                diff  = (x + xdelta - mu) / sigma;
                epow  = 0.5*diff*diff;
                base = exp(-epow) * denom;
                tempy += expA * base;
                /*tempy += A * base;*/
            }
            *y   += tempy/subpixs;
            diff = (x - mu) / sigma;
            epow = 0.5*diff*diff;
            base = exp(-epow) * denom;
            dyda[i+1] = expA * base;
            dyda[i+2] = expA * base * diff / sigma;
            dyda[i+3] = expA * base * (epow*2.0-1)/sigma;
            /*dyda[i+1] = base;*/
            /*dyda[i+2] = A * base * (epow*2.0-1)/sigma;*/
            /*dyda[i+3] = A * base * diff / sigma;*/
            goodSum[(int)(i/3)+1] = base;
        }
    } else {
        *y = 0.0;
        for (i=0; i<tTrace; i+=3){
            sigma = a[i+3];
            mu    = a[i+2];
            expA  = exp(a[i+1]);
            denom = 1.0/sqrt(6.283185307) / sigma;
            diff = (x - mu) / sigma;
            epow = 0.5*diff*diff;
            base = exp(-epow) * denom;
            *y    += expA * base;
            dyda[i+1] = expA * base;
            dyda[i+2] = expA * base * diff / sigma;
            dyda[i+3] = expA * base * (epow*2.0-1)/sigma;
            goodSum[(int)(i/3)+1] = base;
        }
    }
}
void dPoly(double x, double xMin, double xMax, double a[], double *y, double dyda[], int n) {
    double p0=1;
    double xx = 2*(x-xMin)/(xMax-xMin) - 1;
    double p1, p2;
    int i;
    *y = 0;
    p1 = xx;
    *y += a[0] * p0 + a[1] * p1;
    dyda[0] = 1;
    dyda[1] = xx;
    for (i=1; i<=n-2; i++){
        p2 =  ((2*i+1)*xx*p1 - i*p0)/(i+1);
        *y += a[i+1] * p2;
        dyda[i+1] = p2;
        p1=p2;
        p0=p1;
    }
}
void findXLimits(int *xmin, int *xmax, const double *x, const double *xcen,
                 int nTrace, int nx, const double *sigma, double sigmal) {
    int i, t;
    double theMin, theMax;
    double diff;
    for(i=0; i< nTrace; i++) {
        if (xcen[i]<x[0]||xcen[i]>x[nx-1]) {
            printf("init xcen out of x range with xcen[%d]=%g but \
                    x[0]=%g and x[nx-1]=%g\n", i, xcen[i], x[0], x[nx-1]);
            exit(1);
        }
    }
    for(i=0; i< nTrace; i++) {
        diff = sigmal*sigma[i];
        theMin = xcen[i] - diff > x[0]   ? xcen[i] - diff : x[0];
        theMax = xcen[i] + diff < x[nx-1]? xcen[i] + diff : x[nx-1];
        for (t=0;    x[t]<theMin;)
            t++;
        xmin[i] = t--;
        for (t=nx-1; x[t]>theMax;)
            t--;
        xmax[i] = t++;
    }
}
/* void mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
 *   int ma, double **covar, double **alpha, double *chisq,
 *   void (*funcs)(double, double [], double *, double [], int), double *alamda) */
int gaussianFitIter(int ndata, const double x[], const double y[], const double sig[], const double logA[],
        const double xcen[], const double xcenMin[], const double xcenMax[],
        const double sigma[], const double sigmaMin[], const double sigmaMax[],
        int nTrace, int nPoly, int fullreject[], double reject, double sigmal, int ignoreBadFit,
        double a[], double p[], int ia[], int maxiter,
        double gmod[], double bmod[], double chi2[], int bMin[], int bMax[],
        int niters[], double myMindchi2, int myFinishGoodCount, double myDlamda, int debug) { // output
    int i,j,k;
    FILE *fp;
    int niter, finished;
    int nGood=0;
    int nCoeff=3;
    double **covar, **alpha;
    double alamda, chisq, *goodSum;
    // x value minmax,    xind minmax
    int *xmin, *xmax;//, *bMax, *bMin;
    int *blocks, **corMat, *nBands;
    int nBlock, ib;
    int blockMin, blockMax, blockWidth;
    int imax, imin, nFib;
    int *workia, *workfullreject, *fullworkia;
    const double *workx, *worky, *worksig, *workLogA;
    const double *workxcen, *workxcenMin, *workxcenMax;
    const double *worksigma, *worksigmaMin, *worksigmaMax;
    double *worka, *workgoodSum, *workp;
    double *workgmod, *workbmod;
    double lastGoodChi2, *lastGoodAns, *lastGoodP;
    //double *workPolyMin, *workPolyMax;
    //workPolyMin = dvector(1, nCoeff);
    //workPolyMax = dvector(1, nCoeff);

    int ansStart, ansEnd;
    int ma = nCoeff*nTrace + nPoly;
    int temp, temp1, temp2;
    double *tempy;
    fullworkia = ivector(1, ma);

        /*
        if (debug>0){
            printf("x: %x\t",           (unsigned) x);
            printf("y: %x\t",           (unsigned) y);
            printf("sig: %x\t",         (unsigned) sig);
            printf("logA: %x\t",        (unsigned) logA);
            printf("xcen: %x\t",        (unsigned) xcen);
            printf("xcenMin: %x\t",     (unsigned) xcenMin);
            printf("xcenMax: %x\t",     (unsigned) xcenMax);
            printf("sigma: %x\t",       (unsigned) sigma);
            printf("sigmaMin: %x\t",    (unsigned) sigmaMin);
            printf("sigmaMax: %x\t",    (unsigned) sigmaMax);
            printf("fullreject: %x\t",  (unsigned) fullreject);
            printf("gmod: %x\t",        (unsigned) gmod);
            printf("chi2: %x\t",        (unsigned) chi2);
            printf("bMin: %x\t",        (unsigned) bMin);
            printf("bMax: %x\t",        (unsigned) bMax);
            printf("niters: %x\t",      (unsigned) niters);
            printf("\n\n");
        }*/
    
    myDebug = debug;
    gIgnoreBadFit = ignoreBadFit;
    mindchi2 = myMindchi2;
    finishGoodCount = myFinishGoodCount;
    dlamda = myDlamda;
    
    lastGoodAns = dvector(1,ma);
    lastGoodP = dvector(1,ma);
    covar = dmatrix(1, ma, 1, ma);
    alpha = dmatrix(1, ma, 1, ma);
    //for (i=0;i<nTrace*nCoeff;i++) ia[i]=1;
    goodSum = dvector(1,nTrace);
    for (i=0;i<nTrace;i++) fullreject[i]=0;
    xmin = ivector(0, nTrace-1);
    xmax = ivector(0, nTrace-1);
    tempy = dvector(1, ndata);
   // find max and min of profile influence
    findXLimits(xmin, xmax, x, xcen, nTrace, ndata, sigma, sigmal);
        if (myDebug>0){
        fp=fopen("xminmax.dat","wt");
        for(i=0, k=0; i<nTrace; i++)
            fprintf(fp, "%5d %5d %5d %5d %5d\n", i, i*nCoeff, xmin[i], xmax[i], xmax[i]-xmin[i]+1);
        fclose(fp);
        fp=fopen("y.dat","wt");
        for(i=0; i<ndata; i++)
            fprintf(fp, "%20.15e ", y[i]);
        fclose(fp);
        fp=fopen("sig.dat","wt");
        for(i=0; i<ndata; i++)
            fprintf(fp, "%20.15e ", sig[i]);
        fclose(fp);
        }
   // find range of influence and blocks
    blocks = ivector(0, nTrace);
    nBands = ivector(0, nTrace-1);
    corMat = imatrix(0, nTrace-1, 0, nTrace-1);
    findBlocks(xmin, xmax, blocks, bMin, bMax, corMat, nTrace, nBands);
        if (myDebug>0){
        fp=fopen("corMat.dat","wt");
        for(i=0; i<nTrace; i++) {
            for (j=0; j<nTrace; j++) {
                fprintf(fp, "%2d ", corMat[i][j]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
        fp=fopen("blocks.dat","wt");
        temp = 0;
        temp1 = blocks[0];
        for (i=1; i<=temp1; i++) {
            fprintf(fp, "%4d %4d %4d %4d\n", temp, blocks[i], bMin[i-1], bMax[i-1]);
            printf("%4d %4d %4d %4d %4d\n", i, temp, blocks[i], bMin[i-1], bMax[i-1]);
            temp = blocks[i]+1;
        }
        printf("\n");
        fclose(fp);
        }
    nBlock = blocks[0];
    chi2[0] = nBlock;
    niters[0] = nBlock;
    // correct if xcenMin[?] < x[0] or xcenMax > x[ndata-1]
    for (i=0;i<nTrace;i++)
        if (xcenMin[i] < x[0]) {
            printf("xcenMin[%d]<x[0]: %10.5g < %10.5g", i, xcenMin[i], x[0]);
            exit(1);
        }
    for (i=0;i<nTrace;i++)
        if (xcenMax[i] > x[ndata-1]){
            printf("xcenMax[%d]>x[ndata-1]: %10.5g > %10.5g", i, xcenMax[i], x[ndata-1]);
            exit(1);
        }
    bMin--; bMax--; // [0~n-1] to [1~n]
    xcenMin--; xcenMax--; sigmaMin--; sigmaMax--; xcen--; sigma--; logA--;
    x--; y--; sig--; a--; ia--; p--; gmod--; bmod--; fullreject--;// the index of functions in numerical recipe are in range [1~n]
    imax = -1;
   // for each block
    for (ib=1; ib<=nBlock; ib++) {
            if (myDebug>1) {
            for (i=0; i<=nBlock+1; i++) {
                printf("%4d %4d %4d %4d\n", i, blocks[i], bMin[i], bMax[i]);
                temp = blocks[i]+1;
            }}
        blockNum = ib;
        nGood = 0;
        blockMin = bMin[ib]; // not the real index, but the offset
        blockMax = bMax[ib];
        blockWidth = blockMax - blockMin + 1;
        imin = imax + 1; // not the real index, but the offset
        imax = blocks[ib];
        nFib = imax - imin + 1;
        if (imax >=nTrace||imin<0||blockMin<0||blockMax>=ndata) {
            printf("something wrong\n nTrace:%d, ndata:%d, imin:%d,\
                    imax:%d, blockMin:%d, blockMax:%d",
                    nTrace, ndata, imin, imax, blockMin, blockMax);
            exit(1);
        }

        workx   = x   + blockMin;
        worky   = y   + blockMin;
        worksig = sig + blockMin;
        workLogA     = logA     + imin;
        workxcen     = xcen     + imin;
        workxcenMin  = xcenMin  + imin;
        workxcenMax  = xcenMax  + imin;
        worksigma    = sigma    + imin;
        worksigmaMin = sigmaMin + imin;
        worksigmaMax = sigmaMax + imin;
        // variables to modify
        workia  = ia + imin*nCoeff;
        workp   = p  + imin*nCoeff;
        worka   = a  + imin*nCoeff;
        for (i=1; i<=nFib*nCoeff;i++)
            fullworkia[i] = workia[i];
        for (i=nFib*nCoeff+1; i<=nFib*nCoeff+nPoly; i++)
            fullworkia[i] = 1;
        
        workgoodSum    = goodSum    + imin;
        workfullreject = fullreject + imin;
        workgmod = gmod + blockMin;
        workbmod = bmod + blockMin;
        if (myDebug>1) {
            printf("block index: %d\n", ib);
            printf("\t blockMin: %d, blockMax: %d, blockWidth: %d\n\
                   imin:%d, imax: %d, nFib: %d\n", blockMin, blockMax, blockWidth,
                   imin, imax, nFib);
        }
        
        for (i=0;i<nFib;i++) {
            worka[i*nCoeff+1] = workLogA[i+1];
            worka[i*nCoeff+2] = workxcen[i+1];
            worka[i*nCoeff+3] = worksigma[i+1];
        }
        for (i=1; i<=blockWidth; i++)
            tempy[i] = worky[i];
        // init value and range of constant term
        worka[nFib*nCoeff+1]=percentile(tempy, blockWidth, 0.5);
        // init value and range of higher order term
        for (i=2; i<=nPoly; i++) {
            worka[nFib*nCoeff+i]=0;
        }
        if (myDebug>1){
            printf("data and init a:\n");
            for (i=1; i<=blockWidth; i++) {
                printf("%10.5g ",workx[i]);
            }
            printf("\n");
            for (i=1; i<=blockWidth; i++) {
                printf("%10.5g ",worky[i]);
            }
            printf("\n");
            for (i=1; i<=blockWidth; i++) {
                printf("%10.5g ",worksig[i]);
            }
            printf("\n");
            for (i=1; i<=nFib*nCoeff+nPoly; i++) {
                printf("%10.5g ",worka[i]);
            }
            printf("\n");
            for (i=1; i<=nFib*nCoeff+nPoly; i++) {
                printf("%10d ",fullworkia[i]);
            }
            printf("\n");
        }
        // init covar
        for (i=1; i<=nFib*nCoeff+nPoly; i++) {
            for (j=1; j<=nFib*nCoeff+nPoly; j++) {
                covar[i][j]=0.0;
            }
        }
       // main iteration
        alamda = -1;
        for(niter=1; niter<maxiter; niter++) {
            finished=\
            mrqminGB(blockWidth, workx, worky, worksig,
                     workxcenMin, workxcenMax,
                     worksigmaMin, worksigmaMax,
                     worka, fullworkia, workgoodSum, nFib, nPoly,
                     covar, alpha, workgmod, workbmod, reject,
                     &chisq, &alamda);
            if (myDebug>1){
            printf("niter = %d\n", niter);
            printf("\t chi2 = %15.10g\nans = ", chisq);
            //}
            //if (myDebug>2){
            for (i=1; i<=nFib*nCoeff+nPoly; i++) {
                printf("%10.5g ", worka[i]);
            }
            printf("\n\n");
            }
            if (finished)
                break;
        }
       // backup before last call
        for (i=1; i<=nFib*nCoeff; i++) {
            workp[i] = covar[i][i];
            if (workp[i]<0) // bad fit..
                workp[i]=1e99;
        }
        lastGoodChi2 = chisq;
        for (i=1; i<=nFib*nCoeff+nPoly; i++) {
            lastGoodAns[i] = worka[i];
            lastGoodP[i] = workp[i];
        }
        if (myDebug>1) {
            printf("before last call chi2= %g\nans:\n", chisq);
            for (i=1; i<=nFib*nCoeff+nPoly;i++)
                printf("%10.5g ", worka[i]);
            printf("\n");
        }
        // last call
        if (finished<3) {
            if (myDebug>0) {
                printf("last call\n");
            }
            alamda = 0.0;
            finished=\
            mrqminGB(blockWidth, workx, worky, worksig,
                 workxcenMin, workxcenMax,
                 worksigmaMin, worksigmaMax,
                 worka, fullworkia, workgoodSum, nFib, nPoly,
                 covar, alpha, workgmod, workbmod, reject,
                 &chisq, &alamda);
        }
        if (myDebug>1) {
            printf("after last call chi2= %g \nans:\n", chisq);
            for (i=1; i<=nFib*nCoeff+nPoly;i++)
                printf("%10.5g ", worka[i]);
            printf("\n");
        }
        if (chisq<3*lastGoodChi2){
            for (i=1, j=1; i<=nFib; i++,j+=nCoeff) {
                if (workgoodSum[i]<reject || exp(worka[j])/covar[j][j]<0.01 || covar[j][j]<0) {
                    workfullreject[i] = 1;
                    worka[j]=-99;
                } else {
                    nGood++;
                }
            }
            for (i=1; i<=nFib*nCoeff; i++) {
                workp[i] = covar[i][i];
                if (workp[i]<0) // bad fit..
                    workp[i]=1e99;
            }
        } else {
            chisq = lastGoodChi2;
            for (i=1; i<=nFib*nCoeff+nPoly; i++) {// memory leak???
                worka[i] = lastGoodAns[i];
                workp[i] = lastGoodP[i];
            }
        }
        
        lastmrqcofGB(workx,worky,worksig,blockWidth,worka,fullworkia,workgoodSum,nFib,nPoly,workgmod,workbmod,&chisq);
        
        chi2[ib] = chisq/(blockWidth-nGood-1);
        niters[ib] = niter;
        if (myDebug>1) {
            printf("block: %d  chi2= %g reduce chi2= %g\nans:\n", ib, chisq, chisq/(blockWidth-nGood-1));
            for (i=1; i<=nFib*nCoeff+nPoly;i++)
                printf("%10.5g ", worka[i]);
            
            printf("\nans and p:\n");
            for (i=1; i<=nTrace*nCoeff;i++)
                printf("%4d %10.5g %10.5g\n", i, a[i], p[i]);
            printf("\n\n");
        }
    }
    bMin++; bMax++; // [0~n-1] to [1~n]
    free_ivector(fullworkia, 1, ma);
    free_dmatrix(covar, 1, ma, 1, ma);
    free_dmatrix(alpha, 1, ma, 1, ma);
    free_dvector(goodSum, 1, nTrace);
    free_dvector(tempy, 1, ndata);
    free_dvector(lastGoodAns, 1, ma);
    free_dvector(lastGoodP, 1, ma);
    //free_dvector(workPolyMin, 1, nCoeff);
    //free_dvector(workPolyMax, 1, nCoeff);
    free_ivector(xmin, 0, nTrace-1);
    free_ivector(xmax, 0, nTrace-1);
    free_ivector(blocks, 0, nTrace);
    //free_ivector(bMin, 0, nTrace-1);
    //free_ivector(bMax, 0, nTrace-1);
    free_ivector(nBands, 0, nTrace-1);
    free_imatrix(corMat, 0, nTrace-1, 0, nTrace-1);
    xcenMin++; xcenMax++; sigmaMin++; sigmaMax++; xcen++; sigma++;
    x++; y++; sig++; a++; ia++; p++; gmod++; bmod++; fullreject++; logA++;
    return 1;
}

double percentile(double A[], int n, double p) {
    // sub of A ~ [1,...n]
    int i;
    if (p<0 || p>1) {
        printf("p = %10.5g should be in between 0 and 1\n", p);
        exit(1);
    }
    i = (int) (1 + (n-1) * p);
    return randomized_select(A, 1, n, i);
}

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

double dmin(const double A[], int n) {
    // the sub of A is [1, n]
    double theMin=9e99;
    for (int i=1; i<=n; i++) {
        if (A[i]<theMin)
            theMin = A[i];
    }
    return theMin;
}

double dmax(const double A[], int n) {
    // the sub of A is [1, n]
    double theMax=-9e99;
    for (int i=1; i<=n; i++) {
        if (A[i]>theMax)
        theMax = A[i];
    }
    return theMax;
}


void genBlocks(int ndata, const double x[],
              const double xcen[], int nTrace, double sigma[],
              int blocks[], int bMin[], int bMax[], int nBnads[], int debug) { // output
    int i,j,k,l,m,n;
    FILE *fp;
    // x value minmax,    xind minmax
    int *xmin, *xmax;//, *bMax, *bMin;
    int **corMat, *nBands;
    int nBlock, ib;
    int blockMin, blockMax, blockWidth;
    int imax, imin, nFib;
    
    xmin = ivector(0, nTrace-1);
    xmax = ivector(0, nTrace-1);
    // find max and min of profile influence
    findXLimits(xmin, xmax, x, xcen, nTrace, ndata, sigma, 1);

    // find range of influence and blocks
    // blocks = ivector(0, nTrace-1);
    // nBands = ivector(0, nTrace-1);
    corMat = imatrix(0, nTrace-1, 0, nTrace-1);
    findBlocks(xmin, xmax, blocks, bMin, bMax, corMat, nTrace, nBands);
    free_imatrix(corMat, 0, nTrace-1, 0, nTrace-1);
    free_ivector(xmin, 0, nTrace-1);
    free_ivector(xmax, 0, nTrace-1);
}

void variateConvolve(int ndata, int npad,
                     const double x[], const double y[], const double sigma[],
                     double result[], double kernel[]) {
    int i,j,k,l;
    int nB = 2*npad + 1;
    double thisKerSum, thisResult, diff, value;
    
    double *xx, *yy;
    xx = (double *)malloc((2*npad + ndata) * sizeof(double));
    yy = (double *)malloc((2*npad + ndata) * sizeof(double));
    memset(xx, 0, (2*npad + ndata) * sizeof(double));
    memset(yy, 0, (2*npad + ndata) * sizeof(double));
    for (i=npad,j=0; j<ndata; i++, j++) {
        xx[i] = x[j];
        yy[i] = y[j];
    }
    value = x[1] - x[0];
    for (i=0;i<npad;i++)
        xx[i] = x[0] - (npad-i)*value;
    value = x[ndata-1] - x[ndata-2];
    for (i=npad+ndata,j=1;i<2*npad+ndata;i++,j++)
        xx[i] = x[ndata-1] + j*value;
    
    //double *kernel;
    //kernel = (double *)malloc(nB * ndata *sizeof(double));
    
    // build kernels for each pixel
    for (i=npad,j=0; j<ndata; i++, j++) {// i is ind for xx, j is ind for x
        thisKerSum = 0;
        if (sigma[j]>0) {
            for (k=-npad,l=0; k<=npad; k++,l++) {
                diff = (xx[i+k]-xx[i]);
                value= exp( - diff*diff/(sigma[j]*sigma[j] * 2.0));
                thisKerSum += value;
                kernel[nB*j+l] = value;
            }
            for (k=-npad,l=0; k<=npad; k++,l++) {
                kernel[nB*j+l] /= thisKerSum;
            }
        } else {
            for (k=-npad,l=0; k<=npad; k++,l++) {
                kernel[nB*j+l] = 0;
            }
            kernel[nB*j+npad] = 1;
        }
    }
    // do the convolve
    for (i=npad,j=0; j<ndata; i++, j++) {// i is ind for xx, j is ind for x
        thisResult = 0;
        for (k=-npad,l=0; k<=npad; k++,l++) {
            thisResult += kernel[nB*j+l] * yy[i+k];
        }
        result[j] = thisResult;
    }
    
    //free(kernel);
    free(xx);
    free(yy);
}

void variateDeconvolve(int ndata, int npad,
                     const double x[], const double y[], const double sigma[],
                     double result[], double kernel[],
                     int debug) {
    int i,j,k,l;
    int nB = 2*npad + 1;
    int nb = npad + 1;
    double thisKerSum, thisResult, diff, value;
    
    double *xx, *yy;
    xx = (double *)malloc((2*npad + ndata) * sizeof(double));
    yy = (double *)malloc((2*npad + ndata) * sizeof(double));
    memset(xx, 0, (2*npad + ndata) * sizeof(double));
    memset(yy, 0, (2*npad + ndata) * sizeof(double));
    for (i=npad,j=0; j<ndata; i++, j++) {
        xx[i] = x[j];
        yy[i] = y[j];
    }
    value = x[1] - x[0];
    for (i=0;i<npad;i++)
        xx[i] = x[0] - (npad-i)*value;
    value = x[ndata-1] - x[ndata-2];
    for (i=npad+ndata,j=1;i<2*npad+ndata;i++,j++)
        xx[i] = x[ndata-1] + j*value;
    
    // build kernels for each pixel
    //double *kernel;
    //kernel = (double *)malloc(nB * ndata *sizeof(double));
    for (i=npad,j=0; j<ndata; i++, j++) {// i is ind for xx, j is ind for x
        if (sigma[j]>0) {
            thisKerSum = 0;
            for (k=-npad,l=0; k<=npad; k++,l++) {
                diff = (xx[i+k]-xx[i]);
                value= exp( - diff*diff/(sigma[j]*sigma[j] * 2.0));
                thisKerSum += value;
                kernel[nB*j+l] = value;
            }
            for (k=-npad,l=0; k<=npad; k++,l++) {
                kernel[nB*j+l] /= thisKerSum;
            }
        } else {
            for (k=-npad,l=0; k<=npad; k++,l++) {
                kernel[nB*j+l] = 0;
            }
            kernel[nB*j+npad] = 1;
        }
    }

    // do LL decomposition for the kernel matrix
    double *L,*U;
    L = (double *)malloc(nb * ndata * sizeof(double));
    U = (double *)malloc(nb * ndata * sizeof(double));
    memset(L, 0, nb * ndata * sizeof(double));
    memset(U, 0, nb * ndata * sizeof(double));
        if (debug>1) {
        FILE *fp;
        fp=fopen("L.dat","wt");
        for(i=0; i<ndata; i++) {
            for (j=0; j<nb; j++)
                fprintf(fp, "%10.5e ", L[nb*i+j]);
            fprintf(fp, "\n");
        }
        fclose(fp);
        fp=fopen("U.dat","wt");
        for(i=0; i<ndata; i++) {
            for (j=0; j<nb; j++)
                fprintf(fp, "%10.5e ", U[nb*i+j]);
            fprintf(fp, "\n");
        }
        fclose(fp);
        }
    blockLUdecomposition(ndata, npad, kernel, L, U, debug);
        if (debug>0) {
        FILE *fp;
        fp=fopen("L.dat","wt");
        for(i=0; i<ndata; i++) {
            for (j=0; j<nb; j++)
                fprintf(fp, "%10.5e ", L[nb*i+j]);
            fprintf(fp, "\n");
        }
        fclose(fp);
        fp=fopen("U.dat","wt");
        for(i=0; i<ndata; i++) {
            for (j=0; j<nb; j++)
                fprintf(fp, "%10.5e ", U[nb*i+j]);
            fprintf(fp, "\n");
        }
        fclose(fp);
        }
    // solve for result use L and U
    solveUseLU(ndata, npad, L, U, y, result);
    free(xx);
    free(yy);
    free(L);
    free(U);
}

void solveUseLU(int n, int N, double L[], double U[], const double y[], double r[]) {
    int i,j,k,l;
    int nB = 2*N + 1;
    int nb = N + 1;
    double t;

    for (i=0; i<n; i++) {
        t = y[i];
        for (k=1; (k<=i)&&(k<=N); k++)
            t -= L[nb*i + k] * r[i-k];
        r[i] = t/L[nb*i + 0];
        if (!isfinite(r[i]))
            printf("bad value in deconvolve!!!\n");
    }
    for (i=n-1,j=0; j<n; i--, j++) {
        t = r[i];
        for (k=1; (k<=j)&&(k<=N); k++)
            t -= U[nb*i + k] * r[i+k];
        r[i] = t/U[nb*i + 0];
        if (!isfinite(r[i]))
            printf("bad value in deconvolve!!!\n");
    }
}

void blockLUdecomposition(int n, int N, double a[], double L[], double U[], int debug) {
    int i,j,k,l;
    int nB = 2*N + 1;
    int nb = N + 1;
    int aend, ai, ii;
    double t;
    
    // magic!!!!!! don't try to modify this code..
    for (i=0; i<n; i++)
        L[nb*i+0] = 1;
    
    for (i=0; i<n; i++) {
        aend = i+N<n?2*N:N+(n-i-1);
        for (ai=N, k=0; ai<=aend; ai++, k++) {
            t = a[nB*i + ai];
            for (j=1; (j+k<=N)&(i-j>=0); j++)
                t -= L[nb*i + j] * U[nb*(i-j) + (k+j)];
            U[nb*i + k] = t/L[nb*i + 0];
            if (!isfinite(U[nb*i + k])) {
                printf("bad value in deconvolve!!!\n");
            }
        }
        for (ii=1; (ii<=N)&&(i+ii<n); ii++) {
            t = a[nB*(i+ii) + (N-ii)];
            for (j=1; (ii+j<=N)&&(i-j>=0); j++)
                t -= L[nb*(i+ii) + (j+ii)] * U[nb*(i-j) + j];
            L[nb*(i+ii) + ii] = t/U[nb*i + 0];
            if (!isfinite(L[nb*(i+ii) + ii])) {
                printf("bad value in deconvolve!!!\n");
            }
        }
    }
}
