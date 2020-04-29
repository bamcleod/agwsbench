#include "dfs.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

pixtype flattened[NROWS*NCOLS] __attribute__((aligned(64)));
pixtype bias[NROWS*NCOLS] __attribute__((aligned(64)));
pixtype flat[NROWS*NCOLS] __attribute__((aligned(64)));


int main (int argc, char *argv[]) {

    int nsubapx = 6;
    int nsubapy = 3;
    int nx = NCOLS;
    int ny = NROWS;
    int *subaps;
    int nsubaps = nsubapx * nsubapy;
    int subapsize = 64;
    int npixsubap = subapsize * subapsize;
    clock_t start, stop;
    int nloops = 1000;
    double *xcentroids, *ycentroids;
    int ix, iy, iloop;
    double thresh = 0;
    int i;
    int subapspacing = 48;
    int x0=0, y0=0;

    double **flat;
    double **bias;
    double **fftsums;

    pixtype image[NROWS*NCOLS];

    // Create flat and bias arrays //
    flat    = (double **) malloc(sizeof(double *) * nsubaps);
    bias    = (double **) malloc(sizeof(double *) * nsubaps);
    fftsums = (double **) malloc(sizeof(double *) * nsubaps);

    for (i=0;i<nsubaps;i++) {
	flat[i]    = (double *) malloc (sizeof(double) * subapsize * subapsize);
	bias[i]    = (double *) malloc (sizeof(double) * subapsize * subapsize);
	fftsums[i] = (double *) malloc (sizeof(double) * subapsize * (subapsize + 1) / 2);
    }
	
    if (argc != 2) {
	fprintf(stderr, "Usage: dfs_timing pixels_per_subap\n");
	exit(1);
    }

    subapsize = atoi(argv[1]);
	    
    xcentroids = (double *) malloc (nsubaps * sizeof(double));
    ycentroids = (double *) malloc (nsubaps * sizeof(double));

    // Initialize the subaperture list    
    subaps = (int *)malloc(nsubaps * 2 * sizeof(int));
    for (iy=0;iy<nsubapy;iy++) {
        for (ix=0;ix<nsubapx;ix++) {
            subaps[2*(iy*nsubapx+ix)]   = iy * subapspacing + y0;
            subaps[2*(iy*nsubapx+ix)+1] = ix * subapspacing + x0;
        }
    }

    for (iloop=-1; iloop<nloops; iloop++) {
	if (iloop == 0) {
	    start = clock();
	}
	dfs_computeffts(image, bias, flat, subaps, nsubaps, nx, ny, subapsize, fftsums) ;
    }
    stop = clock();
    printf("FFT computation: %f msec per loop\n", (float)(stop-start) / CLOCKS_PER_SEC / nloops * 1000);


    return 0;
}


