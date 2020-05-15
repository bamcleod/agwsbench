#include "sh.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

pixtype debias[1024*1024] __attribute__((aligned(64)));
pixtype flattened[1024*1024] __attribute__((aligned(64)));

pixtype image[1024*1024] __attribute__((aligned(64)));
pixtype bias[1024*1024] __attribute__((aligned(64)));
pixtype flat[1024*1024] __attribute__((aligned(64)));


int main (int argc, char *argv[]) {

    int imagesize = 1024;

    int nsubapx = 48;
    int nsubapy = 48;
    int *subaps;
    int nsubaps = nsubapx * nsubapy;
    int subapsize = 8;
    int x0 = imagesize/2 - nsubapx/2 * subapsize;    
    int y0 = imagesize/2 - nsubapy/2 * subapsize;
    clock_t start, stop;
    int nloops = 1000;
    double *xcentroids, *ycentroids, *intensities;
    int ix, iy, iloop;
    double thresh = 0;
    int kernsize = 4;
    pixtype *kernel;
    int i;

    if (argc != 3) {
	fprintf(stderr, "Usage: shcentroids number_of_subaps pixels_per_subap\n");
	exit(1);
    }
    nsubapx = nsubapy = atoi(argv[1]);
    subapsize = atoi(argv[2]);

	    
    xcentroids = (double *) malloc (nsubaps * sizeof(double));
    ycentroids = (double *) malloc (nsubaps * sizeof(double));
    intensities = (double *) malloc (nsubaps * sizeof(double));

    subaps = (int *)malloc(nsubaps * 2 * sizeof(int));

    kernel = (pixtype *)calloc(kernsize * kernsize , sizeof(pixtype));

    // Initialize the subaperture list    
    for (iy=0;iy<nsubapy;iy++) {
	for (ix=0;ix<nsubapx;ix++) {
	    subaps[2*(iy*nsubapx+ix)]   = iy * subapsize + y0;
	    subaps[2*(iy*nsubapx+ix)+1] = ix * subapsize + x0;
	}
    }

    // Initialize the input images
    for (i=0; i<imagesize*imagesize; i++) {
	image[i] = 1000;
	bias[i] = 0;
	flat[i] = 1.0;
    }

    for (iloop=-1; iloop<nloops; iloop++) {
	if (iloop == 0) {
	    start = clock();
	}
    
	shthresh(subaps, nsubaps, imagesize, imagesize, subapsize, thresh, xcentroids, ycentroids);
    }
    stop = clock();

    printf("Thresholded centroid: %f msec per loop\n", (float)(stop-start) / CLOCKS_PER_SEC / nloops * 1000);

    for (iloop=-1; iloop<nloops; iloop++) {
	if (iloop == 0) {
	    start = clock();
	}
	shcorrelate_pointerstodata(subaps, nsubapx, nsubapy, imagesize, imagesize, subapsize, kernsize, kernsize, kernel, xcentroids, ycentroids);
    }
    stop = clock();
    printf("Cross correlation with pointers: %f msec per loop\n", (float)(stop-start) / CLOCKS_PER_SEC / nloops * 1000);

    for (iloop=-1; iloop<nloops; iloop++) {
	if (iloop == 0) {
	    start = clock();
	}
	shcorrelate_copydata(subaps, nsubapx, nsubapy, imagesize, imagesize, subapsize, kernsize, kernsize, kernel, xcentroids, ycentroids);
    }
    stop = clock();
    printf("Cross correlation with copied data: %f msec per loop\n", (float)(stop-start) / CLOCKS_PER_SEC / nloops * 1000);

    for (iloop=-1; iloop<nloops; iloop++) {
	if (iloop == 0) {
	    start = clock();
	}
	shcorrelate_fft(subaps, nsubapx, nsubapy, imagesize, imagesize, subapsize, kernsize, kernsize, kernel, xcentroids, ycentroids, intensities);
    }
    stop = clock();
    printf("Cross correlation with FFT: %f msec per loop\n", (float)(stop-start) / CLOCKS_PER_SEC / nloops * 1000);


    return 0;
}


