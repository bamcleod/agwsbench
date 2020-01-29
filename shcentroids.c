typedef int pixtype;

#include <stdio.h>
#include <time.h>
#include <stdlib.h>

/* Thresholded centroid */
void shthresh(pixtype *image, 
	      int *subaps, // Array of subaperture locations subap[2*i]=y subap[2*1+1]=x
	      int nsubaps,  // Number of subapertures
	      pixtype *flat,   // Flat image
	      pixtype *bias,   // Bias image
	      int ny,          // Number of pixels in image y-axis
	      int nx,          // Number of pixels in image x-axis
	      int nsubap,      // Number of pixels per subaperture
	      pixtype thresh,   // Threshold value
	      double *xcens,    // Return xcentoids
	      double *ycens     // Return xcentoids

)

{
    double ysum, xsum, sum;
    int isubap, iy, ix;
    int pixoff;
    double pixval;
    int count;
    int yoff, xoff;

    // For each subaperture
    for (isubap=0; isubap < nsubaps; isubap++) {
	count = 0;
	xsum = 0;
	ysum = 0;
	sum = 0;
	// For each pixel in the subaperture
	for (iy=0; iy < nsubap; iy++) {
	    for (ix=0; ix < nsubap; ix++) {
		yoff = subaps[2*isubap];
		xoff = subaps[2*isubap+1];
		pixoff = (yoff + iy) * nx + xoff + iy;

		pixval = (image[pixoff] - bias[pixoff]) / flat[pixoff];

		if (pixval > thresh) {
		    ysum += pixval * iy;
		    xsum += pixval * ix;
		    sum +=  pixval;
		    count++;
		}
	    }
	    if (count>0) {
		xcens[isubap] = xsum/count;
		ycens[isubap] = ysum/count;
	    }
	}
    }
}


int main (int argc, char *argv[]) {

    int imagesize = 1024;
    pixtype *image;
    pixtype *bias;
    pixtype *flat;
    int nsubapx = 48;
    int nsubapy = 48;
    int *subaps;
    int nsubaps = nsubapx * nsubapy;
    int subapsize = 8;
    int x0 = imagesize/2 - nsubapx/2 * subapsize;    
    int y0 = imagesize/2 - nsubapy/2 * subapsize;
    clock_t start, stop;
    int nloops = 1000;
    double *xcentroids, *ycentroids;
    int ix, iy, iloop;
    double thresh = 0;

    image = (pixtype *)calloc(imagesize*imagesize,sizeof(pixtype));
    bias = (pixtype *)calloc(imagesize*imagesize,sizeof(pixtype));
    flat = (pixtype *)calloc(imagesize*imagesize,sizeof(pixtype));
    xcentroids = (double *) malloc (nsubaps * sizeof(double));
    ycentroids = (double *) malloc (nsubaps * sizeof(double));

    subaps = (int *)malloc(nsubaps * 2 * sizeof(int));

    // Initialize the subaperture list    
    for (iy=0;iy<nsubapy;iy++) {
	for (ix=0;ix<nsubapx;ix++) {
	    subaps[iy*nsubapx+ix]   = iy * subapsize + y0;
	    subaps[iy*nsubapx+ix+1] = ix * subapsize + x0;
	}
    }

    // Initialize the input images
    for (iy=0; iy<imagesize; iy++) {
	for (ix=0; ix<imagesize; ix++) {
	    image[iy*imagesize+ix] = 1000;
	    bias[iy*imagesize+ix] = 0;
	    flat[iy*imagesize+ix] = 1.0;
	}
    }

    start = clock();
    for (iloop=0; iloop<nloops; iloop++) {
	shthresh(image, subaps, nsubaps, flat, bias, imagesize, imagesize, subapsize, thresh, xcentroids, ycentroids);
    }
    stop = clock();

    printf("%f msec per loop\n", (float)(stop-start) / CLOCKS_PER_SEC / nloops * 1000);

    return 0;
}


			
		    
		
	    
	
