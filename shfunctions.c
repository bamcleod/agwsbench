// THIS CODE HAS ONLY BEEN EVALUATED FOR TIMING AND NOT FOR CORRECTNESS OF RESULTS

#include "sh.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>

extern pixtype debias[1024*1024] __attribute__((aligned(64)));
extern pixtype flattened[1024*1024] __attribute__((aligned(64)));

extern pixtype image[1024*1024] __attribute__((aligned(64)));
extern pixtype bias[1024*1024] __attribute__((aligned(64)));
extern pixtype flat[1024*1024] __attribute__((aligned(64)));

/* Thresholded centroid */
void shthresh(int *subaps, // Array of subaperture locations subap[2*i]=y subap[2*1+1]=x
	      int nsubaps,  // Number of subapertures
	      int ny,          // Number of pixels in image y-axis
	      int nx,          // Number of pixels in image x-axis
	      int npixsubap,      // Number of pixels per subaperture
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
	for (iy=0; iy < npixsubap; iy++) {
	    for (ix=0; ix < npixsubap; ix++) {
		yoff = subaps[2*isubap];
		xoff = subaps[2*isubap+1];
		pixoff = (yoff + iy) * nx + xoff + ix;
		pixval = (image[pixoff] - bias[pixoff]) * flat[pixoff];

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

/* Cross correlation 
* This version makes multiple copies of each pixe, one for each lag computed but then results in well vectorized computation
* Seems to be slower than the previous version
*/
void shcorrelate_copydata(int *subaps, // Array of subaperture locations subap[2*i]=y subap[2*1+1]=x
		 int nxap,  // Number of subapertures x
		 int nyap,  // Number of subapertures y
		 int nx,          // Number of pixels in image x-axis
		 int ny,          // Number of pixels in image y-axis
		 int npixsubap,      // Number of pixels per subaperture
		 int nxkern,      // Number of pixels in the kernel y
		 int nykern,      // Number of pixels in the kernel x
		 pixtype *kernel, // Kernel
		 double *xcens,    // Return xcentoids
		 double *ycens     // Return xcentoids

)

{
    int ipix;
    int nsubaps = nxap * nyap;  // Total number of subapertures 
    int isubap, iy, ix;
    int pixoff;
    double pixval;
    int count;
    int yoff, xoff;

    int nkern = nxkern * nykern;
    static pixtype **datapointers = NULL;
    int nxlag = (npixsubap - nxkern + 1);
    int nylag = (npixsubap - nykern + 1);
    int nlags = nxlag * nylag;

    int iyap, ixap, iap;   // Subaperture indices
    int iylag, ixlag, ilag; // Lag indices
    int iykern, ixkern, ikern; // Pixel indices
    static double *crosscorr = NULL;;
    double sum;
    double summax;
    int imax;
    double I1, I2, I3;
    int ixmax, iymax;
    double sumclock;
    clock_t start, stop;
    int starty = ny/2 - npixsubap * nyap/2;
    int startx = nx/2 - npixsubap * nxap/2;
    int endy = starty + npixsubap * nyap;
    int endx = startx + npixsubap * nxap;
    static pixtype *buff = NULL;
    static pixtype *sums = NULL;
    static pixtype *mults = NULL;
    static pixtype *bigkern = NULL;
    
    // Allocated buffers for the pointers

    if (datapointers == NULL) {

	datapointers = (pixtype **)malloc(nsubaps * nlags * nkern * sizeof (pixtype*));
	crosscorr = (double *)malloc(nlags * sizeof(double));
	buff = (pixtype *)malloc(nsubaps * nlags * nkern * sizeof(pixtype));
	mults = (pixtype *)malloc(nsubaps * nlags * nkern * sizeof(pixtype));
	sums = (pixtype *)malloc(nsubaps * nlags * sizeof(pixtype));
	bigkern = (pixtype *)malloc(nsubaps * nlags * nkern * sizeof(pixtype));

	// Set up all the pointers

	// Loop over subapertures
	for (int i=0, iyap=0; iyap < nyap; iyap++) {
	    for (ixap=0; ixap < nxap; ixap++) {
		iap = iyap * nxap + ixap;
		// Loop over lags
		for (iylag=0; iylag < nylag; iylag++) {
		    for (ixlag=0; ixlag < nxlag; ixlag++) {
			ilag = iylag * nxlag + ixlag;
			// Loop over pixels in kernel
			for (iykern=0; iykern < nykern; iykern++) {
			    for (ixkern=0; ixkern < nxkern; ixkern++) {
				ikern = iykern * nxkern + ixkern;
				ix = subaps[2*iap]   + ixlag + ixkern;
				iy = subaps[2*iap+1] + iylag + iykern;
				datapointers[i] = & flattened[iy * nx + ix];
				bigkern[i++] = kernel[ikern];
			    }
			}
		    }
		}
	    }
	}
    }
    

    // Bias and flat correction
    
    for (iy=starty;iy<endy;iy++) {
	for (ix=startx;ix<endx;ix++) {
	    flattened[iy*nx+ix] = (image[iy*nx+ix] - bias[iy*nx+ix]) * flat[iy*nx+ix];
	}
    }

    // Copy the data into a linear buffer so it can be well vectorized
    for (int i=0; i < nsubaps * nlags * nkern; i++) {
	buff[i] = *(datapointers[i]);
    }

    // Compute all the cross correlations

    for (int i=0; i < nsubaps * nlags * nkern; i++) {
	mults[i] = buff[i] * bigkern[i];
    }
    for (int i=0; i<nsubaps * nlags; i++) {
	sums[i] = 0;
	for (ikern = 0; ikern < nkern; ikern++) {
	    sums[i] += mults[i*nkern+ikern];
	}
    }

    // For each subap find the peak in its cross correlation
    
    for (int i=0,iap=0; iap < nsubaps; iap++) {
	summax = 0;
	imax = -1;
	for (ilag = 0; ilag < nlags; ilag++) {
	    // For each pixel
	    sum = 0;
	    crosscorr[ilag] = sums[iap*nlags + ilag];
	    if (sum > summax) {
		imax = ilag;
		summax = sum;
	    }
	}

	//stop = clock();
	//sumclock += stop-start;

	ixmax = imax % nxkern;
	iymax = imax / nykern;

	if (ixmax == 0) { ixmax = 1; }
	if (iymax == 0) { iymax = 1; }
	if (ixmax == nxkern - 1) {ixmax--;}
	if (iymax == nykern - 1) {iymax--;}
	

	// Interpolate to find fractional pixel peak
	// y
	I1 = crosscorr[(iymax - 1) * nxkern + ixmax];
	I2 = crosscorr[(iymax    ) * nxkern + ixmax];
	I3 = crosscorr[(iymax + 1) * nxkern + ixmax];
	ycens[iap] = iymax + (I1 - I3) / (I1 + I3 - 2 * I2) / 2;
	// x
	I1 = crosscorr[iymax * nxkern + ixmax - 1];
	I2 = crosscorr[iymax * nxkern + ixmax    ];
	I3 = crosscorr[iymax * nxkern + ixmax + 1];
	xcens[iap] = ixmax + (I1 - I3) / (I1 + I3 - 2 * I2) / 2;

    }
    //  printf("Kernel calc: %f msec \n", (float)(sumclock) / CLOCKS_PER_SEC * 1000);

}

			
/* Cross correlation 
 * This version sets up pointers to the data and accesses each pixel by dereferencing the pointers.  This cannot be well vectorized but minimizes any copying of the data
*/
void shcorrelate_pointerstodata(int *subaps, // Array of subaperture locations subap[2*i]=y subap[2*1+1]=x
		 int nxap,  // Number of subapertures x
		 int nyap,  // Number of subapertures y
		 int nx,          // Number of pixels in image x-axis
		 int ny,          // Number of pixels in image y-axis
		 int npixsubap,      // Number of pixels per subaperture
		 int nxkern,      // Number of pixels in the kernel y
		 int nykern,      // Number of pixels in the kernel x
		 pixtype *kernel, // Kernel
		 double *xcens,    // Return xcentoids
		 double *ycens     // Return xcentoids

)

{
    int ipix;
    int nsubaps = nxap * nyap;  // Total number of subapertures 
    int isubap, iy, ix;
    int pixoff;
    double pixval;
    int count;
    int yoff, xoff;

    int nkern = nxkern * nykern;
    static pixtype ****datapointers = NULL;
    int nxlag = (npixsubap - nxkern + 1);
    int nylag = (npixsubap - nykern + 1);
    int nlags = nxlag * nylag;

    int iyap, ixap, iap;   // Subaperture indices
    int iylag, ixlag, ilag; // Lag indices
    int iykern, ixkern, ikern; // Pixel indices
    static double *crosscorr = NULL;;
    double sum;
    double summax;
    int imax;
    double I1, I2, I3;
    int ixmax, iymax;
    double sumclock;
    clock_t start, stop;
    int starty = ny/2 - npixsubap * nyap/2;
    int startx = nx/2 - npixsubap * nxap/2;
    int endy = starty + npixsubap * nyap;
    int endx = startx + npixsubap * nxap;

    
    // Allocated buffers for the pointers

    if (datapointers == NULL) {
	datapointers = (pixtype ****)malloc(nsubaps * sizeof (pixtype***));
	for (iap=0;iap<nsubaps;iap++) {
	    datapointers[iap] = (pixtype ***)malloc(nlags * sizeof (pixtype **));
	    for (ilag=0; ilag<nlags; ilag++) {
		datapointers[iap][ilag] = (pixtype **)malloc(nkern * sizeof (pixtype *));
	    }
	}

	// Set up all the pointers
	
	// Loop over subapertures
	for (iyap=0; iyap < nyap; iyap++) {
	    for (ixap=0; ixap < nxap; ixap++) {
		iap = iyap * nxap + ixap;
		// Loop over lags
		for (iylag=0; iylag < nylag; iylag++) {
		    for (ixlag=0; ixlag < nxlag; ixlag++) {
			ilag = iylag * nxlag + ixlag;
			// Loop over pixels in kernel
			for (iykern=0; iykern < nykern; iykern++) {
			    for (ixkern=0; ixkern < nxkern; ixkern++) {
				ikern = iykern * nxkern + ixkern;
				ix = subaps[2*iap]   + ixlag + ixkern;
				iy = subaps[2*iap+1] + iylag + iykern;
				datapointers[iap][ilag][ikern] = & flattened[iy * nx + ix];
			    }
			}
		    }
		}
	    }
	}
    }
    
    if (crosscorr == NULL) {
	crosscorr = (double *)malloc(nlags * sizeof(double));
    }

    // Bias and flat correction
    
    for (iy=starty;iy<endy;iy++) {
	for (ix=startx;ix<endx;ix++) {
	    flattened[iy*nx+ix] = (image[iy*nx+ix] - bias[iy*nx+ix]) * flat[iy*nx+ix];
	}
    }
    
    imax = 0;
    sumclock=0;
    // For each subaperture
    for (iap=0; iap < nsubaps; iap++) {


	/*
	yoff = subaps[2*iap];
	xoff = subaps[2*iap+1];
	for (iy=0; iy < npixsubap; iy++) {
	    for (ix=0; ix < npixsubap; ix++) {
		pixoff = (yoff + iy) * nx + xoff + ix;
		flattened[pixoff] = (image[pixoff] - bias[pixoff]) * flat[pixoff];
	    }
	}
	*/
	
	// For each lag compute one value in the cross correlation matrix

	//	start = clock();

	summax = 0;
	imax = -1;
	for (ilag = 0; ilag < nlags; ilag++) {
	    // For each pixel
	    sum = 0;
	    for (ikern = 0; ikern < nkern; ikern++) {
		sum += *(datapointers[iap][ilag][ikern]) * kernel[ikern];
	    }
	    crosscorr[ilag] = sum;
	    if (sum > summax) {
		imax = ilag;
		summax = sum;
	    }
	}

	//stop = clock();
	//sumclock += stop-start;

	ixmax = imax % nxkern;
	iymax = imax / nykern;

	if (ixmax == 0) { ixmax = 1; }
	if (iymax == 0) { iymax = 1; }
	if (ixmax == nxkern - 1) {ixmax--;}
	if (iymax == nykern - 1) {iymax--;}
	

	// Interpolate to find fractional pixel peak
	// y
	I1 = crosscorr[(iymax - 1) * nxkern + ixmax];
	I2 = crosscorr[(iymax    ) * nxkern + ixmax];
	I3 = crosscorr[(iymax + 1) * nxkern + ixmax];
	ycens[iap] = iymax + (I1 - I3) / (I1 + I3 - 2 * I2) / 2;
	// x
	I1 = crosscorr[iymax * nxkern + ixmax - 1];
	I2 = crosscorr[iymax * nxkern + ixmax    ];
	I3 = crosscorr[iymax * nxkern + ixmax + 1];
	xcens[iap] = ixmax + (I1 - I3) / (I1 + I3 - 2 * I2) / 2;

    }
    //  printf("Kernel calc: %f msec \n", (float)(sumclock) / CLOCKS_PER_SEC * 1000);

}

			
		    
