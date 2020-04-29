// THIS CODE HAS ONLY BEEN EVALUATED FOR TIMING AND NOT FOR CORRECTNESS OF RESULTS

#include "dfs.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include <complex.h>
#include <fftw3.h>


		
/* DFS computation of sum(abs(fft)) */

void dfs_computeffts(
		     pixtype *image, // The input image
		     double **bias,  // Bias frames for each subap
		     double **flat, // Flat frames for each subap (includes apodizing function)
		     int *subaps, // Array of subaperture locations subap[2*i]=y subap[2*1+1]=x
		     int nsubaps,  // Number of subapertures
		     int nx,          // Number of pixels in image x-axis
		     int ny,          // Number of pixels in image y-axis
		     int npixsubap,      // Number of pixels per subaperture
		     double **fftsums // Output ffts with dimension [nsubaps, npixsubap * (npixsubap + 1) /2]
		     )

{
    int ipix;
    int isubap, iy, ix;
    int pixoff;
    double pixval;
    int yoff, xoff;

    static pixtype ****datapointers = NULL;

    int iyap, ixap, iap;   // Subaperture indices
    static double *crosscorr = NULL;;
    double sum;
    double summax;
    int imax;
    double I1, I2, I3;
    int ixmax, iymax;
    double sumclock;
    clock_t start, stop;

    static double *flattened = NULL;
    static fftw_complex *fft = NULL;

    fftw_plan plan_forward;
    
    int npixfftin  = npixsubap * npixsubap;
    int npixfftout = npixsubap * (npixsubap + 1) /2;
    int indx;
    
    // Allocated buffers for the pointers

    if (flattened == NULL) {
	flattened = (double *) calloc (npixfftin,  sizeof(double));
	fft     = (fftw_complex *) malloc (npixfftout * sizeof(fftw_complex));

	// Create fftw plan
	plan_forward = fftw_plan_dft_r2c_2d(npixsubap, npixsubap, flattened, fft,       FFTW_MEASURE);
	
    }
    
    // For each subaperture
    for (isubap=0; isubap < nsubaps; isubap++) {
	// For each pixel in the subaperture get pixels from image, debias, and flatten
	for (iy=0; iy < npixsubap; iy++) {
	    for (ix=0; ix < npixsubap; ix++) {
		yoff = subaps[2*isubap];
		xoff = subaps[2*isubap+1];
		pixoff = (yoff + iy) * nx + xoff + ix;  // Offset into parent image
		ipix = iy * npixsubap+ ix;              // Offset into subaperture
		flattened[ipix] = (image[pixoff] - bias[isubap][ipix]) * flat[isubap][ipix];
	    }
	}

	// Compute FFT
	fftw_execute(plan_forward);

	// Compute sum of absolute value of FFT //
	for (int i=0; i < npixfftout; i++) {
	    fftsums[isubap][i] += sqrt(fft[i][0]*fft[i][0] + fft[i][1]*fft[i][1]);
	}
    }
}

			
		    
		
