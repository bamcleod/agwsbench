// THIS CODE HAS ONLY BEEN EVALUATED FOR TIMING AND NOT FOR CORRECTNESS OF RESULTS

#include "sh.h"

#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#include <complex.h>
#include <fftw3.h>


extern pixtype image[1024*1024] __attribute__((aligned(64)));
extern pixtype bias[1024*1024] __attribute__((aligned(64)));
extern pixtype flat[1024*1024] __attribute__((aligned(64)));
			
		
/* Cross correlation 
 * This version uses fftw
*/
void shcorrelate_fft(int *subaps, // Array of subaperture locations subap[2*i]=y subap[2*1+1]=x
		     int nxap,  // Number of subapertures x
		     int nyap,  // Number of subapertures y
		     int nx,          // Number of pixels in image x-axis
		     int ny,          // Number of pixels in image y-axis
		     int npixsubap,      // Number of pixels per subaperture
		     int nxkern,      // Number of pixels in the kernel y
		     int nykern,      // Number of pixels in the kernel x
		     pixtype *kernel, // Kernel
		     double *xcens,    // Return xcentoids
		     double *ycens,     // Return xcentoids
		     double *intensities
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
    static double *flattened = NULL;
    static fftw_complex *fftkern = NULL;
    static fftw_complex *fft = NULL;
    static fftw_complex *prod = NULL;
    fftw_plan plan_forward, plan_reverse;
    
    int npixfftin  = npixsubap * npixsubap;
    int npixfftout = npixsubap * (npixsubap/2 + 1);
    int indx;
    
    // Allocated buffers for the pointers

    if (flattened == NULL) {
	flattened = (double *) calloc (npixfftin,  sizeof(double));
	crosscorr  = (double *) malloc (npixfftin * sizeof(double));
	fftkern = (fftw_complex *) malloc (npixfftout * sizeof(fftw_complex));
	fft     = (fftw_complex *) malloc (npixfftout * sizeof(fftw_complex));
	prod    = (fftw_complex *) malloc (npixfftout * sizeof(fftw_complex));
				       
	// Create fftw plans
	 plan_forward = fftw_plan_dft_r2c_2d(npixsubap, npixsubap, flattened, fft,       FFTW_MEASURE);
	 plan_reverse = fftw_plan_dft_c2r_2d(npixsubap, npixsubap, prod,      crosscorr, FFTW_MEASURE);
	
	 // Compute conj(FFT(kernel))

	 // Permit input kernel to be smaller than npixsubap and then pad it into center of larger array
	 for (iykern=0; iykern < nykern; iykern++) {
	     for (ixkern=0; ixkern < nxkern; ixkern++) {
		 flattened[(iykern + (npixsubap-nykern)/2)*nxkern + ixkern + (npixsubap-nxkern)/2] = kernel[iykern * nxkern + ixkern];
	     }
	 }
	 fftw_execute(plan_forward);
	 for (int i=0; i < npixfftout; i++) {
	     fftkern[i] = conj(fft[i]);
	 }
    }

    // For each subaperture
    for (isubap=0; isubap < nsubaps; isubap++) {
	double edgesum = 0;
	double sum = 0;
	// For each pixel in the subaperture get pixels from image, debias, and flatten
	for (iy=0; iy < npixsubap; iy++) {
	    for (ix=0; ix < npixsubap; ix++) {
		yoff = subaps[2*isubap];
		xoff = subaps[2*isubap+1];
		pixoff = (yoff + iy) * nx + xoff + ix;
		flattened[iy*npixsubap+ix] = (image[pixoff] - bias[pixoff]) * flat[pixoff];
		if (iy==0 || iy==npixsubap-1 || ix==0 || ix==npixsubap-1) {
		    edgesum +=1;
		} else {
		    sum += 1;
		}
	    }
	    /* Compute local background using mean of edge pixels */
	    int np1 = npixsubap - 1;
	    sum -= edgesum * np1 * np1 / (4 * np1);  /* Scale background sum by (pixels in box) / (pixels in edge) */
	}

	// Compute FFT
	fftw_execute(plan_forward);

	// Multiply by conj(FFT(kernel))
	for (int i=0; i < npixfftout; i++) {
	    prod[i] = fft[i] * fftkern[i];
	}

	// Take the inverse transform
	fftw_execute(plan_reverse);

	summax = 0;
	imax = 0;
	// Find the peak pixel in the cross correlation
	for (int i = 0; i < npixsubap * npixsubap; i++) {
	    if (crosscorr[i] > summax) {
		imax = i;
		summax = crosscorr[i];
	    }
	}

	// Interpolate to find fractional pixel peak
	// y
	indx = ((iymax - 1) * nxkern + ixmax);
	if ( indx<0 ) indx += npixfftin;
	I1 = crosscorr[ indx % npixfftin ];
	indx += npixsubap;
	I2 = crosscorr[ indx % npixfftin ];
	indx += npixsubap;
	I3 = crosscorr[ indx % npixfftin ];
	ycens[iap] = iymax + (I1 - I3) / (I1 + I3 - 2 * I2) / 2;

	// x
	indx = (iymax * nxkern + ixmax - 1 );
	if ( indx<0 ) indx += npixfftin;
	I1 = crosscorr[ indx % npixfftin ];
	indx += 1;
	I2 = crosscorr[ indx % npixfftin ];
	indx += 1;
	I3 = crosscorr[ indx % npixfftin ];
	xcens[iap] = ixmax + (I1 - I3) / (I1 + I3 - 2 * I2) / 2;

	// This should give the centroid relative to the center of the image because after the reverse fft we didn't shift.


    }
    //  printf("Kernel calc: %f msec \n", (float)(sumclock) / CLOCKS_PER_SEC * 1000);

}

			
		    
		
