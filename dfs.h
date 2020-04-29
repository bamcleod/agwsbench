#define NROWS 256
#define NCOLS 320

typedef float pixtype;

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
		     );
