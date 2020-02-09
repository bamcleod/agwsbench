typedef float pixtype;

void shthresh(
	      int *subaps, // Array of subaperture locations subap[2*i]=y subap[2*1+1]=x
	      int nsubaps,  // Number of subapertures
	      int ny,          // Number of pixels in image y-axis
	      int nx,          // Number of pixels in image x-axis
	      int npixsubap,      // Number of pixels per subaperture
	      pixtype thresh,   // Threshold value
	      double *xcens,    // Return xcentoids
	      double *ycens     // Return xcentoids
	      );

void shcorrelate_pointerstodata(
		 int *subaps, // Array of subaperture locations subap[2*i]=y subap[2*1+1]=x
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
		 );

void shcorrelate_copydata(
		 int *subaps, // Array of subaperture locations subap[2*i]=y subap[2*1+1]=x
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
		 );

void shcorrelate_fft(
		 int *subaps, // Array of subaperture locations subap[2*i]=y subap[2*1+1]=x
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
		 );

