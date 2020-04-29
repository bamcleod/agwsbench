// g++ -I/home/bmcleod/include -O lstsq_eigen.cpp -o lstsq_eigen

// or, even better:

// g++ -I /home/bmcleod/include -O2 -mavx512f -mavx512cd -fopt-info-vec-all lstsq_eigen.cpp -o lstsq_eigen

// similar time:

// g++ -I /home/bmcleod/include -O3 -mavx512f -mavx512cd -fopt-info-vec-all -march=native lstsq_eigen.cpp -o lstsq_eigen


#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

#define NFRAMES 3000
#define NROWS    256
#define NCOLS    320
#define NFOWLER    7
#define NSKIP      1

void preprocess(unsigned short int *data, int nframes, int nrows, int ncols, float *slopes)
{
    int i;
    int f;
    int r;
    int c;
    
    MatrixXf  A(NFOWLER-NSKIP,2);
    VectorXf  b(NFOWLER-NSKIP);

    for ( i=0; i<NFOWLER-NSKIP; i++) {
        A(i,0) = 1; //intercept
        A(i,1) = i; // slope
    }
    MatrixXf At    = A.transpose();
    MatrixXf AtA   = At * A;
    MatrixXf Ainv  = AtA.inverse() * At;

    // Now fit the full matrix/
    for ( f=0; f<nframes; f++) {
	for ( r=0; r<nrows; r++) { 
	    for ( c=0; c<ncols; c++) {
		for ( i=0; i<NFOWLER-NSKIP; i++) {
		    b(i) = data[(f*NFOWLER+i) * nrows * ncols + r * ncols + c];
		}
		slopes[f * nrows * ncols + r * ncols + c] = (Ainv * b)[1];
	    }
	}
    }
}

int main()
{
    MatrixXf  A(NFOWLER-NSKIP,2);
    VectorXf  b(NFOWLER-NSKIP);
    int       i, f, r ,c;
    unsigned short int *data;
    float     *slopes;

    data   = (unsigned short int *)malloc( NFRAMES * NROWS * NCOLS * NFOWLER * sizeof(short int));
    slopes =     (float *)malloc( NFRAMES * NROWS * NCOLS *           sizeof(float));

    /*    

    for ( i=0; i<NFOWLER-NSKIP; i++) {
        A(i,0) = 1; //intercept
        A(i,1) = i; // slope
    }
    MatrixXf At    = A.transpose();
    MatrixXf AtA   = At * A;
    MatrixXf Ainv  = AtA.inverse() * At;

    // Make some example fake data
    for ( i=0; i<NFOWLER-NSKIP; i++) {
	b(i) += i * 100 + 30 ; 
    }
    b += VectorXf::Random(NFOWLER-NSKIP);
    
    cout << "Here is the matrix A:\n" << A << endl;
    cout << "Here is the right hand side b:\n" << b << endl;
    cout << "The solution is:\n" <<  Ainv * b;
    cout << "\n";
    */
    
    preprocess(data, NFRAMES, NROWS, NCOLS, slopes);
	
}

