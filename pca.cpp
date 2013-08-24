#include "pca.h"
#include <iostream>
#include <iomanip>
#include "stdafx.h"
#include "test.h"

#include <Accelerate/Accelerate.h>
#include <vecLib/vBLAS.h>

using namespace Halide;

void SamplePCA();
Image<float> Covariance(Image<float> data);
void TestBLAS();
void TestLAPACK();
void dumpImage(std::string hdr, Image<float> img)
{
	std::cout << std::endl << "___" << hdr << "__" << std::endl;
	
	for (int i=0; i<img.extent(0); i++) {
		for(int j=0; j<img.extent(1); j++) {
			std::cout << img(i,j) << " ";
			
		}
		std::cout << std::endl;
	}
}
Image<float> transpose(Image<float> m) 
{
	Var x, y;
	Func f;
	f(x,y) = m(y,x);
	return f.realize( m.extent(1), m.extent(0) );
}

//
// Principal component analysis
// Based on http://www.cs.otago.ac.nz/cosc453/student_tutorials/principal_components.pdf
//

void TestPCA()
{
	// specify the dataset
	/*
	int N = 12; // samples
	int DIMS = 2;
	float rawdata[][12] = {
		{9, 15,  25, 14, 10, 18, 0,  16,  5, 19, 16, 20}, // "hours studied"
		{39, 56, 93, 61, 50, 75, 32, 85, 42, 70, 66, 80}  // "marks"
	};
	*/
	
	// Exercises (2.1)
	int N = 5;
	int DIMS = 2;
	float rawdata[][5] = {
		{10, 39, 19, 23, 28}, // x
		{43, 13, 32, 21, 20}  // y
	};
	
	// create an DIMxN Image<float> for the data set
	Image<float> data(DIMS,N);
	for(int i=0; i<N; i++) {
		for(int dim=0; dim<DIMS; dim++) {
			data(dim,i) = rawdata[dim][i];
		}
	}
	
	Covariance(data);

	// Test BLAS
	TestBLAS();
	TestLAPACK();
	
	// Run PCA from the PDF
	SamplePCA();
	
	
}

void TestBLAS()
{
	std::cout << "__TestBLAS__" << std::endl;
	
	float *X, *Y;
	int N = 10;
	X = new float[N];
	Y = new float[N];
	
	for(int i=0; i<N; i++) {
		X[i] = 1.0;
		Y[i] = 2.0;
	}
	
	float result = cblas_sdot( N, X, 1, Y, 1);
	std::cout << "BLAS dot product is: " << result << std::endl;
	
	delete[] X;
	delete[] Y;
}

void TestLAPACK()
{
	std::cout << "__TestLAPACK__" << std::endl;
	
	int VERS_MAJOR, VERS_MINOR, VERS_PATCH;
	ilaver_ (&VERS_MAJOR, &VERS_MINOR, &VERS_PATCH);
	
	std::cout << "LAPACK VERSION " << VERS_MAJOR << "." << VERS_MINOR << "." << VERS_PATCH << std::endl;
	
	// compute the eignvalues of a matrix... (COLUMN MAJOR)
	double cov[4] = { /* c0 */ 0.616555556, 0.615444444, /*c1*/ 0.615444444, 0.716555556 };
	
	// eigenvalues should be 0.0490833989 and 1.28402771
	
	// LAPACK makes me want to vomit
	char JOBZ = 'V'; // compute eigenvalues only; 'V' computes vectors too
	char UPLO = 'L'; // U: upper triangle of A is stored; 'L': lower...
	int N = 2; // order of matrix
	double* A = cov; // real arr 
	// NOTE: if JOBZ = 'N' the upper or lower triangle is destroyed
	int LDA = 2; // leading dimension of the array
	double W[2]; // *out* -- eignvalues
	double WORK[3*N-1]; // optimal LWORK ...
	int LWORK = 3*N-1;
	int INFO; // *out* -- info for failure
	
	// d="Double" sy="Symmetric" ev_=eigenvalues
	dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);
	
	std::cout << "Info from dsyev_ returned " << INFO << std::endl;
	
	if(INFO == 0) {
		std::cout << "Optimal LWORK = " << WORK[0] << std::endl;
		std::cout << "Eigenvalues are " << W[0] << ", " << W[1] << std::endl;
		std::cout << "Eigenvectors: " << std::endl;
		for(int i=0; i<2; i++) {
			for (int j=0; j<2; j++) {
				std::cout << A[N*j+i] << " ";
			}
			std::cout << std::endl;
		}
	}
}

// compute eigenvectors and values for a 2x2 covariance matrix
struct EigenVector2 {
public:
	float value;
	float x, y;
};
void ComputeEigenvectors2(Image<float> covarianceMatrix, EigenVector2 eigenVectors [2])
{
	if(covarianceMatrix.extent(1) != covarianceMatrix.extent(0))
		throw "covariance matrix must be square";
	if(covarianceMatrix.extent(0) != 2)
		throw "covariance matrix (for this routine) must only have order 2";
		
	int order = covarianceMatrix.extent(0);
		
	// compute the eignvalues of a matrix... (COLUMN MAJOR)
	double cov[4] = { 
		/* c0 */ covarianceMatrix(0,0), covarianceMatrix(1,0),
		/*c1*/ covarianceMatrix(1,0), covarianceMatrix(1,1)
	};
	
	// eigenvalues should be 0.0490833989 and 1.28402771
	
	// LAPACK makes me want to vomit
	char JOBZ = 'V'; // compute eigenvalues only; 'V' computes vectors too
	char UPLO = 'L'; // U: upper triangle of A is stored; 'L': lower... (doesn't matter here, it will be full)
	int N = order; // order of matrix
	double* A = cov; // real arr
	// NOTE: if JOBZ = 'N' the upper or lower triangle is destroyed
	int LDA = 2; // leading dimension of the array
	double W[2]; // *out* -- eignvalues
	double WORK[3*N-1]; // optimal LWORK ...
	int LWORK = 3*N-1;
	int INFO; // *out* -- info for failure
	
	// d="Double" sy="Symmetric" ev_=eigenvalues
	dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);
	
	std::cout << "Info from dsyev_ returned " << INFO << std::endl;
	
	// if successful, INFO contains a zero return code
	if(INFO == 0) {
		std::cout << "Optimal LWORK = " << WORK[0] << std::endl;
		std::cout << "Eigenvalues are " << W[0] << ", " << W[1] << std::endl;
		std::cout << "Eigenvectors: " << std::endl;
		for(int i=0; i<2; i++) {
			for (int j=0; j<2; j++) {
				std::cout << A[N*j+i] << " ";
			}
			std::cout << std::endl;
		}
		
		eigenVectors[0].value = (float)W[0];
		eigenVectors[1].value = (float)W[1];
		
		// COLUMN major, vectors are in the columns
		eigenVectors[0].x = (float)A[0];
		eigenVectors[0].y = (float)A[1];
		eigenVectors[1].x = (float)A[2];
		eigenVectors[1].y = (float)A[3]; // this MAY NOT BE RIGHT
	}	
	else {
		std::cout << "Error from LAPACK::dsyev_ " << INFO << std::endl;
		throw "Eigenvalue/vector computation failed (call to dsyev_ returned non-zero)";
	}
}

extern "C" float dump(float x)
{
	std::cout << x << std::endl;
	return x;
}
HalideExtern_1(float, dump, float);

// data is DIMxSAMPLES
Image<float> Covariance(Image<float> data)
{
	int DIMS = data.extent(0);
	int N = data.extent(1);
	
	std::cout << "DIMS=" << DIMS << " SAMPLES(N)=" << N << std::endl;
	
	// compute the mean
	Func f_mean; 
	Var dim; 
	RDom r(0,N);
	f_mean(dim) = sum( data(dim, r.x) )/(float)N;
	Image<float> mean = f_mean.realize(DIMS); // mean of both dimensions
	
	std::cout << "mean(H,M) = " << mean(0) << ", " << mean(1) << std::endl;
	
	// compute the covariance of the two samples
	/*
	Func f_covar;
	f_covar() = sum( (data(_h,r.x)-mean(_h)) * (data(_m,r.x)-mean(_m)));
	Image<float> covar = f_covar.realize();
	
	std::cout << "covar(H,M) = " << covar(0) << std::endl;
	*/
	
	// compute a covariance matrix
	Func f_covarmx;
	Var x, y; // x = first dimension, y = second dimension... 
	// NOTE calculation is suplicated since (x,y) is the same as (y,x)
	f_covarmx(x,y) = sum( (data(x,r.x)-mean(x)) * (data(y,r.x)-mean(y)) ) / Halide::cast<float>(N-1);
	
	Image<float> covarmx = f_covarmx.realize(DIMS,DIMS);
	
	std::cout << "Covariance Matrix" << std::endl;
	for(int i=0; i<DIMS; i++) {
		for(int j=0; j<DIMS; j++) {
			std::cout << covarmx(i,j) << " ";
		}
		std::cout << std::endl;
	}
	
	return covarmx;
}

void SamplePCA()
{
	float data_x[] = {2.5f, 0.5f, 2.2f, 1.9f, 3.1f, 2.3f, 2,    1,    1.5f, 1.1f};
	float data_y[] = {2.4f, 0.7f, 2.9f, 2.2f, 3,    2.7f, 1.6f, 1.1f, 1.6f, 0.9f};
	int N = 10;
	
	// samples in COLUMNS (ROW MAJOR)
	// 2xN
	Image<float> data(2,N);
	for(int i=0; i<N; i++) {
		data(0,i) = data_x[i];
		data(1,i) = data_y[i];
	}
	
	RDom r(0,N);
	Var x, y, dim;

	// sample dot product (for kicks)
	/*
	Func dot;
	dot() = sum( data(0,r.x) * data(1,r.x) );
	Image<float> result = dot.realize();
	std::cout << result(0) << std::endl;
	*/
	
	Func mean;
	mean(dim) = sum( data(dim,r.x) ) / (float)N;

	Image<float> data_mean = mean.realize(2);
	std::cout << "Mean(x,y) = " << data_mean(0) << ", " << data_mean(1) << std::endl;

	// data with mean subtracted
	Func subtract_mean;
	subtract_mean(x,y) = data(x,y) - mean(x);
	Image<float> data_less_mean = subtract_mean.realize( 2, data.extent(1));
	
	dumpImage( "data_less_mean", data_less_mean);
	
	std::cout << "___ Covariance of Data ___" << std::endl;
	
	// ALL THE ABOVE is done (duplicated) in the Covariance calculation
	Image<float> covarmx = Covariance( data );
	
	// Eigenvalue/Eigenvector calculation is beyond the scope of implementation
	// so I'm using BLAS/LAPACK as shipped with Mac OSX Accellerate framework

	EigenVector2 ev[2];
	ComputeEigenvectors2(covarmx, ev);
	
	for(int i=0; i<2; i++) {
		std::cout << "EigenValue[" << i << "] " << ev[i].value << " : (" << ev[i].x << "," << ev[i].y << ")" << std::endl;
	}
	
	// reproduce the data with only the principal component
	// easy since we only have 2 components
	int principal = (ev[0].value > ev[1].value) ? 0 : 1;
	int other = 1-principal;
	
	// feature vector, eigenvectors in columns
	//  2xF
	Image<float> fv(2,2); // ROW MAJOR
	fv(0,0) = ev[principal].x;	// eig1
	fv(1,0) = ev[principal].y;
	
	// leave this out for demonstrating reduced dimensionality
	fv(0,1) = ev[other].x;		// eig2
	fv(1,1) = ev[other].y;
	
	dumpImage( "feature vector", fv );
	
	// Fx2 (F = 1, number of chosen features)
	Func f_transpose_fv;
	f_transpose_fv(x,y) = fv(y,x); // transpose the feature vector to get it in Fx2 (F = number of chosen features)
	
	std::cout << "_TRANSPOSED FEATURE VECTOR" << std::endl;
	dumpImage( "feature_vector(transposed)", f_transpose_fv.realize(  fv.extent(1), fv.extent(0)));
	
	dumpImage( "DataAdjust:", data_less_mean);
	
	// 2xN -- NOTE data_less_mean is already 2xN
	//Func f_transpose_dataAdjust;
	//f_transpose_dataAdjust(x,y) = data_less_mean(y,x);
	
	//dumpImage( "DataAdjust(transposed)", f_transpose_dataAdjust.realize( data_less_mean.extent(1), data_less_mean.extent(0)));
	
	Func f_mul; 
	RDom ri(0,fv.extent(0)); 				// prior to transposition, this is (0,2)
	
	Var i, j;
	// multiple reduction domains found
	f_mul(i,j) = sum( f_transpose_fv(i, ri.x) * data_less_mean(ri.x, j) ); // matrix multiplication
	
	
	Image<float> finalData = f_mul.realize( fv.extent(1), data_less_mean.extent(1) );  // (1,N)
	
	dumpImage( "Final Data (reduced dimensions)", finalData);
	
	dumpImage( "Final Data (in table form for spreadsheet)", transpose(finalData));
}












































