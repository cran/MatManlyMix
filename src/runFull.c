
#include "array.h"
#include "MatManlyMix.h"


void run_Mstep_Manly_Full(double *y, double *misc_double, int *misc_int, double *gamma1, double *la1, double *nu1, double *invS1, double *Mu1, double *invPsi1, double *detS, double *detPsi, double *tau){
		
	double ***Y, **la, **nu, ***Mu, ***invS, ***invPsi, **gamma;

	int p, n, T, K, Mu_type;


	p = misc_int[0];
	T = misc_int[1];
	n = misc_int[2];
	K = misc_int[3];	
	Mu_type = misc_int[4];	

	MAKE_3ARRAY(Y, p, T, n);
	MAKE_3ARRAY(Mu, p, T, K);
	MAKE_3ARRAY(invS, p, p, K);
	MAKE_3ARRAY(invPsi, T, T, K);
	MAKE_MATRIX(gamma, n, K);
	MAKE_MATRIX(la, K, p);
	MAKE_MATRIX(nu, K, T);
	
	array1to3(p, T, n, y, Y);
	array1to3(p, T, K, Mu1, Mu);
	array1to3(p, p, K, invS1, invS);
	array1to3(T, T, K, invPsi1, invPsi);
	array1to2(n, K, gamma1, gamma);
	array1to2(K, p, la1, la);
	array1to2(K, T, nu1, nu);


	Mstep_Manly_Full(p, T, n, K, misc_double, Y, la, nu, gamma, invS, Mu, invPsi, detS, detPsi, tau, Mu_type);


	//printf("run M step la %lf %lf\n", la[0][0], la[0][1]);	


	array3to1(p, T, n, y, Y);
	array3to1(p, T, K, Mu1, Mu);
	array3to1(p, p, K, invS1, invS);
	array3to1(T, T, K, invPsi1, invPsi);
	array2to1(n, K, gamma1, gamma);
	array2to1(K, p, la1, la);
	array2to1(K, T, nu1, nu);

		
	FREE_3ARRAY(Y);
	FREE_3ARRAY(Mu);
	FREE_3ARRAY(invS);
	FREE_3ARRAY(invPsi);
	FREE_MATRIX(gamma);
	FREE_MATRIX(la);
	FREE_MATRIX(nu);		
}




void run_Mstep_Manly_Reg(double *y, double *misc_double, int *misc_int, double *gamma1, double *la1, double *nu1, double *invS1, double *x, double *beta1, double *invPsi1, double *detS, double *detPsi, double *tau){
		
	double ***Y, ***X, ***beta, **la, **nu, ***invS, ***invPsi, **gamma;

	int p, n, T, q, K;


	p = misc_int[0];
	T = misc_int[1];
	n = misc_int[2];
	q = misc_int[3];
	K = misc_int[4];	

	MAKE_3ARRAY(Y, p, T, n);
	MAKE_3ARRAY(X, T, q, n);
	MAKE_3ARRAY(beta, q, p, K);
	MAKE_3ARRAY(invS, p, p, K);
	MAKE_3ARRAY(invPsi, T, T, K);
	MAKE_MATRIX(gamma, n, K);
	MAKE_MATRIX(la, K, p);
	MAKE_MATRIX(nu, K, T);
	
	array1to3(p, T, n, y, Y);
	array1to3(T, q, n, x, X);
	array1to3(q, p, K, beta1, beta);
	array1to3(p, p, K, invS1, invS);
	array1to3(T, T, K, invPsi1, invPsi);
	array1to2(n, K, gamma1, gamma);
	array1to2(K, p, la1, la);
	array1to2(K, T, nu1, nu);

	Mstep_Manly_Reg(p, T, n, q, K, misc_double, Y, la, nu, gamma, invS, X, beta, invPsi, detS, detPsi, tau);

	//printf("run M step la %lf %lf\n", la[0][0], la[0][1]);	


	array3to1(p, T, n, y, Y);
	array3to1(T, q, n, x, X);
	array3to1(q, p, K, beta1, beta);
	array3to1(p, p, K, invS1, invS);
	array3to1(T, T, K, invPsi1, invPsi);
	array2to1(n, K, gamma1, gamma);
	array2to1(K, p, la1, la);
	array2to1(K, T, nu1, nu);
		
	FREE_3ARRAY(Y);
	FREE_3ARRAY(beta);
	FREE_3ARRAY(invS);
	FREE_3ARRAY(invPsi);
	FREE_MATRIX(gamma);
	FREE_MATRIX(la);
	FREE_MATRIX(nu);
	FREE_3ARRAY(X);		
}




void run_Mat_Manly_Full(double *y, int *misc_int, double *misc_double, double *tau, double *la1, double *nu1, double *Mu1, double *invS1, double *invPsi1, double *detS, double *detPsi, double *gamma1, int *id, double *ll, int *conv){
		
	double ***Y, **la, **nu, ***Mu, ***invS, ***invPsi, **gamma;

	int p, n, T, K, max_iter, Mu_type;


	p = misc_int[0];
	T = misc_int[1];
	n = misc_int[2];
	K = misc_int[3];	
	max_iter = misc_int[4];	
	Mu_type = misc_int[5];	

	MAKE_3ARRAY(Y, p, T, n);
	MAKE_3ARRAY(Mu, p, T, K);
	MAKE_3ARRAY(invS, p, p, K);
	MAKE_3ARRAY(invPsi, T, T, K);
	MAKE_MATRIX(gamma, n, K);
	MAKE_MATRIX(la, K, p);
	MAKE_MATRIX(nu, K, T);
	
	array1to3(p, T, n, y, Y);
	array1to3(p, T, K, Mu1, Mu);
	array1to3(p, p, K, invS1, invS);
	array1to3(T, T, K, invPsi1, invPsi);
	array1to2(n, K, gamma1, gamma);
	array1to2(K, p, la1, la);
	array1to2(K, T, nu1, nu);


	EM_Manly_Full(p, T, n, K, Y, la, nu, max_iter, misc_double, tau, Mu, invS, invPsi, detS, detPsi, gamma, id, ll, conv, Mu_type);


	array3to1(p, T, n, y, Y);
	array3to1(p, T, K, Mu1, Mu);
	array3to1(p, p, K, invS1, invS);
	array3to1(T, T, K, invPsi1, invPsi);
	array2to1(n, K, gamma1, gamma);
	array2to1(K, p, la1, la);
	array2to1(K, T, nu1, nu);

		
	FREE_3ARRAY(Y);
	FREE_3ARRAY(Mu);
	FREE_3ARRAY(invS);
	FREE_3ARRAY(invPsi);
	FREE_MATRIX(gamma);
	FREE_MATRIX(la);
	FREE_MATRIX(nu);
		
}



void run_Mat_Manly_Reg(double *y, int *misc_int, double *misc_double, double *tau, double *la1, double *nu1, double *x, double *beta1, double *invS1, double *invPsi1, double *detS, double *detPsi, double *gamma1, int *id, double *ll, int *conv){
		
	double ***Y, ***X, ***beta, **la, **nu, ***invS, ***invPsi, **gamma;

	int p, n, T, q, K, max_iter;


	p = misc_int[0];
	T = misc_int[1];
	n = misc_int[2];
	q = misc_int[3];
	K = misc_int[4];	
	max_iter = misc_int[5];	
	

	MAKE_3ARRAY(Y, p, T, n);
	MAKE_3ARRAY(X, T, q, n);
	MAKE_3ARRAY(beta, q, p, K);
	MAKE_3ARRAY(invS, p, p, K);
	MAKE_3ARRAY(invPsi, T, T, K);
	MAKE_MATRIX(gamma, n, K);
	MAKE_MATRIX(la, K, p);
	MAKE_MATRIX(nu, K, T);
	
	array1to3(p, T, n, y, Y);
	array1to3(T, q, n, x, X);
	array1to3(q, p, K, beta1, beta);
	array1to3(p, p, K, invS1, invS);
	array1to3(T, T, K, invPsi1, invPsi);
	array1to2(n, K, gamma1, gamma);
	array1to2(K, p, la1, la);
	array1to2(K, T, nu1, nu);


	EM_Manly_Reg(p, T, n, q, K, Y, la, nu, max_iter, misc_double, tau, X, beta, invS, invPsi, detS, detPsi, gamma, id, ll, conv);



	array3to1(p, T, n, y, Y);
	array3to1(T, q, n, x, X);
	array3to1(q, p, K, beta1, beta);
	array3to1(p, p, K, invS1, invS);
	array3to1(T, T, K, invPsi1, invPsi);
	array2to1(n, K, gamma1, gamma);
	array2to1(K, p, la1, la);
	array2to1(K, T, nu1, nu);

		
	FREE_3ARRAY(Y);
	FREE_3ARRAY(beta);
	FREE_3ARRAY(invS);
	FREE_3ARRAY(invPsi);
	FREE_MATRIX(gamma);
	FREE_MATRIX(la);
	FREE_MATRIX(nu);
	FREE_3ARRAY(X);		
}

