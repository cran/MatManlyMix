

#include "array.h"
#include "math.h"
#include "MatManlyMix.h"
 
#define Inf 1e+140


#ifdef __HAVE_R_
	#include <R.h>
	#include <Rmath.h>
#endif


// Manly transformation for matrix Y
void Manly_trans(int p, int T, double *la, double *nu, double **Y, double **MY){

  int j, t;

  for (j=0; j<p; j++){
    for(t=0; t<T; t++) {
      if (fabs(la[j] + nu[t])<1e-12){
	MY[j][t] = Y[j][t];
      }
      else{
     	MY[j][t] = (exp(Y[j][t] * (la[j] + nu[t])) - 1) / (la[j] + nu[t]);
      }
    }
      
  }

}

void Manly_trans_whole(int n, int p, int T, double *la, double *nu, double ***Y, double ***MY){

  int i, j, t;
 for(i=0; i<n; i++){


  for (j=0; j<p; j++){
    for(t=0; t<T; t++) {
      if (fabs(la[j] + nu[t])<1e-12){
	MY[j][t][i] = Y[j][t][i];
      }
      else{
	MY[j][t][i] = (exp(Y[j][t][i] * (la[j] + nu[t])) - 1) / (la[j] + nu[t]);
      }
    }
      
  }
 }
  
}


void cpyk(double ***a, int nrows, int ncols, int k, double **b)
{
  int i,j;
  for(i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      b[i][j]=a[i][j][k];
    }
  }

}


/*copies matrix A to matrix B*/

void cpy(double **a, int nrows, int ncols, double **b)
{
  int i,j;
  for(i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      b[i][j]=a[i][j];
    }
  }

}




/* Multiplies matrix a and vector x and puts the result in y which should be
 pre-allocated */

void matxvec(double **a, int arows, int acols,
		double *x, int xrows, double *y)
{
  int i, k;
  
  for (i=0; i<arows; i++){
    y[i] = 0;
    for (k=0; k<acols; k++){
      y[i] += a[i][k] * x[k];
    }
  }

}



void cpyk2(double **a, int nrows, int ncols, double ***b, int k)
{
  int i,j;
  for(i=0;i<nrows;i++) {
    for (j=0;j<ncols;j++) {
      b[i][j][k]=a[i][j];
    }
  }

}


double mGpdf_Manly_Full(int p, int T, double *la, double *nu, double **Y, double **Mu, double **invS, double **invPsi, double detS, double detPsi){

	int j, t;
	double trace = 0.0, dens, jac = 0.0;
	double **MY, **tMY, **temp1, **temp2, **maha;

	MAKE_MATRIX(maha, p, p);
	MAKE_MATRIX(temp1, p, T);
	MAKE_MATRIX(temp2, p, T);
	MAKE_MATRIX(tMY, T, p);	
	MAKE_MATRIX(MY, p, T);


	
  	Manly_trans(p, T, la, nu, Y, MY);
	mat_(p, T, MY, Mu);

	tA(MY, T, p, tMY);

	multiply(invS, p, p, MY, p, T, temp1);
	
	multiply(temp1, p, T, invPsi, T, T, temp2);	

	multiply(temp2, p, T, tMY, T, p, maha);	


	for(j=0; j<p; j++){	 

		trace += maha[j][j];
	}
	
	dens = 1.0 / pow((2*PI), p*T/2.0) / pow(detS, T/2.0) / pow(detPsi, p/2.0) * exp(-1.0 / 2.0 * trace);


	for(j=0; j<p; j++){

		for(t=0; t<T; t++){

			jac = jac + (la[j] + nu[t]) * Y[j][t];
	
		}	
	}

	dens = dens * exp(jac);


	FREE_MATRIX(maha);
	FREE_MATRIX(temp1);
	FREE_MATRIX(temp2);
	FREE_MATRIX(tMY);
	FREE_MATRIX(MY);
	return dens;
}




double mGloglik_Manly_Full(int p, int T, int n, int K, double ***Y, double **la, double **nu, double *tau, double ***Mu, double ***invS, double ***invPsi, double *detS, double *detPsi){

	int i,k;
	double loglik = 0.0;
	double ll, dens = 0.0;
	double **Yi, **Muk, **invSk, **invPsik;
	
	MAKE_MATRIX(Yi, p, T);
	MAKE_MATRIX(Muk, p, T);
	MAKE_MATRIX(invSk, p, p);
	MAKE_MATRIX(invPsik, T, T);


	for(i=0; i<n; i++){	 

		ll = 0;

		for(k=0; k<K; k++){

			cpyk(Y, p, T, i, Yi);
			cpyk(Mu, p, T, k, Muk);
			cpyk(invS, p, p, k, invSk);
			cpyk(invPsi, T, T, k, invPsik);

			dens = mGpdf_Manly_Full(p, T, la[k], nu[k], Yi, Muk, invSk, invPsik, detS[k], detPsi[k]);


			ll += tau[k] * dens;


		}

		loglik += log(ll);
	}

	FREE_MATRIX(Yi);
	FREE_MATRIX(Muk);
	FREE_MATRIX(invSk);
	FREE_MATRIX(invPsik);

	return loglik;
}






void Estep_Manly_Full(int p, int T, int n, int K, double ***Y, double **la, double **nu, double *tau, double ***Mu, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma){

	int i,k;		
	double dens =0.0, *sum_gamma, **Yi, **Muk, **invSk, **invPsik;

	MAKE_MATRIX(Yi, p, T);
	MAKE_MATRIX(Muk, p, T);
	MAKE_MATRIX(invSk, p, p);
	MAKE_MATRIX(invPsik, T, T);

	MAKE_VECTOR(sum_gamma, n);

	anull(sum_gamma, n);

	if(K == 1){

		for(i=0; i<n; i++){	 
			gamma[i][0] = 1.0;

		}
	}
	else{

		for(i=0; i<n; i++){	 



			for(k=0; k<K; k++){
				cpyk(Y, p, T, i, Yi);
				cpyk(Mu, p, T, k, Muk);
				cpyk(invS, p, p, k, invSk);
				cpyk(invPsi, T, T, k, invPsik);

				dens = mGpdf_Manly_Full(p, T, la[k], nu[k], Yi, Muk, invSk, invPsik, detS[k], detPsi[k]);

				gamma[i][k] = tau[k] * dens;

				sum_gamma[i] += gamma[i][k];


			}

	
			for(k=0; k<K; k++){

				gamma[i][k] = gamma[i][k] / sum_gamma[i];


			}


		}

	}
	FREE_MATRIX(Yi);
	FREE_MATRIX(Muk);
	FREE_MATRIX(invSk);
	FREE_MATRIX(invPsik);

	FREE_VECTOR(sum_gamma);

}








// Q-function
double Q1(int n, int p, int T, double *la_nonzero, double *nu, int *index, double ***Y, double *gamma_k, double **invPsik, int Mu_type){
  
	int i, j, t, j1, j2, count;
	double res, jac, sum_gamma, det1, det2, det;
	double *la;
	double **Sk, **temp1, **temp2;
	double ***MY, **Muk, *Eig1, *Eig2, **MYi, **tMYi, **invPsik0, *Muk_vec, *MY_vec, *Eig, **L, **matconst, **invconst;


	MAKE_VECTOR(Muk_vec, p+T-1);
	MAKE_VECTOR(MY_vec, p+T-1);
	MAKE_MATRIX(matconst, p+T-1, p+T-1);
	MAKE_MATRIX(invconst, p+T-1, p+T-1);
	MAKE_MATRIX(L, p+T-1, p+T-1);
	MAKE_VECTOR(Eig, p+T-1);


		
	MAKE_VECTOR(la, p);
	MAKE_VECTOR(Eig1, T);
	MAKE_VECTOR(Eig2, p);
	MAKE_3ARRAY(MY, p,T,n);
	MAKE_MATRIX(MYi, p, T);
	MAKE_MATRIX(tMYi, T, p);
	MAKE_MATRIX(Sk, p, p);
	MAKE_MATRIX(invPsik0, T, T);
	MAKE_MATRIX(temp1, p, T);
	MAKE_MATRIX(temp2, p, p);

	MAKE_MATRIX(Muk, p, T);

	Anull(matconst, p+T-1, p+T-1);
	for(i=0; i<p; i++){	
		for(j=0; j<p; j++){

			if(i == j){
				matconst[i][j] = T*1.0;
			}

		}
	}

	for(i=p; i<p+T-1; i++){	
		for(j=p; j<p+T-1; j++){

			if(i == j){
				matconst[i][j] = p*1.0;
			}

		}
	}

	for(i=0; i<p; i++){	
		for(j=p; j<p+T-1; j++){
			matconst[i][j] = 1.0;

		}
	}
	for(i=p; i<p+T-1; i++){	
		for(j=0; j<p; j++){
			matconst[i][j] = 1.0;

		}
	}
	anull(Eig, p+T-1);

	#ifdef __HAVE_R_
		EigValDec(p+T-1, Eig, matconst, &det);
	#else
		cephes_symmeigens_down(p+T-1, Eig, matconst, &det);
	#endif

	Anull(L, p+T-1, p+T-1);

	for (t=0; t<p+T-1; t++){
		L[t][t] = 1.0 / Eig[t];
		
	}
	
	XAXt(matconst, p+T-1, L, invconst);



	count = 0;



	cpy(invPsik, T, T, invPsik0);

	for(j=0; j<p; j++){
		if(index[j]==1){
			la[j] = la_nonzero[count];
			count += 1;
		}
		else{
			la[j] = 0; 

		}
	}	


	sum_gamma = 0;
			
	for(i=0; i<n; i++){


		sum_gamma += gamma_k[i];
	}


	Manly_trans_whole(n, p, T, la, nu, Y, MY);
	res = 0;
	Anull(Muk, p, T);
	if(Mu_type == 0){

		for(i=0; i<n; i++){
			for(j=0; j<p; j++){

				for(t=0; t<T; t++){
			
					Muk[j][t] += gamma_k[i] * MY[j][t][i] / sum_gamma;


				}
			}
		}


	}
	else if(Mu_type == 1){

		anull(MY_vec, p+T-1);
	

		for(j=0; j<p; j++){
			for(i=0; i<n; i++){
				for(t=0; t<T; t++){		
					MY_vec[j] += gamma_k[i] * MY[j][t][i] / sum_gamma;
	//printf(" Q  %lf\n",MY_vec[j]);


				}
			}
		}


		for(t=p; t<p+T-1; t++){
			for(i=0; i<n; i++){
				for(j=0; j<p; j++){		
					MY_vec[t] += gamma_k[i] * MY[j][t-p][i] / sum_gamma;
//printf(" Q  %lf\n",MY_vec[t]);

				}
			}
		}


		matxvec(invconst, p+T-1, p+T-1, MY_vec, p+T-1, Muk_vec);

		//printf(" invconst %lf  %lf  %lf  %lf\n",invconst[0][0],invconst[1][1],invconst[0][2],invconst[2][0]);

		for(j=0; j<p; j++){

			for(t=0; t<T; t++){
			
				if(t<T-1){
					Muk[j][t] = Muk_vec[j] +  Muk_vec[p+t];

				}
				else{
					Muk[j][t] = Muk_vec[j];
				}
			}
		}


	}



	
	Anull(Sk, p, p);
	
	for(i=0; i<n; i++){


		cpyk(MY, p, T, i, MYi);
		mat_(p, T, MYi, Muk);
	
		tA(MYi, T, p, tMYi);


		multiply(MYi, p, T, invPsik0, T, T, temp1);

		multiply(temp1, p, T, tMYi, T, p, temp2);


	
		for(j1=0; j1<p; j1++){			
	
			for(j2=0; j2<p; j2++){
	
				Sk[j1][j2] += gamma_k[i] * temp2[j1][j2] / T /sum_gamma;


			}
		}

	}
	



	anull(Eig1, T);

	#ifdef __HAVE_R_
		EigValDec(T, Eig1, invPsik0, &det1);
	#else
		cephes_symmeigens_down(T, Eig1, invPsik0, &det1);
	#endif



	anull(Eig2, p);

	#ifdef __HAVE_R_
		EigValDec(p, Eig2, Sk, &det2);
	#else
		cephes_symmeigens_down(p, Eig2, Sk, &det2);
	#endif



	res = -p / 2.0 * sum_gamma * log(1.0 / det1) -T / 2.0 * sum_gamma * log(det2); 

	
	for(i=0; i<n; i++){

		jac = 0;

		for(j=0; j<p; j++){

			for(t=0; t<T; t++){

				jac += Y[j][t][i] * (la[j] + nu[t]); 
			}

		}


		res = res + gamma_k[i] * jac;
	}



	FREE_VECTOR(la);
	FREE_MATRIX(MYi);
	FREE_3ARRAY(MY);
	FREE_MATRIX(tMYi);
	FREE_VECTOR(Eig1);
	FREE_VECTOR(Eig2);

	FREE_MATRIX(Sk);
	FREE_MATRIX(invPsik0);
	FREE_MATRIX(temp1);
	FREE_MATRIX(temp2);
	FREE_MATRIX(Muk);
	FREE_MATRIX(matconst);
	FREE_MATRIX(invconst);
	FREE_MATRIX(L);
	FREE_VECTOR(Eig);
	FREE_VECTOR(Muk_vec);
	FREE_VECTOR(MY_vec);
	return(-res);

}







// Q-function
double Q2(int n, int p, int T, double *nu_nonzero, double *la, int *index, double ***Y, double *gamma_k, double **invPsik, int Mu_type){
  
	int i, j, t, j1, j2, count;
	double res, jac, sum_gamma, det1, det2, det;
	double *nu;
	double **Sk, **temp1, **temp2;
	double ***MY, **Muk, *Eig1, *Eig2, **MYi, **tMYi, **invPsik0, *Muk_vec, *MY_vec, *Eig, **L, **matconst, **invconst;

	MAKE_VECTOR(Muk_vec, p+T-1);
	MAKE_VECTOR(MY_vec, p+T-1);
	MAKE_MATRIX(matconst, p+T-1, p+T-1);
	MAKE_MATRIX(invconst, p+T-1, p+T-1);
	MAKE_MATRIX(L, p+T-1, p+T-1);
	MAKE_VECTOR(Eig, p+T-1);


		
	MAKE_VECTOR(nu, T);
	MAKE_VECTOR(Eig1, T);
	MAKE_VECTOR(Eig2, p);
	MAKE_3ARRAY(MY, p,T,n);
	MAKE_MATRIX(MYi, p, T);
	MAKE_MATRIX(tMYi, T, p);
	MAKE_MATRIX(Sk, p, p);
	MAKE_MATRIX(invPsik0, T, T);
	MAKE_MATRIX(temp1, p, T);
	MAKE_MATRIX(temp2, p, p);

	MAKE_MATRIX(Muk, p, T);

	//printf(" Q2 %lf %lf\n", nu_nonzero[0], nu_nonzero[1]);
	//printf(" Q2 invPsik %lf  %lf  %lf  %lf\n",invPsik[0][0],invPsik[1][1],invPsik[0][2],invPsik[2][0]);



	Anull(matconst, p+T-1, p+T-1);
	for(i=0; i<p; i++){	
		for(j=0; j<p; j++){

			if(i == j){
				matconst[i][j] = T*1.0;
			}

		}
	}

	for(i=p; i<p+T-1; i++){	
		for(j=p; j<p+T-1; j++){

			if(i == j){
				matconst[i][j] = p*1.0;
			}

		}
	}

	for(i=0; i<p; i++){	
		for(j=p; j<p+T-1; j++){
			matconst[i][j] = 1.0;

		}
	}
	for(i=p; i<p+T-1; i++){	
		for(j=0; j<p; j++){
			matconst[i][j] = 1.0;

		}
	}
	anull(Eig, p+T-1);

	#ifdef __HAVE_R_
		EigValDec(p+T-1, Eig, matconst, &det);
	#else
		cephes_symmeigens_down(p+T-1, Eig, matconst, &det);
	#endif

	Anull(L, p+T-1, p+T-1);

	for (t=0; t<p+T-1; t++){
		L[t][t] = 1.0 / Eig[t];
		
	}
	
	XAXt(matconst, p+T-1, L, invconst);



	count = 0;


	cpy(invPsik, T, T, invPsik0);

	for(t=0; t<T; t++){
		if(index[t]==1){
			nu[t] = nu_nonzero[count];
			count += 1;
		}
		else{
			nu[t] = 0; 

		}
	}	


	sum_gamma = 0;
			
	for(i=0; i<n; i++){


		sum_gamma += gamma_k[i];
	}


	Manly_trans_whole(n, p, T, la, nu, Y, MY);
	res = 0;
	Anull(Muk, p, T);

	if(Mu_type == 0){

		for(i=0; i<n; i++){
			for(j=0; j<p; j++){

				for(t=0; t<T; t++){
			
					Muk[j][t] += gamma_k[i] * MY[j][t][i] / sum_gamma;


				}
			}
		}


	}
	else if(Mu_type == 1){

		anull(MY_vec, p+T-1);
	

		for(j=0; j<p; j++){
			for(i=0; i<n; i++){
				for(t=0; t<T; t++){		
					MY_vec[j] += gamma_k[i] * MY[j][t][i] / sum_gamma;
	//printf(" Q  %lf\n",MY_vec[j]);


				}
			}
		}


		for(t=p; t<p+T-1; t++){
			for(i=0; i<n; i++){
				for(j=0; j<p; j++){		
					MY_vec[t] += gamma_k[i] * MY[j][t-p][i] / sum_gamma;
//printf(" Q  %lf\n",MY_vec[t]);

				}
			}
		}


		matxvec(invconst, p+T-1, p+T-1, MY_vec, p+T-1, Muk_vec);

		//printf(" invconst %lf  %lf  %lf  %lf\n",invconst[0][0],invconst[1][1],invconst[0][2],invconst[2][0]);

		for(j=0; j<p; j++){

			for(t=0; t<T; t++){
			
				if(t<T-1){
					Muk[j][t] = Muk_vec[j] +  Muk_vec[p+t];

				}
				else{
					Muk[j][t] = Muk_vec[j];
				}
			}
		}


	}




	
	Anull(Sk, p, p);
	
	for(i=0; i<n; i++){


		cpyk(MY, p, T, i, MYi);
		mat_(p, T, MYi, Muk);
	
		tA(MYi, T, p, tMYi);

		multiply(MYi, p, T, invPsik0, T, T, temp1);


		multiply(temp1, p, T, tMYi, T, p, temp2);

	
		for(j1=0; j1<p; j1++){			
	
			for(j2=0; j2<p; j2++){
	
				Sk[j1][j2] += gamma_k[i] * temp2[j1][j2] / T /sum_gamma;


			}
		}

	}
	



	anull(Eig1, T);

	#ifdef __HAVE_R_
		EigValDec(T, Eig1, invPsik0, &det1);
	#else
		cephes_symmeigens_down(T, Eig1, invPsik0, &det1);
	#endif

	//printf(" Q Psi %lf\n", det1);


	anull(Eig2, p);

	#ifdef __HAVE_R_
		EigValDec(p, Eig2, Sk, &det2);
	#else
		cephes_symmeigens_down(p, Eig2, Sk, &det2);
	#endif

	//printf(" Q S %lf\n", det2);
	//printf(" Q sum_gamma %lf\n", sum_gamma);


	res = -p / 2.0 * sum_gamma * log(1.0 / det1) -T / 2.0 * sum_gamma * log(det2); 

	
	for(i=0; i<n; i++){

		jac = 0;

		for(j=0; j<p; j++){

			for(t=0; t<T; t++){

				jac += Y[j][t][i] * (la[j] + nu[t]); 
			}

		}


		res = res + gamma_k[i] * jac;
	}



	FREE_VECTOR(nu);
	FREE_MATRIX(MYi);
	FREE_3ARRAY(MY);
	FREE_MATRIX(tMYi);
	FREE_VECTOR(Eig1);
	FREE_VECTOR(Eig2);
	FREE_MATRIX(Sk);
	FREE_MATRIX(invPsik0);
	FREE_MATRIX(temp1);
	FREE_MATRIX(temp2);
	FREE_MATRIX(matconst);
	FREE_MATRIX(invconst);
	FREE_MATRIX(L);
	FREE_VECTOR(Eig);
	FREE_VECTOR(Muk_vec);
	FREE_VECTOR(MY_vec);
	FREE_MATRIX(Muk);
	//printf(" Q value %lf\n", -res);

	return(-res);


}




double Mstep_Manly_Full(int p, int T, int n, int K, double *misc_double, double ***Y, double **la, double **nu, double **gamma, double ***invS, double ***Mu, double ***invPsi, double *detS, double *detPsi, double *tau, int Mu_type){

	int i,j,j1,j2,k,t,t1,t2,sum_index1, count, *index1, sum_index2, *index2;
	double *Q_value, Q_value0, min, min_value, eps, det, *sum_gamma, *gamma_k, **Psi, **S, **temp1, **temp2, **temp3, **temp4, **invPsik, *Eig, **L;
	double **Muk, **invSk, **matconst, **invconst, **MYi, **tMYi, ***MY;
	double *Eig1, *Eig2;
	double **L1, **L2, *Muk_vec, *MY_vec;

	MAKE_VECTOR(Muk_vec, p+T-1);
	MAKE_VECTOR(MY_vec, p+T-1);
	MAKE_MATRIX(matconst, p+T-1, p+T-1);
	MAKE_MATRIX(invconst, p+T-1, p+T-1);
	MAKE_MATRIX(L, p+T-1, p+T-1);
	MAKE_VECTOR(Eig, p+T-1);


	MAKE_3ARRAY(MY, p,T,n);
	MAKE_MATRIX(MYi, p, T);
	MAKE_VECTOR(sum_gamma, K);
	MAKE_VECTOR(gamma_k, n);
	MAKE_VECTOR(Eig1, T);
	MAKE_VECTOR(Eig2, p);
	MAKE_VECTOR(index1, p);
	MAKE_VECTOR(index2, T);
	MAKE_MATRIX(Psi, T, T);
	MAKE_MATRIX(S, p, p);
	MAKE_MATRIX(tMYi, T, p);
	MAKE_MATRIX(temp1, p, T);
	MAKE_MATRIX(temp2, p, p);
	MAKE_MATRIX(temp3, T, p);
	MAKE_MATRIX(temp4, T, T);
	MAKE_MATRIX(L1, T, T);
	MAKE_MATRIX(L2, p, p);
	MAKE_MATRIX(invPsik, T, T);
	MAKE_MATRIX(invSk, p, p);
	MAKE_MATRIX(Muk, p, T);
	MAKE_VECTOR(Q_value, K);






	Anull(matconst, p+T-1, p+T-1);
	for(i=0; i<p; i++){	
		for(j=0; j<p; j++){

			if(i == j){
				matconst[i][j] = T*1.0;
			}

		}
	}

	for(i=p; i<p+T-1; i++){	
		for(j=p; j<p+T-1; j++){

			if(i == j){
				matconst[i][j] = p*1.0;
			}

		}
	}

	for(i=0; i<p; i++){	
		for(j=p; j<p+T-1; j++){
			matconst[i][j] = 1.0;

		}
	}
	for(i=p; i<p+T-1; i++){	
		for(j=0; j<p; j++){
			matconst[i][j] = 1.0;

		}
	}
	anull(Eig, p+T-1);

	#ifdef __HAVE_R_
		EigValDec(p+T-1, Eig, matconst, &det);
	#else
		cephes_symmeigens_down(p+T-1, Eig, matconst, &det);
	#endif

	Anull(L, p+T-1, p+T-1);

	for (t=0; t<p+T-1; t++){
		L[t][t] = 1.0 / Eig[t];
		
	}
	
	XAXt(matconst, p+T-1, L, invconst);



	eps = misc_double[0];

	anull(sum_gamma, K);
	Anull3(Mu, p, T, K);



	for(k=0; k<K; k++){

		for(i=0; i<n; i++){	
	
			sum_gamma[k] += gamma[i][k];
	
		}

		tau[k] = sum_gamma[k] / n; 

	}

	Q_value0 = 0;
	//optimize with respect to la;	
	for(k=0; k<K; k++){
		sum_index1 = 0;
	
		cpyv(gamma, k, n, gamma_k);		
		cpyk(invPsi, T, T, k, invPsik);


		for(j=0; j<p; j++){
			index1[j] = (la[k][j] != 0.0);
			sum_index1 += index1[j];
		}


		if(sum_index1 > 0){

			double *la_nonzero;

			MAKE_VECTOR(la_nonzero, sum_index1);
			count = 0;
			for(j=0; j<p; j++){
				if(index1[j] == 1){
					la_nonzero[count] = la[k][j];
					count += 1;
				}
			}
			


			min_value = simplex1(Q1, n, p, T, nu[k], index1, Y, gamma_k, invPsik, la_nonzero, eps, 0.1, Mu_type);

			count = 0;
			for(j=0; j<p; j++){
				if(index1[j] == 1){
					la[k][j] = la_nonzero[count];
					count += 1;
				}	
				
				else{
					la[k][j] = 0.0;
				}
			}			
			
			FREE_VECTOR(la_nonzero);
			Q_value[k] = min_value;
		
		} 


		else {
			double *la_nonzero;

			MAKE_VECTOR(la_nonzero, p);

			anull(la_nonzero, p);

			Q_value[k] = Q1(n, p, T, la_nonzero, nu[k], index1, Y, gamma_k, invPsik, Mu_type);


			FREE_VECTOR(la_nonzero);
		}

		Q_value0 += Q_value[k];
	}




	Q_value0 = 0;
	
	for(k=0; k<K; k++){
		sum_index2 = 0;
	
		cpyv(gamma, k, n, gamma_k);		
		cpyk(invPsi, T, T, k, invPsik);


		for(t=0; t<T; t++){
			index2[t] = (nu[k][t] != 0.0);
			sum_index2 += index2[t];
		}
	//printf(" Q  %d\n",sum_index2);



		if(sum_index2 > 0){

			double *nu_nonzero;

			MAKE_VECTOR(nu_nonzero, sum_index2);
			count = 0;
			for(t=0; t<T; t++){
				if(index2[t] == 1){
					nu_nonzero[count] = nu[k][t];
					count += 1;
				}
			}
			


			min_value = simplex2(Q2, n, p, T, la[k], index2, Y, gamma_k, invPsik, nu_nonzero, eps, 0.1, Mu_type);

			//printf(" nu_nonzero %lf %lf\n", nu_nonzero[0], nu_nonzero[1]);

			count = 0;
			for(t=0; t<T; t++){
				if(index2[t] == 1){
					nu[k][t] = nu_nonzero[count];
					count += 1;
				}	
				
				else{
					nu[k][t] = 0.0;
				}
			}			
			
			FREE_VECTOR(nu_nonzero);
			Q_value[k] = min_value;
		
		} 


		else {
			double *nu_nonzero;

			MAKE_VECTOR(nu_nonzero, T);

			anull(nu_nonzero, T);

			Q_value[k] = Q2(n, p, T, nu_nonzero, la[k], index2, Y, gamma_k, invPsik, Mu_type);


			FREE_VECTOR(nu_nonzero);
		}

		Q_value0 += Q_value[k];

	}


	for(k=0; k<K; k++){

		Manly_trans_whole(n, p, T, la[k], nu[k], Y, MY);

		if(Mu_type == 0){
			for(i=0; i<n; i++){

				cpyk(MY, p, T, i, MYi);
	
				for(j=0; j<p; j++){

					for(t=0; t<T; t++){
			
						Mu[j][t][k] += gamma[i][k] * MYi[j][t] / sum_gamma[k];


					}
				}

			}
		}
		else if(Mu_type == 1){

			anull(MY_vec, p+T-1);


	
			for(j=0; j<p; j++){
				for(i=0; i<n; i++){
					for(t=0; t<T; t++){		
						MY_vec[j] += gamma[i][k] * MY[j][t][i] / sum_gamma[k];
					}
				}
			}


			for(t=p; t<p+T-1; t++){
				for(i=0; i<n; i++){
					for(j=0; j<p; j++){		
						MY_vec[t] += gamma[i][k] * MY[j][t-p][i] / sum_gamma[k];
					}
				}
			}


			matxvec(invconst, p+T-1, p+T-1, MY_vec, p+T-1, Muk_vec);


			for(j=0; j<p; j++){
		
				for(t=0; t<T; t++){
			
					if(t<T-1){
						Mu[j][t][k] = Muk_vec[j] +  Muk_vec[p+t];

					}
					else{
						Mu[j][t][k] = Muk_vec[j];
					}
				}


			}
	

		}

		cpyk(Mu, p, T, k, Muk);
		cpyk(invPsi, T, T, k, invPsik);

		Anull(S, p, p);


		for(i=0; i<n; i++){

			cpyk(MY, p, T, i, MYi);

			mat_(p, T, MYi, Muk);
	
			tA(MYi, T, p, tMYi);


			multiply(MYi, p, T, invPsik, T, T, temp1);

			multiply(temp1, p, T, tMYi, T, p, temp2);

			for(j1=0; j1<p; j1++){			

				for(j2=0; j2<p; j2++){
	
					S[j1][j2] += gamma[i][k] * temp2[j1][j2] / T /sum_gamma[k];

				}
			}


		}



		anull(Eig2, p);

		#ifdef __HAVE_R_
			EigValDec(p, Eig2, S, &det);
		#else
			cephes_symmeigens_down(p, Eig2, S, &det);
		#endif


		min = INFINITY;
		for(j=0; j<p; j++){
			if(Eig2[j] < min){
				min = Eig2[j];
			}

		}


		if(min < pow(10, -300)){
			break;
		}else{
			Anull(L2, p, p);


			for (j=0; j<p; j++){
				L2[j][j] = 1.0 / Eig2[j];
			
			}
	
			XAXt(S, p, L2, invSk);

			detS[k] = det;

			cpyk2(invSk, p, p, invS, k);

		}



		Anull(Psi, T, T);


		for(i=0; i<n; i++){


			cpyk(MY, p, T, i, MYi);

			mat_(p, T, MYi, Muk);

			tA(MYi, T, p, tMYi);

			multiply(tMYi, T, p, invSk, p, p, temp3);

			multiply(temp3, T, p, MYi, p, T, temp4);

			for(t1=0; t1<T; t1++){			

				for(t2=0; t2<T; t2++){
	
					Psi[t1][t2] += gamma[i][k] * temp4[t1][t2] / p /sum_gamma[k];

				}
			}


		}


		anull(Eig1, T);

		#ifdef __HAVE_R_
			EigValDec(T, Eig1, Psi, &det);
		#else
			cephes_symmeigens_down(T, Eig1, Psi, &det);
		#endif

		min = INFINITY;
		for(t=0; t<T; t++){
			if(Eig1[t] < min){
				min = Eig1[t];
			}

		}

		if(min < pow(10, -300)){
			break;
		}else{
		
			Anull(L1, T, T);

			for (t=0; t<T; t++){
				L1[t][t] = 1.0 / Eig1[t];
		
			}
	
			XAXt(Psi, T, L1, invPsik);

			detPsi[k] = det;
			cpyk2(invPsik, T, T, invPsi, k);

		}



	}
	FREE_MATRIX(matconst);
	FREE_MATRIX(invconst);
	FREE_MATRIX(L);
	FREE_VECTOR(Eig);
	FREE_VECTOR(Muk_vec);
	FREE_VECTOR(MY_vec);

	FREE_3ARRAY(MY);
	FREE_VECTOR(Q_value);
	FREE_MATRIX(MYi);
	FREE_MATRIX(Muk);
	FREE_VECTOR(sum_gamma);
	FREE_VECTOR(gamma_k);
	FREE_VECTOR(Eig1);
	FREE_VECTOR(index1);
	FREE_VECTOR(index2);
	FREE_VECTOR(Eig2);
	FREE_MATRIX(Psi);	
	FREE_MATRIX(S);
	FREE_MATRIX(tMYi);
	FREE_MATRIX(temp1);
	FREE_MATRIX(temp2);
	FREE_MATRIX(temp3);
	FREE_MATRIX(temp4);
	FREE_MATRIX(L1);
	FREE_MATRIX(L2);
	FREE_MATRIX(invPsik);
	FREE_MATRIX(invSk);
	return Q_value0;
}




void EM_Manly_Full(int p, int T, int n, int K, double ***Y, double **la, double **nu, int max_iter, double *misc_double, double *tau, double ***Mu, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma, int *id, double *ll, int *conv, int Mu_type){
	int i,k,iter;
	double eps,loglik_old,loglik,max;

 	eps = misc_double[0];
	loglik_old = -INFINITY;
	iter = 0;

	do{
		loglik = loglik_old; 
		
		iter += 1;


		Estep_Manly_Full(p, T, n, K, Y, la, nu, tau, Mu, invS, invPsi, detS, detPsi, gamma);

		loglik_old = Mstep_Manly_Full(p, T, n, K, misc_double, Y, la, nu, gamma, invS, Mu, invPsi, detS, detPsi, tau, Mu_type);

							
	}

	while ((iter < max_iter) && (fabs(loglik - loglik_old) / fabs(loglik_old) > eps));

	ll[0] = mGloglik_Manly_Full(p, T, n, K, Y, la, nu, tau, Mu, invS, invPsi, detS, detPsi);
	
	conv[0] = iter;
	if(fabs(loglik - loglik_old) / fabs(loglik_old) <= eps){
		conv[1] = 0;
	} else{
		conv[1] = 1;
	}
	
	anulli(id, n);
	
	for(i=0; i<n; i++){
		max = -INFINITY;
		for(k=0; k<K; k++){
			if(gamma[i][k] > max){
				id[i] = k+1;
				max = gamma[i][k];
			}

		}

	}




}
