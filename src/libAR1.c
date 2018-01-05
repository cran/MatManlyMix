
#include "array.h"
#include "math.h"
#include "MatManlyMix.h" 


#define Inf 1e+140


#ifdef __HAVE_R_
	#include <R.h>
	#include <Rmath.h>
#endif


// Manly transformation for matrix Y
void Manly_trans_AR(int p, int T, double *la, double **Y, double **MY){

  int j, t;

  for (j=0; j<p; j++){
    if (fabs(la[j])<1e-12){
      for(t=0; t<T; t++) {
	MY[j][t] = Y[j][t];
      }
    } else{
      for(t=0; t<T; t++) {
	MY[j][t] = (exp(Y[j][t] * la[j]) - 1) / la[j];
      }
    }
      
  }

}

void Manly_trans_whole_AR(int n, int p, int T, double *la, double ***Y, double ***MY){

  int i, j, t;
 for(i=0; i<n; i++){


  for (j=0; j<p; j++){
    if (fabs(la[j])<1e-12){
      for(t=0; t<T; t++) {
	MY[j][t][i] = Y[j][t][i];
      }
    } else{
      for(t=0; t<T; t++) {
	MY[j][t][i] = (exp(Y[j][t][i] * la[j]) - 1) / la[j];
      }
    }
      
  }

 }
  
}


double mGpdf_Manly_AR(int p, int T, double *la, double **Y, double **Mu, double **invS, double **invPsi, double detS, double detPsi){

	int j, t;
	double trace = 0.0, dens, jac = 0.0;
	double **MY, **tMY, **temp1, **temp2, **maha;

	MAKE_MATRIX(maha, p, p);
	MAKE_MATRIX(temp1, p, T);
	MAKE_MATRIX(temp2, p, T);
	MAKE_MATRIX(tMY, T, p);	
	MAKE_MATRIX(MY, p, T);


	
  	Manly_trans_AR(p, T, la, Y, MY);
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

			jac = jac + la[j] * Y[j][t];
	
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






double mGloglik_Manly_AR(int p, int T, int n, int K, double ***Y, double **la, double *tau, double ***Mu, double ***invS, double ***invPsi, double *detS, double *detPsi){

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

			dens = mGpdf_Manly_AR(p, T, la[k], Yi, Muk, invSk, invPsik, detS[k], detPsi[k]);


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




void Estep_Manly_AR(int p, int T, int n, int K, double ***Y, double **la, double *tau, double ***Mu, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma){

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

				dens = mGpdf_Manly_AR(p, T, la[k], Yi, Muk, invSk, invPsik, detS[k], detPsi[k]);

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
double Q_AR(int n, int p, int T, double *la_nonzero, int *index, double ***Y, double *gamma_k, double **invPsik, int Mu_type){
  
	int i, j, j1, j2, t, count;
	double res, jac, sum_gamma, det1, det2, det;
	double *la;
	double **Sk, **temp3, **temp4;
	double ***MY, **Muk, *Eig1, *Eig2, **MYi, **tMYi, **invPsik0, *Muk_vec, *MY_vec, *Eig, **L, **matconst, **invconst;

	MAKE_VECTOR(Muk_vec, p+T-1);
	MAKE_VECTOR(MY_vec, p+T-1);
	MAKE_MATRIX(matconst, p+T-1, p+T-1);
	MAKE_MATRIX(invconst, p+T-1, p+T-1);
	MAKE_MATRIX(L, p+T-1, p+T-1);
	MAKE_VECTOR(Eig, p+T-1);


		
	MAKE_VECTOR(la, p);
	MAKE_VECTOR(Eig1, p);
	MAKE_VECTOR(Eig2, T);
	MAKE_3ARRAY(MY, p,T,n);
	MAKE_MATRIX(MYi, p, T);
	MAKE_MATRIX(tMYi, T, p);
	MAKE_MATRIX(Sk, p, p);

	MAKE_MATRIX(invPsik0, T, T);
	MAKE_MATRIX(temp3, p, T);
	MAKE_MATRIX(temp4, p, p);

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


	Manly_trans_whole_AR(n, p, T, la, Y, MY);
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

		multiply(MYi, p, T, invPsik0, T, T, temp3);

		multiply(temp3, p, T, tMYi, T, p, temp4);

		for(j1=0; j1<p; j1++){			

			for(j2=0; j2<p; j2++){
	
				Sk[j1][j2] += gamma_k[i] * temp4[j1][j2] / T /sum_gamma;

			}
		}


	}





	anull(Eig1, p);

	#ifdef __HAVE_R_
		EigValDec(p, Eig1, Sk, &det1);
	#else
		cephes_symmeigens_down(p, Eig1, Sk, &det1);
	#endif



	anull(Eig2, T);

	#ifdef __HAVE_R_
		EigValDec(T, Eig2, invPsik0, &det2);
	#else
		cephes_symmeigens_down(T, Eig2, invPsik0, &det2);
	#endif

	//printf("Q det1 det2 %lf %lf \n", det1, det2);


	res = -p / 2.0 * sum_gamma * log(1.0 / det2) -T / 2.0 * sum_gamma * log(det1); 

	
	for(i=0; i<n; i++){

		jac = 0;

		for(j=0; j<p; j++){

			for(t=0; t<T; t++){

				jac += Y[j][t][i] * la[j]; 
			}

		}


		res = res + gamma_k[i] * jac;

	}

	//printf("Q res %lf \n", res);


	FREE_VECTOR(la);
	FREE_MATRIX(MYi);
	FREE_3ARRAY(MY);
	FREE_MATRIX(tMYi);
	FREE_VECTOR(Eig1);
	FREE_VECTOR(Eig2);
	FREE_MATRIX(Sk);

	FREE_MATRIX(temp3);
	FREE_MATRIX(temp4);

	FREE_MATRIX(matconst);
	FREE_MATRIX(invconst);
	FREE_MATRIX(L);
	FREE_VECTOR(Eig);
	FREE_VECTOR(Muk_vec);
	FREE_VECTOR(MY_vec);


	FREE_MATRIX(Muk);

	FREE_MATRIX(invPsik0);

	return(-res);

}







double Mstep_Manly_AR(int p, int T, int n, int K, double *misc_double, double ***Y, double **la, double **gamma, double ***invS, double ***Mu, double ***invPsi, double *detS, double *detPsi, double *tau, int Mu_type){

	int i,j,j1,j2,k,t,t1,t2, sum_index, count, *index;
	double *Q_value, Q_value0, min, min_value, eps, det, *sum_gamma, *gamma_k, **Psi, **S, **temp1, **temp2, **temp3, **temp4, **invPsik, *Eig, **L;
	double **Muk, **invSk, **MYi, **tMYi, ***MY, *par, Psi2, phi, **A1, **A2, **I;
	double *Eig1, *Eig2;
	double **L1, **L2, *Muk_vec, *MY_vec, **matconst, **invconst;

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
	MAKE_VECTOR(index, p);

	MAKE_MATRIX(Psi, T, T);
	MAKE_MATRIX(S, p, p);
	MAKE_MATRIX(tMYi, T, p);
	MAKE_MATRIX(temp1, T, p);
	MAKE_MATRIX(temp2, T, T);
	MAKE_MATRIX(temp3, p, T);
	MAKE_MATRIX(temp4, p, p);
	MAKE_MATRIX(L1, T, T);
	MAKE_MATRIX(L2, p, p);
	MAKE_MATRIX(invPsik, T, T);
	MAKE_MATRIX(invSk, p, p);
	MAKE_MATRIX(Muk, p, T);
	MAKE_VECTOR(Q_value, K);
	MAKE_VECTOR(par,2);
	MAKE_MATRIX(A1, T, T);
	MAKE_MATRIX(A2, T, T);
	MAKE_MATRIX(I, T, T);


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
	anull(par, 2);

	anull(sum_gamma, K);
	Anull3(Mu, p, T, K);



	for(k=0; k<K; k++){

		for(i=0; i<n; i++){	
	
			sum_gamma[k] += gamma[i][k];
	
		}

		tau[k] = sum_gamma[k] / n; 

	}

	//optimize with respect to la;
	Q_value0 = 0;
	
	for(k=0; k<K; k++){
		sum_index = 0;
	
		cpyv(gamma, k, n, gamma_k);		
		cpyk(invPsi, T, T, k, invPsik);


		for(j=0; j<p; j++){
			index[j] = (la[k][j] != 0.0);
			sum_index += index[j];
		}


		if(sum_index > 0){

			double *la_nonzero;

			MAKE_VECTOR(la_nonzero, sum_index);
			count = 0;
			for(j=0; j<p; j++){
				if(index[j] == 1){
					la_nonzero[count] = la[k][j];
					count += 1;
				}
			}
			

			//Q_value = Q(n, p, T, la_nonzero, index, Y, gamma_k, invSk);

			min_value = simplex_AR(Q_AR, n, p, T, index, Y, gamma_k, invPsik, la_nonzero, eps, 0.1, Mu_type);

			count = 0;
			for(j=0; j<p; j++){
				if(index[j] == 1){
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

			Q_value[k] = Q_AR(n, p, T, la_nonzero, index, Y, gamma_k, invPsik, Mu_type);


			FREE_VECTOR(la_nonzero);
		}

		Q_value0 += Q_value[k];

	}


	

	for(k=0; k<K; k++){

		Manly_trans_whole_AR(n, p, T, la[k], Y, MY);

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

			multiply(MYi, p, T, invPsik, T, T, temp3);

			multiply(temp3, p, T, tMYi, T, p, temp4);

			for(j1=0; j1<p; j1++){			

				for(j2=0; j2<p; j2++){
	
					S[j1][j2] += gamma[i][k] * temp4[j1][j2] / T /sum_gamma[k];

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


		cpyv(gamma, k, n, gamma_k);	
	
		findPsi2phi(n, p, T, par, MY, Muk, gamma_k, invSk, eps/10.0);

		Psi2 = par[0];
		phi = par[1];

		Anull(A1, T, T);
		Anull(A2, T, T);
		Anull(I, T, T);
		Anull(invPsik, T, T);


		for(t1=0; t1<T; t1++){

			for(t2=0; t2<T; t2++){

				if(abs(t1-t2) == 1){

					A1[t1][t2] = - phi / Psi2;
 

				}

			}

		}

		for(t=1; t<T-1; t++){


			A2[t][t] = pow(phi, 2) / Psi2;
 

		}


		for(t=0; t<T; t++){


			I[t][t] = 1.0 / Psi2;
 

		}

		for(t1=0; t1<T; t1++){

			for(t2=0; t2<T; t2++){

				invPsik[t1][t2] = A1[t1][t2]+A2[t1][t2]+I[t1][t2];


			}

		}


		detPsi[k] = pow(Psi2, T) /  (1.0 - pow(phi, 2));
		cpyk2(invPsik, T, T, invPsi, k);


	}

	FREE_MATRIX(matconst);
	FREE_MATRIX(invconst);
	FREE_MATRIX(L);
	FREE_VECTOR(Eig);
	FREE_VECTOR(Muk_vec);
	FREE_VECTOR(MY_vec);



	FREE_VECTOR(par);

	FREE_3ARRAY(MY);
	FREE_VECTOR(Q_value);
	FREE_MATRIX(MYi);
	FREE_MATRIX(Muk);
	FREE_VECTOR(sum_gamma);
	FREE_VECTOR(gamma_k);
	FREE_VECTOR(Eig1);
	FREE_VECTOR(index);
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
	FREE_MATRIX(A1);
	FREE_MATRIX(A2);
	FREE_MATRIX(I);

	return Q_value0;
}

void findPsi2phi(int n, int p, int T, double *par, double ***MY, double **Muk, double *gamma_k, double **invSk, double eps){

	int t1, t2, i, t, j;

	double phi, Psi2, traceA1 = 0.0, traceA2 = 0.0, traceI = 0.0, sum_gamma = 0.0, trA1, trA2, trI, *coeff;
	
	double **A1, **A2, **I, **MYi, **tMYi, **temp1, **temp2, **temp3, **maha1, **maha2, **maha3, *try, *start;



	MAKE_MATRIX(A1, T, T);
	MAKE_MATRIX(A2, T, T);
	MAKE_MATRIX(I, T, T);
	MAKE_MATRIX(MYi, p, T);
	MAKE_MATRIX(tMYi, T, p);
	MAKE_MATRIX(temp1, p, T);
	MAKE_MATRIX(temp2, p, T);
	MAKE_MATRIX(temp3, p, T);
	MAKE_MATRIX(maha1, p, p);
	MAKE_MATRIX(maha2, p, p);
	MAKE_MATRIX(maha3, p, p);

	MAKE_VECTOR(try, 3);
	MAKE_VECTOR(coeff, 3);
	MAKE_VECTOR(start, 3);

	
 

	Psi2 = par[0];
	phi = par[1];
	
	for(i=0; i<n; i++){	
	
		sum_gamma += gamma_k[i];
	
	}
	
	Anull(A1, T, T);
	Anull(A2, T, T);
	Anull(I, T, T);

	//printf("Q psi2 phi %lf %lf \n", Psi2, phi);
	for(t1=0; t1<T; t1++){

		for(t2=0; t2<T; t2++){

			if(abs(t1-t2) == 1){

				A1[t1][t2] = 1;
 

			}

		}

	}

	for(t=1; t<T-1; t++){


		A2[t][t] = 1;
 

	}


	for(t=0; t<T; t++){


		I[t][t] = 1;
 

	}



	for(i=0; i<n; i++){


		cpyk(MY, p, T, i, MYi);

		mat_(p, T, MYi, Muk);

		tA(MYi, T, p, tMYi);

		multiply(invSk, p, p, MYi, p, T, temp1);

		multiply(temp1, p, T, A1, T, T, temp2);

		multiply(temp2, p, T, tMYi, T, p, maha1);

		multiply(temp1, p, T, A2, T, T, temp3);

		multiply(temp3, p, T, tMYi, T, p, maha2);

		multiply(temp1, p, T, tMYi, T, p, maha3);
		
		trA1 = 0;
		trA2 = 0;
		trI = 0;

		for(j=0; j<p; j++){			
	
			trA1 += maha1[j][j];
			trA2 += maha2[j][j];
			trI += maha3[j][j];
		}
		traceA1 += gamma_k[i] * trA1;
		traceA2 += gamma_k[i] * trA2;
		traceI += gamma_k[i] * trI;

	}

	coeff[0] = (1.0 - 2.0 / T) * traceA1 / (2.0 / T - 2) / traceA2;
	coeff[1] = (2 * traceA2 + 2.0 / T * traceI) / (2.0 / T - 2) / traceA2;
	coeff[2] = -traceA1 / (2.0 / T - 2) / traceA2;

	start[0] = 0.5;
	start[1] = 0;
	start[2] = -0.5;

	rootfinding(eq3, start, coeff, eps);

	
	try[0] = (-start[0] * traceA1 + pow(start[0], 2) * traceA2 + traceI) / p / T / sum_gamma;
	try[1] = (-start[1] * traceA1 + pow(start[1], 2) * traceA2 + traceI) / p / T / sum_gamma;
	try[2] = (-start[2] * traceA1 + pow(start[2], 2) * traceA2 + traceI) / p / T / sum_gamma;
	
	if(fabs(start[0])<1 && try[0]>0) {
		phi = start[0];
		Psi2 = try[0];

	}
	else if(fabs(start[1])<1 && try[1]>0) {
		phi = start[1];
		Psi2 = try[1];

	}
	else if(fabs(start[2])<1 && try[2]>0) {
		phi = start[2];
		Psi2 = try[2];

	}
	par[0] = Psi2;
	par[1] = phi;
	

	FREE_VECTOR(coeff);

	FREE_MATRIX(A1);
	FREE_MATRIX(A2);
	FREE_MATRIX(I);
	FREE_MATRIX(MYi);
	FREE_MATRIX(tMYi);
	FREE_MATRIX(temp1);
	FREE_MATRIX(temp2);
	FREE_MATRIX(temp3);
	FREE_MATRIX(maha1);
	FREE_MATRIX(maha2);
	FREE_VECTOR(start);
	FREE_MATRIX(maha3);

	FREE_VECTOR(try);


}

double eq3(double x, double *coeff)
{
    // the function we are interested in

    return pow(x,3) + coeff[0] * pow(x,2) + coeff[1] * x + coeff[2];
}




void rootfinding(double (*func)(double, double *), double *start, double *coeff, double eps){

	int i=0;
	int max_iterations = 1000;
	int done = 0;

	double p0,q0,r0,p,q,r;

	p = start[0];
	q = start[1];
	r = start[2];

	while (i<max_iterations && done == 0){   
		p0 = p;
		q0 = q;
		r0 = r;


		p = p0 - func(p0, coeff)/((p0-q0)*(p0-r0));
		q = q0 - func(q0, coeff)/((q0-p)*(q0-r0));
		r = r0 - func(r0, coeff)/((r0-p)*(r0-q));


    		if (fabs(p-p0)<eps && fabs(q-q0)<eps && fabs(r-r0)<eps){
        		done = 1;
			start[0] = p;
			start[1] = q;
			start[2] = r;

		}
		i++;
	}

}








void EM_Manly_AR(int p, int T, int n, int K, double ***Y, double **la, int max_iter, double *misc_double, double *tau, double ***Mu, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma, int *id, double *ll, int *conv, int Mu_type){
	int i,k,iter;
	double eps,loglik_old,loglik,max;

 	eps = misc_double[0];
	loglik_old = -INFINITY;
	iter = 0;

	do{
		loglik = loglik_old; 
		
		iter += 1;


		Estep_Manly_AR(p, T, n, K, Y, la, tau, Mu, invS, invPsi, detS, detPsi, gamma);

		loglik_old = Mstep_Manly_AR(p, T, n, K, misc_double, Y, la, gamma, invS, Mu, invPsi, detS, detPsi, tau, Mu_type);

							
	}

	while ((iter < max_iter) && (fabs(loglik - loglik_old) / fabs(loglik_old) > eps));

	ll[0] = mGloglik_Manly_AR(p, T, n, K, Y, la, tau, Mu, invS, invPsi, detS, detPsi);
	
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

