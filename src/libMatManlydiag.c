

#include "array.h"
#include "math.h"
#include "MatManlyMix.h"
 
#define Inf 1e+140


#ifdef __HAVE_R_
	#include <R.h>
	#include <Rmath.h>
#endif


double mGpdf_Manly_diag(int p, int T, double *la, double **Y, double **Mu, double **invS, double Psi2, double detS){

	int j, t;
	double trace = 0.0, dens, jac = 0.0;
	double **MY, **tMY, **temp1, **maha;

	MAKE_MATRIX(maha, p, p);
	MAKE_MATRIX(temp1, p, T);
	MAKE_MATRIX(tMY, T, p);	
	MAKE_MATRIX(MY, p, T);


	
  	Manly_trans(p, T, la, Y, MY);
	mat_(p, T, MY, Mu);

	tA(MY, T, p, tMY);

	multiply(invS, p, p, MY, p, T, temp1);

	multiply(temp1, p, T, tMY, T, p, maha);	


	for(j=0; j<p; j++){	 

		trace += maha[j][j];
	}

	
	dens = 1.0 / pow((2*PI), p*T/2.0) / pow(detS, T/2.0) / pow(Psi2, p*T/2.0) * exp(-1.0 / 2.0 * trace / Psi2);


	for(j=0; j<p; j++){

		for(t=0; t<T; t++){

			jac = jac + la[j] * Y[j][t];
	
		}	
	}

	dens = dens * exp(jac);


	FREE_MATRIX(maha);
	FREE_MATRIX(temp1);
	FREE_MATRIX(tMY);
	FREE_MATRIX(MY);
	return dens;
}




double mGloglik_Manly_diag(int p, int T, int n, int K, double ***Y, double **la, double *tau, double ***Mu, double ***invS, double *Psi2, double *detS){

	int i,k;
	double loglik = 0.0;
	double ll, dens = 0.0;
	double **Yi, **Muk, **invSk;
	
	MAKE_MATRIX(Yi, p, T);
	MAKE_MATRIX(Muk, p, T);
	MAKE_MATRIX(invSk, p, p);


	for(i=0; i<n; i++){	 

		ll = 0;

		for(k=0; k<K; k++){

			cpyk(Y, p, T, i, Yi);
			cpyk(Mu, p, T, k, Muk);
			cpyk(invS, p, p, k, invSk);

			dens = mGpdf_Manly_diag(p, T, la[k], Yi, Muk, invSk, Psi2[k], detS[k]);


			ll += tau[k] * dens;


		}

		loglik += log(ll);
	}

	
	FREE_MATRIX(Yi);
	FREE_MATRIX(Muk);
	FREE_MATRIX(invSk);

	return loglik;
}



void Estep_Manly_diag(int p, int T, int n, int K, double ***Y, double **la, double *tau, double ***Mu, double ***invS, double *Psi2, double *detS, double **gamma){

	int i,k;		
	double dens =0.0, *sum_gamma, **Yi, **Muk, **invSk;

	MAKE_MATRIX(Yi, p, T);
	MAKE_MATRIX(Muk, p, T);
	MAKE_MATRIX(invSk, p, p);

	MAKE_VECTOR(sum_gamma, n);

	anull(sum_gamma, n);

	//printf(" E step %lf %lf\n", Psi2[0], Psi2[1]);


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

				dens = mGpdf_Manly_diag(p, T, la[k], Yi, Muk, invSk, Psi2[k], detS[k]);

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

	FREE_VECTOR(sum_gamma);

}








// Q-function
double Q_diag(int n, int p, int T, double *la_nonzero, int *index, double ***Y, double *gamma_k, double **invSk){
  
	int i, j, t, count;
	double res, jac, sum_gamma, det1, det2, trace;
	double *la;
	double Psi2 = 0.0, **temp1;
	double ***MY, **Muk, *Eig2, **MYi, **tMYi, **invSk0, **maha;

		
	MAKE_VECTOR(la, p);

	MAKE_VECTOR(Eig2, p);
	MAKE_3ARRAY(MY, p,T,n);
	MAKE_MATRIX(MYi, p, T);
	MAKE_MATRIX(tMYi, T, p);
	MAKE_MATRIX(maha, p, p);

	MAKE_MATRIX(invSk0, p, p);
	MAKE_MATRIX(temp1, p, T);

	MAKE_MATRIX(Muk, p, T);

	count = 0;


	cpy(invSk, p, p, invSk0);

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


	Manly_trans_whole(n, p, T, la, Y, MY);
	res = 0;
	Anull(Muk, p, T);
	for(i=0; i<n; i++){


		for(j=0; j<p; j++){

			for(t=0; t<T; t++){
			
				Muk[j][t] += gamma_k[i] * MY[j][t][i] / sum_gamma;


			}
		}
	}
	
	
	for(i=0; i<n; i++){

		cpyk(MY, p, T, i, MYi);

		mat_(p, T, MYi, Muk);
	
		tA(MYi, T, p, tMYi);

		multiply(invSk0, p, p, MYi, p, T, temp1);

		multiply(temp1, p, T, tMYi, T, p, maha);
		trace = 0;	

		for(j=0; j<p; j++){	 

			trace += maha[j][j];
		}


	
		Psi2 += gamma_k[i] * trace / p / T / sum_gamma;



	}
	
	//printf(" Q func %lf %lf\n", Psi2);

	det1 = pow(Psi2, T);

	anull(Eig2, p);

	#ifdef __HAVE_R_
		EigValDec(p, Eig2, invSk0, &det2);
	#else
		cephes_symmeigens_down(p, Eig2, invSk0, &det2);
	#endif



	res = -p / 2.0 * sum_gamma * log(det1) -T / 2.0 * sum_gamma * log(1.0 / det2); 

	
	for(i=0; i<n; i++){

		jac = 0;

		for(j=0; j<p; j++){

			for(t=0; t<T; t++){

				jac += Y[j][t][i] * la[j]; 
			}

		}


		res = res + gamma_k[i] * jac;
	}



	FREE_VECTOR(la);
	FREE_MATRIX(MYi);
	FREE_3ARRAY(MY);
	FREE_MATRIX(tMYi);

	FREE_VECTOR(Eig2);
	FREE_MATRIX(maha);
	FREE_MATRIX(temp1);



	FREE_MATRIX(Muk);

	FREE_MATRIX(invSk0);

	return(-res);

}







double Mstep_Manly_diag(int p, int T, int n, int K, double *misc_double, double ***Y, double **la, double **gamma, double ***invS, double ***Mu, double *Psi2, double *detS, double *tau){

	int i,j,j1,j2,k,t, sum_index, count, *index;
	double *Q_value, Q_value0, min, min_value, trace, eps, det, *sum_gamma, *gamma_k, **S, **temp1, **temp4, **invPsik;
	double **Muk, **invSk, **maha, **MYi, **tMYi, ***MY;
	double *Eig1, *Eig2;
	double **L1, **L2;

	MAKE_MATRIX(MYi, p, T);
	MAKE_VECTOR(sum_gamma, K);
	MAKE_VECTOR(gamma_k, n);
	MAKE_VECTOR(Eig1, T);
	MAKE_VECTOR(Eig2, p);
	MAKE_VECTOR(index, p);
	MAKE_MATRIX(S, p, p);
	MAKE_MATRIX(maha, p, p);
	MAKE_MATRIX(tMYi, T, p);
	MAKE_MATRIX(temp1, p, T);
	MAKE_MATRIX(temp4, p, p);
	MAKE_MATRIX(L1, T, T);
	MAKE_MATRIX(L2, p, p);
	MAKE_MATRIX(invPsik, T, T);
	MAKE_MATRIX(invSk, p, p);
	MAKE_MATRIX(Muk, p, T);
	MAKE_VECTOR(Q_value, K);
	MAKE_3ARRAY(MY, p,T,n);

	eps = misc_double[0];

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
		cpyk(invS, p, p, k, invSk);


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

			min_value = simplex(Q_diag, n, p, T, index, Y, gamma_k, invSk, la_nonzero, eps, 0.1);

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

			Q_value[k] = Q_diag(n, p, T, la_nonzero, index, Y, gamma_k, invSk);


			FREE_VECTOR(la_nonzero);
		}

		Q_value0 += Q_value[k];

	}

	anull(Psi2, K);

	for(k=0; k<K; k++){


		Manly_trans_whole(n, p, T, la[k], Y, MY);


		for(i=0; i<n; i++){

			cpyk(MY, p, T, i, MYi);

			for(j=0; j<p; j++){

				for(t=0; t<T; t++){
			
					Mu[j][t][k] += gamma[i][k] * MYi[j][t] / sum_gamma[k];


				}
			}

		}



		cpyk(Mu, p, T, k, Muk);
		cpyk(invS, p, p, k, invSk);


		for(i=0; i<n; i++){


 			cpyk(MY, p, T, i, MYi);

			mat_(p, T, MYi, Muk);
	
			tA(MYi, T, p, tMYi);

			multiply(invSk, p, p, MYi, p, T, temp1);

			multiply(temp1, p, T, tMYi, T, p, maha);	

			trace = 0;
			for(j=0; j<p; j++){	 

				trace += maha[j][j];
			}


	
			Psi2[k] += gamma[i][k] * trace / p / T / sum_gamma[k];


		}


		//printf(" M step %lf \n", Psi2[k]);


		Anull(S, p, p);


		for(i=0; i<n; i++){


			cpyk(MY, p, T, i, MYi);

			mat_(p, T, MYi, Muk);

			tA(MYi, T, p, tMYi);

			multiply(MYi, p, T, tMYi, T, p, temp4);

			for(j1=0; j1<p; j1++){			

				for(j2=0; j2<p; j2++){
	
					S[j1][j2] += gamma[i][k] * temp4[j1][j2] / T /sum_gamma[k] / Psi2[k];

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


	}

	FREE_MATRIX(maha);
	FREE_VECTOR(Q_value);
	FREE_MATRIX(MYi);
	FREE_MATRIX(Muk);
	FREE_VECTOR(sum_gamma);
	FREE_VECTOR(gamma_k);
	FREE_VECTOR(Eig1);
	FREE_VECTOR(index);
	FREE_VECTOR(Eig2);	
	FREE_MATRIX(S);
	FREE_MATRIX(tMYi);
	FREE_MATRIX(temp1);
	FREE_MATRIX(temp4);
	FREE_MATRIX(L1);
	FREE_MATRIX(L2);
	FREE_MATRIX(invPsik);
	FREE_MATRIX(invSk);
	FREE_3ARRAY(MY);

	return Q_value0;
}




void EM_Manly_diag(int p, int T, int n, int K, double ***Y, double **la, int max_iter, double *misc_double, double *tau, double ***Mu, double ***invS, double *Psi2, double *detS, double **gamma, int *id, double *ll, int *conv){
	int i,k,iter;
	double eps,loglik_old,loglik,max;

 	eps = misc_double[0];
	loglik_old = -INFINITY;
	iter = 0;

	do{
		loglik = loglik_old; 
		
		iter += 1;


		Estep_Manly_diag(p, T, n, K, Y, la, tau, Mu, invS, Psi2, detS, gamma);


		loglik_old = Mstep_Manly_diag(p, T, n, K, misc_double, Y, la, gamma, invS, Mu, Psi2, detS, tau);
			
						
	}

	while ((iter < max_iter) && (fabs(loglik - loglik_old) / fabs(loglik_old) > eps));

	ll[0] = mGloglik_Manly_diag(p, T, n, K, Y, la, tau, Mu, invS, Psi2, detS);

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





double mGpdf_Manly_diag_Reg(int p, int T, int q, double *la, double **Y, double **X, double **beta, double **invS, double Psi2, double detS){

	int j, t;
	double trace = 0.0, dens, jac = 0.0;
	double **MY, **tMY, **temp1, **maha, **Mu, **tMu;

	MAKE_MATRIX(maha, p, p);
	MAKE_MATRIX(temp1, p, T);
	MAKE_MATRIX(tMY, T, p);	
	MAKE_MATRIX(MY, p, T);
	MAKE_MATRIX(Mu, p, T);
	MAKE_MATRIX(tMu, T, p);

	

	multiply(X, T, q, beta, q, p, tMu);

	tA(tMu, p, T, Mu);

  	Manly_trans(p, T, la, Y, MY);

	mat_(p, T, MY, Mu);

	tA(MY, T, p, tMY);

	multiply(invS, p, p, MY, p, T, temp1);
	
	multiply(temp1, p, T, tMY, T, p, maha);	


	for(j=0; j<p; j++){	 

		trace += maha[j][j];
	}
	
	dens = 1.0 / pow((2*PI), p*T/2.0) / pow(detS, T/2.0) / pow(Psi2, p*T/2.0) * exp(-1.0 / 2.0 * trace / Psi2);

	for(j=0; j<p; j++){

		for(t=0; t<T; t++){

			jac = jac + la[j] * Y[j][t];
	
		}	
	}

	dens = dens * exp(jac);

	FREE_MATRIX(tMu);
	FREE_MATRIX(maha);
	FREE_MATRIX(temp1);
	FREE_MATRIX(tMY);
	FREE_MATRIX(MY);
	FREE_MATRIX(Mu);
	return dens;
}




double mGloglik_Manly_diag_Reg(int p, int T, int n, int q, int K, double ***Y, double **la, double *tau, double ***X, double ***beta, double ***invS, double *Psi2, double *detS){

	int i,k;
	double loglik = 0.0;
	double ll, dens = 0.0;
	double **Yi, **Xi, **betak, **tMuk, **Muk, **invSk;
	
	MAKE_MATRIX(Yi, p, T);
	MAKE_MATRIX(Xi, T, q);
	MAKE_MATRIX(betak, q, p);
	MAKE_MATRIX(invSk, p, p);
	MAKE_MATRIX(Muk, T, p);
	MAKE_MATRIX(tMuk, p, T);


	for(i=0; i<n; i++){	 

		ll = 0;

		for(k=0; k<K; k++){

			cpyk(Y, p, T, i, Yi);

			cpyk(X, T, q, i, Xi);

			cpyk(beta, q, p, k, betak);

			cpyk(invS, p, p, k, invSk);

			dens = mGpdf_Manly_diag_Reg(p, T, q, la[k], Yi, Xi, betak, invSk, Psi2[k], detS[k]);

			ll += tau[k] * dens;


		}

		loglik += log(ll);
	}

	FREE_MATRIX(Yi);
	FREE_MATRIX(Xi);
	FREE_MATRIX(betak);
	FREE_MATRIX(invSk);
	FREE_MATRIX(Muk);
	FREE_MATRIX(tMuk);

	return loglik;
}



void Estep_Manly_diag_Reg(int p, int T, int n, int q, int K, double ***Y, double **la, double *tau, double ***X, double ***beta, double ***invS, double *Psi2, double *detS,  double **gamma){

	int i,k;		
	double dens =0.0, *sum_gamma, **betak, **Yi, **Xi, **invSk;

	MAKE_MATRIX(Yi, p, T);
	MAKE_MATRIX(Xi, T, q);
	MAKE_MATRIX(betak, q, p);
	MAKE_MATRIX(invSk, p, p);
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

				cpyk(X, T, q, i, Xi);

				cpyk(beta, q, p, k, betak);


				cpyk(invS, p, p, k, invSk);

				dens = mGpdf_Manly_diag_Reg(p, T, q, la[k], Yi, Xi, betak, invSk, Psi2[k], detS[k]);

				gamma[i][k] = tau[k] * dens;

				sum_gamma[i] += gamma[i][k];


			}

	
			for(k=0; k<K; k++){

				gamma[i][k] = gamma[i][k] / sum_gamma[i];

			}


		}	 

	}




	FREE_MATRIX(Yi);
	FREE_MATRIX(Xi);
	FREE_MATRIX(betak);
	FREE_MATRIX(invSk);

	FREE_VECTOR(sum_gamma);

}








// Q-function
double Q_diag_Reg(int n, int p, int T, int q, double *la_nonzero, int *index, double ***Y, double *gamma_k, double **invSk, double ***X){
  
	int i, j, l1, l2,l, t, count;
	double res, jac, sum_gamma, det1, det2, det;
	double *la;
	double ***MY, *Eig2, *Eig, **maha, **MYi, **Xi, **tXi, **tMYi, **temp1, **temp2, **temp3, Psi2, trace, **invSk0;
	double **betak, **L, **inv_part, **part2, **part1, **Muk, **tMuk;

		
	MAKE_VECTOR(la, p);
	MAKE_VECTOR(Eig, q);

	MAKE_VECTOR(Eig2, p);
	MAKE_3ARRAY(MY, p,T,n);
	MAKE_MATRIX(MYi, p, T);
	MAKE_MATRIX(tMYi, T, p);
	MAKE_MATRIX(Xi, T, q);
	MAKE_MATRIX(tXi, q, T);
	MAKE_MATRIX(L, q, q);

	MAKE_MATRIX(inv_part, q, q);
	MAKE_MATRIX(betak, q, p);
	MAKE_MATRIX(part2, q, p);
	MAKE_MATRIX(part1, q, q);
	MAKE_MATRIX(Muk, p, T);
	MAKE_MATRIX(tMuk, T, p);

	MAKE_MATRIX(maha, p, p);
	MAKE_MATRIX(invSk0, p, p);

	MAKE_MATRIX(temp1, p, T);

	MAKE_MATRIX(temp2, q, q);
	MAKE_MATRIX(temp3, q, p);



	count = 0;

	cpy(invSk, p, p, invSk0);

	for(j=0; j<p; j++){
		if(index[j]==1){
			la[j] = la_nonzero[count];
			count += 1;
		}
		else{
			la[j] = 0; 

		}
	}	

	//printf(" Q la %lf %lf\n",la[0], la[1]);
	sum_gamma = 0;
			
	for(i=0; i<n; i++){


		sum_gamma += gamma_k[i];
	}


	Manly_trans_whole(n, p, T, la, Y, MY);

	res = 0;

	Anull(inv_part, q, q);
	Anull(part2, q, p);

	for(i=0; i<n; i++){

	
		cpyk(X, T, q, i, Xi);

		tA(Xi, q, T, tXi);

		multiply(tXi, q, T, Xi, T, q, temp2);

		cpyk(MY, p, T, i, MYi);

		tA(MYi, T, p, tMYi);


		multiply(tXi, q, T, tMYi, T, p, temp3);

		for(l1=0; l1<q; l1++){

			for(l2=0; l2<q; l2++){

				inv_part[l1][l2] += gamma_k[i] * temp2[l1][l2]; 

			}
		}


		for(l=0; l<q; l++){

			for(j=0; j<p; j++){

				part2[l][j] += gamma_k[i] * temp3[l][j]; 

			}
		}
	}


	anull(Eig, q);

	#ifdef __HAVE_R_
		EigValDec(q, Eig, inv_part, &det);
	#else
		cephes_symmeigens_down(q, Eig, inv_part, &det);
	#endif

	Anull(L, q, q);

	for (j=0; j<q; j++){
		L[j][j] = 1.0 / Eig[j];
		
	}
	
	XAXt(inv_part, q, L, part1);

	multiply(part1, q, q, part2, q, p, betak);

	Psi2 = 0;
	for(i=0; i<n; i++){

		cpyk(MY, p, T, i, MYi);

		cpyk(X, T, q, i, Xi);

		multiply(Xi, T, q, betak, q, p, tMuk);

		tA(tMuk, p, T, Muk);

		mat_(p, T, MYi, Muk);
	
		tA(MYi, T, p, tMYi);

		multiply(invSk0, p, p, MYi, p, T, temp1);

		multiply(temp1, p, T, tMYi, T, p, maha);

		trace = 0;	

		for(j=0; j<p; j++){	 

			trace += maha[j][j];
		}


	
		Psi2 += gamma_k[i] * trace / p / T / sum_gamma;



	}

	//printf(" Q  %lf\n",Psi2);


	
	det1 = pow(Psi2, T);
	anull(Eig2, p);

	#ifdef __HAVE_R_
		EigValDec(p, Eig2, invSk0, &det2);
	#else
		cephes_symmeigens_down(p, Eig2, invSk0, &det2);
	#endif

	//printf(" det2 %lf\n", det2);
	//printf(" det1 %lf\n", det1);


	res = -p / 2.0 * sum_gamma * log(det1) -T / 2.0 * sum_gamma * log(1.0 / det2); 

	
	for(i=0; i<n; i++){

		jac = 0;

		for(j=0; j<p; j++){

			for(t=0; t<T; t++){

				jac += Y[j][t][i] * la[j]; 
			}

		}


		res = res + gamma_k[i] * jac;
	}
	//printf(" Q  %lf\n",res);

	
	FREE_VECTOR(la);
	FREE_VECTOR(Eig);

	FREE_VECTOR(Eig2);
	FREE_3ARRAY(MY);
	FREE_MATRIX(MYi);
	FREE_MATRIX(tMYi);
	FREE_MATRIX(betak);
	FREE_MATRIX(L);
	FREE_MATRIX(maha);
	FREE_MATRIX(inv_part);
	FREE_MATRIX(part2);
	FREE_MATRIX(part1);
	FREE_MATRIX(Muk);
	FREE_MATRIX(tMuk);
	FREE_MATRIX(Xi);
	FREE_MATRIX(tXi);
	FREE_MATRIX(temp2);
	FREE_MATRIX(temp3);
	FREE_MATRIX(temp1);

	return(-res);

}







double Mstep_Manly_diag_Reg(int p, int T, int n, int q, int K, double *misc_double, double ***Y, double **la, double **gamma, double ***invS, double ***X, double ***beta, double *Psi2, double *detS, double *tau){

	int i,j,j1,j2,l1, l2, l, k,sum_index, count, *index;
	double *Q_value, Q_value0, min, min_value, eps, det, *sum_gamma, *gamma_k, **S,  **temp2, **temp3, **temp5, **temp6, **maha;
	double trace, ***MY, ***Muk, **Muki, **tMuki, **betak, **part1, **part2, **inv_part, **invSk, **MYi, **tMYi, **Xi, **tXi;
	double  *Eig2, *Eig;
	double **L2, **L;

	MAKE_3ARRAY(MY, p,T,n);	
	MAKE_MATRIX(MYi, p, T);
	MAKE_VECTOR(sum_gamma, K);
	MAKE_VECTOR(gamma_k, n);
	MAKE_VECTOR(Eig, q);

	MAKE_VECTOR(Eig2, p);
	MAKE_VECTOR(index, p);
	MAKE_MATRIX(S, p, p);
	MAKE_MATRIX(tMYi, T, p);
	MAKE_MATRIX(temp2, q, q);
	MAKE_MATRIX(temp3, q, p);
	MAKE_MATRIX(maha, p, p);

	MAKE_MATRIX(temp5, p, T);
	MAKE_MATRIX(temp6, p, p);
	MAKE_MATRIX(L, q, q);

	MAKE_MATRIX(L2, p, p);
	MAKE_MATRIX(Xi, T, q);
	MAKE_MATRIX(tXi, q, T);
	MAKE_MATRIX(invSk, p, p);

	MAKE_MATRIX(tMuki, T, p);
	MAKE_MATRIX(Muki, p, T);
	MAKE_3ARRAY(Muk, p,T,n);	

	MAKE_VECTOR(Q_value, K);
	MAKE_MATRIX(inv_part, q, q);
	MAKE_MATRIX(betak, q, p);
	MAKE_MATRIX(part2, q, p);
	MAKE_MATRIX(part1, q, q);


	eps = misc_double[0];

	anull(sum_gamma, K);


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
		cpyk(invS, p, p, k, invSk);


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
			

			min_value = simplex_Reg(Q_diag_Reg, n, p, T, q, index, Y, gamma_k, invSk, X, la_nonzero, eps, 0.1);

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

			Q_value[k] = Q_diag_Reg(n, p, T, q, la_nonzero, index, Y, gamma_k, invSk, X);

			FREE_VECTOR(la_nonzero);
		}

		Q_value0 += Q_value[k];

	}



	anull(Psi2, K);


	for(k=0; k<K; k++){


		Anull(inv_part, q, q);
		Anull(part2, q, p);

		cpyv(gamma, k, n, gamma_k);		

		Manly_trans_whole(n, p, T, la[k], Y, MY);


		for(i=0; i<n; i++){

			cpyk(MY, p, T, i, MYi);

			tA(MYi, T, p, tMYi);

			cpyk(X, T, q, i, Xi);

			tA(Xi, q, T, tXi);

			multiply(tXi, q, T, Xi, T, q, temp2);

			multiply(tXi, q, T, tMYi, T, p, temp3);
	
			for(l1=0; l1<q; l1++){

				for(l2=0; l2<q; l2++){

					inv_part[l1][l2] += gamma_k[i] * temp2[l1][l2]; 

				}
			}


			for(l=0; l<q; l++){

				for(j=0; j<p; j++){

					part2[l][j] += gamma_k[i] * temp3[l][j]; 

				}
			}
		}

		anull(Eig, q);

		#ifdef __HAVE_R_
			EigValDec(q, Eig, inv_part, &det);
		#else
			cephes_symmeigens_down(q, Eig, inv_part, &det);
		#endif

		Anull(L, q, q);

		for (j=0; j<q; j++){
			L[j][j] = 1.0 / Eig[j];
		
		}
	
		XAXt(inv_part, q, L, part1);

		multiply(part1, q, q, part2, q, p, betak);

		cpyk2(betak, q, p, beta, k);





		for(i=0; i<n; i++){

			cpyk(X, T, q, i, Xi);

			multiply(Xi, T, q, betak, q, p, tMuki);

			tA(tMuki, p, T, Muki);

			cpyk2(Muki, p, T, Muk, i);


		}

		cpyk(invS, p, p, k, invSk);

	//printf(" M step invSk %d %d %lf %lf %lf %lf\n", k,K, invSk[0][0], invSk[0][1], invSk[1][0], invSk[1][1]);
	//printf(" M step betak %d %d %lf %lf %lf %lf\n", k,K, betak[0][0], betak[0][1], betak[1][0], betak[1][1]);

		for(i=0; i<n; i++){

			cpyk(MY, p, T, i, MYi);

			cpyk(Muk, p, T, i, Muki);

			mat_(p, T, MYi, Muki);
	
			tA(MYi, T, p, tMYi);



			multiply(invSk, p, p, MYi, p, T, temp5);

			multiply(temp5, p, T, tMYi, T, p, maha);
			trace = 0;	

			for(j=0; j<p; j++){	 

				trace += maha[j][j];
			}


	
			Psi2[k] += gamma[i][k] * trace / p / T / sum_gamma[k];



		}


//printf(" M step Psi2 %d %d %lf\n", k,K, Psi2[k]);


		Anull(S, p, p);


		for(i=0; i<n; i++){


			cpyk(MY, p, T, i, MYi);

			cpyk(Muk, p, T, i, Muki);

			mat_(p, T, MYi, Muki);

			tA(MYi, T, p, tMYi);

			multiply(MYi, p, T, tMYi, T, p, temp6);

			for(j1=0; j1<p; j1++){			

				for(j2=0; j2<p; j2++){
	
					S[j1][j2] += gamma[i][k] * temp6[j1][j2] / Psi2[k] / T /sum_gamma[k];

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

		
	}


	FREE_3ARRAY(MY);	
	FREE_VECTOR(Eig);
	FREE_VECTOR(Q_value);
	FREE_MATRIX(MYi);
	FREE_MATRIX(tMuki);
	FREE_MATRIX(Muki);
	FREE_3ARRAY(Muk);
	FREE_VECTOR(sum_gamma);
	FREE_VECTOR(gamma_k);
	FREE_VECTOR(index);
	FREE_VECTOR(Eig2);	
	FREE_MATRIX(S);
	FREE_MATRIX(tMYi);
	FREE_MATRIX(temp2);
	FREE_MATRIX(temp3);
	FREE_MATRIX(temp5);
	FREE_MATRIX(temp6);
	FREE_MATRIX(maha);


	FREE_MATRIX(L2);
	FREE_MATRIX(L);
	FREE_MATRIX(invSk);
	FREE_MATRIX(inv_part);
	FREE_MATRIX(betak);
	FREE_MATRIX(part2);
	FREE_MATRIX(part1);
	FREE_MATRIX(Xi);
	FREE_MATRIX(tXi);

	return Q_value0;
}




void EM_Manly_diag_Reg(int p, int T, int n, int q, int K, double ***Y, double **la, int max_iter, double *misc_double, double *tau, double ***X, double ***beta, double ***invS, double *Psi2, double *detS, double **gamma, int *id, double *ll, int *conv){
	int i,k,iter;
	double eps,loglik_old,loglik,max;

 	eps = misc_double[0];
	loglik_old = -INFINITY;
	iter = 0;

	do{
		loglik = loglik_old; 
		
		iter += 1;


		Estep_Manly_diag_Reg(p, T, n, q, K, Y, la, tau, X, beta, invS, Psi2, detS, gamma);

		loglik_old = Mstep_Manly_diag_Reg(p, T, n, q, K, misc_double, Y, la, gamma, invS, X, beta, Psi2, detS, tau);
						
	}

	while ((iter < max_iter) && (fabs(loglik - loglik_old) / fabs(loglik_old) > eps));

	ll[0] = mGloglik_Manly_diag_Reg(p, T, n, q, K, Y, la, tau, X, beta, invS, Psi2, detS);

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



