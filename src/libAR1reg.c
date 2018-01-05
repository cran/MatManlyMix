

#include "array.h"
#include "math.h"
#include "MatManlyMix.h" 


#define Inf 1e+140


#ifdef __HAVE_R_
	#include <R.h>
	#include <Rmath.h>
#endif


double mGpdf_Manly_AR_Reg(int p, int T, int q, double *la, double **Y, double **X, double **beta, double **invS, double **invPsi, double detS, double detPsi){

	int j, t;
	double trace = 0.0, dens, jac = 0.0;
	double **MY, **tMY, **temp1, **temp2, **maha, **Mu, **tMu;

	MAKE_MATRIX(maha, p, p);
	MAKE_MATRIX(temp1, p, T);
	MAKE_MATRIX(temp2, p, T);
	MAKE_MATRIX(tMY, T, p);	
	MAKE_MATRIX(MY, p, T);
	MAKE_MATRIX(Mu, p, T);
	MAKE_MATRIX(tMu, T, p);

	multiply(X, T, q, beta, q, p, tMu);

	tA(tMu, p, T, Mu);
	
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
	FREE_MATRIX(Mu);
	FREE_MATRIX(tMu);
	return dens;
}






double mGloglik_Manly_AR_Reg(int p, int T, int n, int q, int K, double ***Y, double **la, double *tau, double ***X, double ***beta, double ***invS, double ***invPsi, double *detS, double *detPsi){

	int i,k;
	double loglik = 0.0;
	double ll, dens = 0.0;
	double **Yi, **Xi, **betak, **invSk, **invPsik;
	
	MAKE_MATRIX(Yi, p, T);
	MAKE_MATRIX(Xi, T, q);
	MAKE_MATRIX(betak, q, p);
	MAKE_MATRIX(invSk, p, p);
	MAKE_MATRIX(invPsik, T, T);


	for(i=0; i<n; i++){	 

		ll = 0;

		for(k=0; k<K; k++){

			cpyk(Y, p, T, i, Yi);

			cpyk(X, T, q, i, Xi);

			cpyk(beta, q, p, k, betak);
			cpyk(invS, p, p, k, invSk);
			cpyk(invPsi, T, T, k, invPsik);

			dens = mGpdf_Manly_AR_Reg(p, T, q, la[k], Yi, Xi, betak, invSk, invPsik, detS[k], detPsi[k]);


			ll += tau[k] * dens;


		}

		loglik += log(ll);
	}

	FREE_MATRIX(Yi);
	FREE_MATRIX(Xi);
	FREE_MATRIX(betak);
	FREE_MATRIX(invSk);
	FREE_MATRIX(invPsik);

	return loglik;
}




void Estep_Manly_AR_Reg(int p, int T, int n, int q, int K, double ***Y, double **la, double *tau, double ***X, double ***beta, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma){

	int i,k;		
	double dens =0.0, *sum_gamma, **Yi, **betak, **Xi, **invSk, **invPsik;

	MAKE_MATRIX(Yi, p, T);
	MAKE_MATRIX(Xi, T, q);
	MAKE_MATRIX(betak, q, p);

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
				cpyk(X, T, q, i, Xi);

				cpyk(beta, q, p, k, betak);
				cpyk(invS, p, p, k, invSk);
				cpyk(invPsi, T, T, k, invPsik);

				dens = mGpdf_Manly_AR_Reg(p, T, q, la[k], Yi, Xi, betak, invSk, invPsik, detS[k], detPsi[k]);

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
	FREE_MATRIX(invPsik);

	FREE_VECTOR(sum_gamma);

}







// Q-function
double Q_AR_Reg(int n, int p, int T, int q, double *la_nonzero, int *index, double ***Y, double *gamma_k, double **invPsik, double ***X){
  
	int i, j, j1,j2,l1, l2,l, t, count;
	double res, jac, sum_gamma, det1, det2, det;
	double *la;
	double **Sk, **temp1, **temp2, **temp3, **temp4, **temp5;
	double ***MY, *Eig1, *Eig2, *Eig, **MYi, **Xi, **tXi, **tMYi, **invPsik0;
	double **betak, **L, **inv_part, **part2, **part1, **Muk, **tMuk;

		
	MAKE_VECTOR(la, p);
	MAKE_VECTOR(Eig, q);
	MAKE_VECTOR(Eig1, T);
	MAKE_VECTOR(Eig2, p);
	MAKE_3ARRAY(MY, p,T,n);
	MAKE_MATRIX(MYi, p, T);
	MAKE_MATRIX(tMYi, T, p);

	MAKE_MATRIX(tXi, q, T);
	MAKE_MATRIX(Sk, p, p);
	MAKE_MATRIX(L, q, q);
	MAKE_MATRIX(invPsik0, T, T);
	MAKE_MATRIX(temp1, q, T);
	MAKE_MATRIX(temp2, q, q);
	MAKE_MATRIX(temp3, q, p);
	MAKE_MATRIX(temp4, p, T);
	MAKE_MATRIX(temp5, p, p);
	MAKE_MATRIX(inv_part, q, q);
	MAKE_MATRIX(betak, q, p);
	MAKE_MATRIX(part2, q, p);
	MAKE_MATRIX(part1, q, q);
	MAKE_MATRIX(Muk, p, T);
	MAKE_MATRIX(tMuk, T, p);
	MAKE_MATRIX(Xi, T, q);
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

	//printf(" Q la %lf %lf\n",la[0], la[1]);
	sum_gamma = 0;
			
	for(i=0; i<n; i++){


		sum_gamma += gamma_k[i];
	}


	Manly_trans_whole_AR(n, p, T, la, Y, MY);

	res = 0;

	Anull(inv_part, q, q);
	Anull(part2, q, p);

	for(i=0; i<n; i++){

	
		cpyk(X, T, q, i, Xi);

		tA(Xi, q, T, tXi);

		multiply(tXi, q, T, invPsik0, T, T, temp1);

		multiply(temp1, q, T, Xi, T, q, temp2);

		cpyk(MY, p, T, i, MYi);

		tA(MYi, T, p, tMYi);


		multiply(temp1, q, T, tMYi, T, p, temp3);

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


	//printf(" Q beta %lf %lf\n",betak[0][0], betak[0][1]);
	//printf(" Q MY %lf %lf\n",MY[0][0][0], MY[0][1][0]);

	
	Anull(Sk, p, p);
	
	for(i=0; i<n; i++){


		cpyk(MY, p, T, i, MYi);

		cpyk(X, T, q, i, Xi);

		multiply(Xi, T, q, betak, q, p, tMuk);

		tA(tMuk, p, T, Muk);

		mat_(p, T, MYi, Muk);

	
		tA(MYi, T, p, tMYi);

		multiply(MYi, p, T, invPsik0, T, T, temp4);


		multiply(temp4, p, T, tMYi, T, p, temp5);

	
		for(j1=0; j1<p; j1++){			
	
			for(j2=0; j2<p; j2++){
	
				Sk[j1][j2] += gamma_k[i] * temp5[j1][j2] / T /sum_gamma;


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

				jac += Y[j][t][i] * la[j]; 
			}

		}


		res = res + gamma_k[i] * jac;
	}
	//printf(" Q  %lf\n",res);

	
	FREE_VECTOR(la);
	FREE_VECTOR(Eig);
	FREE_VECTOR(Eig1);
	FREE_VECTOR(Eig2);
	FREE_3ARRAY(MY);
	FREE_MATRIX(MYi);
	FREE_MATRIX(tMYi);
	FREE_MATRIX(Sk);
	FREE_MATRIX(betak);
	FREE_MATRIX(L);
	FREE_MATRIX(invPsik0);
	FREE_MATRIX(temp1);
	FREE_MATRIX(temp2);
	FREE_MATRIX(temp3);
	FREE_MATRIX(temp4);
	FREE_MATRIX(temp5);
	FREE_MATRIX(inv_part);
	FREE_MATRIX(part2);
	FREE_MATRIX(part1);
	FREE_MATRIX(Muk);
	FREE_MATRIX(tMuk);
	FREE_MATRIX(Xi);
	FREE_MATRIX(tXi);

	return(-res);

}







double Mstep_Manly_AR_Reg(int p, int T, int n, int q, int K, double *misc_double, double ***Y, double **la, double **gamma, double ***invS, double ***X, double ***beta, double ***invPsi, double *detS, double *detPsi, double *tau){

	int sum_index,k,j,*index,count,i,l1,l2,l,j1,j2,t1,t2,t;    
	double eps, *par, *sum_gamma, Q_value0, *gamma_k, **invPsik, min_value, *Q_value;
	double **inv_part, **part2, **part1, ***MY, **MYi, **tMYi, **Xi, **tXi, **temp1, **temp2, **temp3, *Eig, det, **L, **betak, **Muki, **tMuki, ***Muk, **S, **temp4, **temp5, *Eig2, min, **L2, **invSk, **A1, **A2, **I, Psi2, phi;

	MAKE_VECTOR(par, 2);
	MAKE_VECTOR(sum_gamma, K);
	MAKE_VECTOR(gamma_k, n);
	MAKE_MATRIX(invPsik, T, T);
	MAKE_VECTOR(index, p);
	MAKE_VECTOR(Q_value, K);
	MAKE_MATRIX(inv_part, q, q);
	MAKE_MATRIX(part2, q, p);
	MAKE_MATRIX(part1, q, q);
	MAKE_3ARRAY(MY, p,T,n);
	MAKE_MATRIX(MYi, p, T);
	MAKE_MATRIX(tMYi, T, p);
	MAKE_MATRIX(Xi, T, q);
	MAKE_MATRIX(tXi, q, T);
	MAKE_MATRIX(temp1, q, T);
	MAKE_MATRIX(temp2, q, q);
	MAKE_MATRIX(temp3, q, p);
	MAKE_VECTOR(Eig, q);
	MAKE_MATRIX(L, q, q);
	MAKE_MATRIX(betak, q, p);
 	MAKE_3ARRAY(Muk, p,T,n);
	MAKE_MATRIX(Muki, p, T);
	MAKE_MATRIX(tMuki, T, p);
	MAKE_MATRIX(S, p, p);
	MAKE_MATRIX(temp4, p, T);
	MAKE_MATRIX(temp5, p, p);
	MAKE_VECTOR(Eig2, p);
	MAKE_MATRIX(L2, p, p);
	MAKE_MATRIX(invSk, p, p);
	MAKE_MATRIX(A1, T, T);
	MAKE_MATRIX(A2, T, T);
	MAKE_MATRIX(I, T, T);


	eps = misc_double[0];
	anull(par, 2);

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

			min_value = simplex_Reg(Q_AR_Reg, n, p, T, q, index, Y, gamma_k, invPsik, X, la_nonzero, eps, 0.1);

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

			Q_value[k] = Q_AR_Reg(n, p, T, q, la_nonzero, index, Y, gamma_k, invPsik, X);


			FREE_VECTOR(la_nonzero);
		}

		Q_value0 += Q_value[k];

	}


	

	for(k=0; k<K; k++){


		Anull(inv_part, q, q);
		Anull(part2, q, p);

		cpyv(gamma, k, n, gamma_k);		
		cpyk(invPsi, T, T, k, invPsik);

		Manly_trans_whole_AR(n, p, T, la[k], Y, MY);

		for(i=0; i<n; i++){

			cpyk(MY, p, T, i, MYi);

			tA(MYi, T, p, tMYi);

			cpyk(X, T, q, i, Xi);

			tA(Xi, q, T, tXi);

			multiply(tXi, q, T, invPsik, T, T, temp1);

			multiply(temp1, q, T, Xi, T, q, temp2);

			multiply(temp1, q, T, tMYi, T, p, temp3);
	
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

		Anull(S, p, p);


		for(i=0; i<n; i++){


			cpyk(MY, p, T, i, MYi);

			cpyk(Muk, p, T, i, Muki);

			mat_(p, T, MYi, Muki);

			tA(MYi, T, p, tMYi);

			multiply(MYi, p, T, invPsik, T, T, temp4);

			multiply(temp4, p, T, tMYi, T, p, temp5);

			for(j1=0; j1<p; j1++){			

				for(j2=0; j2<p; j2++){
	
					S[j1][j2] += gamma[i][k] * temp5[j1][j2] / T /sum_gamma[k];

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
	
	
		findPsi2phi_Reg(n, p, T, q, par, MY, X, betak, gamma_k, invSk, eps/10.0);

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



	return Q_value0;



	FREE_VECTOR(par);
	FREE_VECTOR(sum_gamma);
	FREE_VECTOR(gamma_k);
	FREE_MATRIX(invPsik);
	FREE_VECTOR(index);
	FREE_VECTOR(Q_value);
	FREE_MATRIX(inv_part);
	FREE_MATRIX(part2);
	FREE_MATRIX(part1);
	FREE_3ARRAY(MY);
	FREE_MATRIX(MYi);
	FREE_MATRIX(tMYi);
	FREE_MATRIX(Xi);
	FREE_MATRIX(tXi);
	FREE_MATRIX(temp1);
	FREE_MATRIX(temp2);
	FREE_MATRIX(temp3);
	FREE_VECTOR(Eig);
	FREE_MATRIX(L);
	FREE_MATRIX(betak);
 	FREE_3ARRAY(Muk);
	FREE_MATRIX(Muki);
	FREE_MATRIX(tMuki);
	FREE_MATRIX(S);
	FREE_MATRIX(temp4);
	FREE_MATRIX(temp5);
	FREE_VECTOR(Eig2);
	FREE_MATRIX(L2);
	FREE_MATRIX(invSk);
	FREE_MATRIX(A1);
	FREE_MATRIX(A2);
	FREE_MATRIX(I);

}



void findPsi2phi_Reg(int n, int p, int T, int q, double *par, double ***MY, double ***X, double **betak, double *gamma_k, double **invSk, double eps){

	int t1, t2, i, t, j;

	double phi, Psi2, traceA1 = 0.0, traceA2 = 0.0, traceI = 0.0, sum_gamma = 0.0, trA1, trA2, trI, *coeff;
	
	double **A1, **A2, **I, **MYi, **tMYi, **temp1, **temp2, **temp3, **maha1, **maha2, **maha3, *try, *start, **Muk, **tMuk, **Xi;



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
	MAKE_MATRIX(Muk, p, T);
	MAKE_MATRIX(tMuk, T, p);
	MAKE_MATRIX(Xi, T, q);

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

		cpyk(X, T, q, i, Xi);

		multiply(Xi, T, q, betak, q, p, tMuk);

		tA(tMuk, p, T, Muk);

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
	FREE_MATRIX(Muk);
	FREE_MATRIX(tMuk);
	FREE_MATRIX(Xi);
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

void EM_Manly_AR_Reg(int p, int T, int n, int q, int K, double ***Y, double **la, int max_iter, double *misc_double, double *tau, double ***X, double ***beta, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma, int *id, double *ll, int *conv){
	int i,k,iter;
	double eps,loglik_old,loglik,max;

 	eps = misc_double[0];
	loglik_old = -INFINITY;
	iter = 0;

	do{
		loglik = loglik_old; 
		
		iter += 1;


		Estep_Manly_AR_Reg(p, T, n, q, K, Y, la, tau, X, beta, invS, invPsi, detS, detPsi, gamma);

		loglik_old = Mstep_Manly_AR_Reg(p, T, n, q, K, misc_double, Y, la, gamma, invS, X, beta, invPsi, detS, detPsi, tau);

							
	}

	while ((iter < max_iter) && (fabs(loglik - loglik_old) / fabs(loglik_old) > eps));

	ll[0] = mGloglik_Manly_AR_Reg(p, T, n, q, K, Y, la, tau, X, beta, invS, invPsi, detS, detPsi);
	
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

