#ifndef MATMANLYMIX_H
#define MATMANLYMIX_H


void cpyk(double ***a, int nrows, int ncols, int k, double **b);
void cpyk2(double **a, int nrows, int ncols, double ***b, int k);
void cpy(double **a, int nrows, int ncols, double **b);

void array1to2(int a, int b, double *y, double **x);
void array2to1(int a, int b, double *y, double **x);
void array1to3(int a, int b, int c, double *y, double ***x);
void array3to1(int a, int b, int c, double *y, double ***x);



void Manly_trans(int p, int T, double *la, double *nu, double **Y, double **MY);
void Manly_trans_whole(int n, int p, int T, double *la, double *nu, double ***Y, double ***MY);
double mGpdf_Manly_Full(int p, int T, double *la, double *nu, double **Y, double **Mu, double **invS, double **invPsi, double detS, double detPsi);
double mGloglik_Manly_Full(int p, int T, int n, int K, double ***Y, double **la, double **nu, double *tau, double ***Mu, double ***invS, double ***invPsi, double *detS, double *detPsi);
void Estep_Manly_Full(int p, int T, int n, int K, double ***Y, double **la, double **nu, double *tau, double ***Mu, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma);
double Q1(int n, int p, int T, double *la_nonzero, double *nu, int *index, double ***Y, double *gamma_k, double **invPsik, int Mu_type);
double Q2(int n, int p, int T, double *nu_nonzero, double *la, int *index, double ***Y, double *gamma_k, double **invPsik, int Mu_type);
double Mstep_Manly_Full(int p, int T, int n, int K, double *misc_double, double ***Y, double **la, double **nu, double **gamma, double ***invS, double ***Mu, double ***invPsi, double *detS, double *detPsi, double *tau, int Mu_type);
void EM_Manly_Full(int p, int T, int n, int K, double ***Y, double **la, double **nu, int max_iter, double *misc_double, double *tau, double ***Mu, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma, int *id, double *ll, int *conv, int Mu_type);


void EigValDec(int size, double *W, double **A, double (*determinant));
void Anull3(double ***X, int ax, int bx, int cx);
void Anull(double **X, int ax, int bx);
void anull(double *x, int p);
void anulli(int *x, int p);
void XAXt(double **X, int p, double **A, double **Res);
void cpy1(double ***a, int k, int nrows, int ncols, double **b);
void dsyev_(char *JOBZp, char *UPLOp,int *Np, double *A, int *LDAp, double *Wp, double *WORK, int *LWORK, int *INFOp);
int mat_(int a, int b,double **Res, double **Y);
void tA(double **A, int a, int b, double **Res);
void multiply(double **a, int arows, int acols, double **b, int brows, int bcols, double **c);
void cpyv(double **A, int col, int nrows, double *V);
void matxvec(double **a, int arows, int acols, double *x, int xrows, double *y);



double simplex1(double (*func)(int, int, int, double *, double *, int *, double ***, double *, double **, int), int n1, int p, int T, double *nu, int *index, double ***X, double *gamma_k, double **invPsik, double *start, double EPSILON, double scale, int Mu_type);

double simplex2(double (*func)(int, int, int, double *, double *, int *, double ***, double *, double **, int), int n1, int p, int T, double *la, int *index, double ***X, double *gamma_k, double **invPsik, double *start, double EPSILON, double scale, int Mu_type);


double mGpdf_Manly_Reg(int p, int T, int q, double *la, double *nu, double **Y, double **X, double **beta, double **invS, double **invPsi, double detS, double detPsi);
double mGloglik_Manly_Reg(int p, int T, int n, int q, int K, double ***Y, double **la, double **nu, double *tau, double ***X, double ***beta, double ***invS, double ***invPsi, double *detS, double *detPsi);
void Estep_Manly_Reg(int p, int T, int n, int q, int K, double ***Y, double **la, double **nu, double *tau, double ***X, double ***beta, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma);
double Q_Reg1(int n, int p, int T, int q, double *la_nonzero, double *nu, int *index, double ***Y, double *gamma_k, double **invPsik, double ***X);
double Q_Reg2(int n, int p, int T, int q, double *nu_nonzero, double *la, int *index, double ***Y, double *gamma_k, double **invPsik, double ***X);
double Mstep_Manly_Reg(int p, int T, int n, int q, int K, double *misc_double, double ***Y, double **la, double **nu, double **gamma, double ***invS, double ***X, double ***beta, double ***invPsi, double *detS, double *detPsi, double *tau);
void EM_Manly_Reg(int p, int T, int n, int q, int K, double ***Y, double **la, double **nu, int max_iter, double *misc_double, double *tau, double ***X, double ***beta, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma, int *id, double *ll, int *conv);


double simplex_Reg1(double (*func)(int, int, int, int, double *, double *, int *, double ***, double *, double **, double ***), int n1, int p, int T, int q, double *nu, int *index, double ***X, double *gamma_k, double **invPsik, double ***X1, double *start, double EPSILON, double scale);
double simplex_Reg2(double (*func)(int, int, int, int, double *, double *, int *, double ***, double *, double **, double ***), int n1, int p, int T, int q, double *la, int *index, double ***X, double *gamma_k, double **invPsik, double ***X1, double *start, double EPSILON, double scale);



void Manly_trans_AR(int p, int T, double *la, double **Y, double **MY);
void Manly_trans_whole_AR(int n, int p, int T, double *la, double ***Y, double ***MY);
double mGpdf_Manly_AR(int p, int T, double *la, double **Y, double **Mu, double **invS, double **invPsi, double detS, double detPsi);
double mGloglik_Manly_AR(int p, int T, int n, int K, double ***Y, double **la, double *tau, double ***Mu, double ***invS, double ***invPsi, double *detS, double *detPsi);
void Estep_Manly_AR(int p, int T, int n, int K, double ***Y, double **la, double *tau, double ***Mu, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma);
double Q_AR(int n, int p, int T, double *la_nonzero, int *index, double ***Y, double *gamma_k, double **invPsik, int Mu_type);
double Mstep_Manly_AR(int p, int T, int n, int K, double *misc_double, double ***Y, double **la, double **gamma, double ***invS, double ***Mu, double ***invPsi, double *detS, double *detPsi, double *tau, int Mu_type);
void findPsi2phi(int n, int p, int T, double *par, double ***MY, double **Muk, double *gamma_k, double **invSk, double eps);
double eq3(double x, double *coeff);
void rootfinding(double (*func)(double, double *), double *start, double *coeff, double eps);
void EM_Manly_AR(int p, int T, int n, int K, double ***Y, double **la, int max_iter, double *misc_double, double *tau, double ***Mu, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma, int *id, double *ll, int *conv, int Mu_type);


double simplex_AR(double (*func)(int, int, int, double *, int *, double ***, double *, double **, int), int n1, int p, int T, int *index, double ***X, double *gamma_k, double **invPsik, double *start, double EPSILON, double scale, int Mu_type);


double mGpdf_Manly_AR_Reg(int p, int T, int q, double *la, double **Y, double **X, double **beta, double **invS, double **invPsi, double detS, double detPsi);
double mGloglik_Manly_AR_Reg(int p, int T, int n, int q, int K, double ***Y, double **la, double *tau, double ***X, double ***beta, double ***invS, double ***invPsi, double *detS, double *detPsi);
void Estep_Manly_AR_Reg(int p, int T, int n, int q, int K, double ***Y, double **la, double *tau, double ***X, double ***beta, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma);
double Q_AR_Reg(int n, int p, int T, int q, double *la_nonzero, int *index, double ***Y, double *gamma_k, double **invPsik, double ***X);
double Mstep_Manly_AR_Reg(int p, int T, int n, int q, int K, double *misc_double, double ***Y, double **la, double **gamma, double ***invS, double ***X, double ***beta, double ***invPsi, double *detS, double *detPsi, double *tau);
void findPsi2phi_Reg(int n, int p, int T, int q, double *par, double ***MY, double ***X, double **betak, double *gamma_k, double **invSk, double eps);
void EM_Manly_AR_Reg(int p, int T, int n, int q, int K, double ***Y, double **la, int max_iter, double *misc_double, double *tau, double ***X, double ***beta, double ***invS, double ***invPsi, double *detS, double *detPsi, double **gamma, int *id, double *ll, int *conv);

double simplex_Reg(double (*func)(int, int, int, int, double *, int *, double ***, double *, double **, double ***), int n1, int p, int T, int q, int *index, double ***X, double *gamma_k, double **invPsik, double ***X1, double *start, double EPSILON, double scale);

/* WCC: "libEVD.c" and "libEVD_LAPACK.c" */
#ifndef __HAVE_R_
	void cephes_symmeigens_down(int p, double *eval, double **A, double (*determinant));
#else
	void dsyev_(char *JOBZp, char *UPLOp,int *Np, double *A, int *LDAp, double *Wp, double *WORK, int *LWORK, int *INFOp);
	void EigValDec(int size, double *W, double **A, double (*determinant));
#endif

#endif /* MATMANLYMIX_H */
