/* 
 * Program: rosen.c
 * Author : Michael F. Hutt
 * http://www.mikehutt.com
 * 11/3/97
 * $Id: crosen.c,v 1.4 2007/07/10 12:42:02 mike Exp $
 *
 * Copyright (c) 1997-2004 <Michael F. Hutt>
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *
 * An implementation of the Nelder-Mead simplex method.
 */

#include "array.h"
#include <math.h>

#define MAX_IT      1000      /* maximum number of iterations */
#define ALPHA       1.0       /* reflection coefficient */
#define BETA        0.5       /* contraction coefficient */
#define GAMMA       2.0       /* expansion coefficient */


double simplex_Reg2(double (*func)(int, int, int, int, double *, double *, int *, double ***, double *, double **, double ***), int n1, int p, int T, int q, double *la, int *index, double ***X, double *gamma_k, double **invPsik, double ***X1, double *start, double EPSILON, double scale)
{
  
  int vs;        	/* vertex with smallest value */
  int vh;        	/* vertex with next smallest value */
  int vg;        	/* vertex with largest value */
  
  int i,j,m,row, n;
  int k;		/* track the number of function evaluations */
  int itr;		/* track the number of iterations */
  
  double **v;           /* holds vertices of simplex */
  double pn,qn;         /* values used to create initial simplex */
  double *f;            /* value of function at each vertex */
  double fr;            /* value of function at reflection point */
  double fe;            /* value of function at expansion point */
  double fc;            /* value of function at contraction point */
  double *vr;           /* reflection - coordinates */
  double *ve;           /* expansion - coordinates */
  double *vc;           /* contraction - coordinates */
  double *vm;           /* centroid - coordinates */
  double min;
  
  double fsum,favg,s,cent;
  n = 0;
  for(j=0;j<T;j++){
    if(index[j]==1){
	n +=1;
    }
  }
  
  /* dynamically allocate arrays */
  
  /* allocate the rows of the arrays */

  /* changed from original */

  MAKE_MATRIX(v, n+1, n);
  MAKE_VECTOR(f, n+1);
  MAKE_VECTOR(vr, n);
  MAKE_VECTOR(ve, n);
  MAKE_VECTOR(vc, n);
  MAKE_VECTOR(vm, n);
  
  /* original memory allocation
  v =  (double **) malloc ((n+1) * sizeof(double *));
  f =  (double *) malloc ((n+1) * sizeof(double));
  vr = (double *) malloc (n * sizeof(double));
  ve = (double *) malloc (n * sizeof(double));  
  vc = (double *) malloc (n * sizeof(double));  
  vm = (double *) malloc (n * sizeof(double));  
  
  // allocate the columns of the arrays 
  for (i=0;i<=n;i++) {
    v[i] = (double *) malloc (n * sizeof(double));
  }
  */
  
  /* create the initial simplex */
  /* assume one of the vertices is 0,0 */
  
  pn = scale*(sqrt(n+1)-1+n)/(n*sqrt(2));
  qn = scale*(sqrt(n+1)-1)/(n*sqrt(2));
  
  for (i=0;i<n;i++) {
    v[0][i] = start[i];
  }
  
  for (i=1;i<=n;i++) {
    for (j=0;j<n;j++) {
      if (i-1 == j) {
	v[i][j] = pn + start[j];
      }
      else {
	v[i][j] = qn + start[j];
      }
    }
  }
  
  /* find the initial function values */
  for (j=0;j<=n;j++) {
    f[j] = func(n1, p, T, q, v[j], la, index, X, gamma_k, invPsik, X1);
  }
  
  k = n+1;
  
  
  /* begin the main loop of the minimization */
  for (itr=1;itr<=MAX_IT;itr++) {     
    /* find the index of the largest value */
    vg=0;
    for (j=0;j<=n;j++) {
      if (f[j] > f[vg]) {
	vg = j;
      }
    }
    
    /* find the index of the smallest value */
    vs=0;
    for (j=0;j<=n;j++) {
      if (f[j] < f[vs]) {
	vs = j;
      }
    }
    
    /* find the index of the second largest value */
    vh=vs;
    for (j=0;j<=n;j++) {
      if (f[j] > f[vh] && f[j] < f[vg]) {
	vh = j;
      }
    }
    
    /* calculate the centroid */
    for (j=0;j<=n-1;j++) {
      cent=0.0;
      for (m=0;m<=n;m++) {
	if (m!=vg) {
	  cent += v[m][j];
	}
      }
      vm[j] = cent/n;
    }
    
    /* reflect vg to new vertex vr */
    for (j=0;j<=n-1;j++) {
      vr[j] = (1+ALPHA)*vm[j] - ALPHA*v[vg][j];
    }
    fr = func(n1, p, T, q, vr, la, index, X, gamma_k, invPsik, X1);
    k++;
    
    /* added <= */
    if (fr <= f[vh] && fr > f[vs]) {
      for (j=0;j<=n-1;j++) {
	v[vg][j] = vr[j];
      }
      f[vg] = fr;
    }
    
    /* investigate a step further in this direction */
    /* added <= */
    if ( fr <=  f[vs]) {
      for (j=0;j<=n-1;j++) {
	ve[j] = GAMMA*vr[j] + (1-GAMMA)*vm[j];
      }
      fe = func(n1, p, T, q, ve, la, index, X, gamma_k, invPsik, X1);

      k++;
      
   
      if (fe < fr) {
	for (j=0;j<=n-1;j++) {
	  v[vg][j] = ve[j];
	}
	f[vg] = fe;
      }
      else {
	for (j=0;j<=n-1;j++) {
	  v[vg][j] = vr[j];
	}
	f[vg] = fr;
      }
    }
    
    /* check to see if a contraction is necessary */
    if (fr > f[vh]) {
      for (j=0;j<=n-1;j++) {
	vc[j] = BETA*v[vg][j] + (1-BETA)*vm[j];
      }
      fc = func(n1, p, T, q, vc, la, index, X, gamma_k, invPsik, X1);
      k++;
      if (fc < f[vg]) {
	for (j=0;j<=n-1;j++) {
	  v[vg][j] = vc[j];
	}
	f[vg] = fc;
      }
      /* at this point the contraction is not successful,
	 we must halve the distance from vs to all the 
	 vertices of the simplex and then continue.
	 10/31/97 - modified to account for ALL vertices. 
      */
      else {
	for (row=0;row<=n;row++) {
	  if (row != vs) {
	    for (j=0;j<=n-1;j++) {
	      v[row][j] = v[vs][j]+(v[row][j]-v[vs][j])/2.0;
	    }
	  }
	}
	f[vg] = func(n1, p, T, q, v[vg], la, index, X, gamma_k, invPsik, X1);

	k++;
	f[vh] = func(n1, p, T, q, v[vh], la, index, X, gamma_k, invPsik, X1);
	k++;
	
	
      }
    }
    
    
    /* test for convergence */
    fsum = 0.0;
    for (j=0;j<=n;j++) {
      fsum += f[j];
    }
    favg = fsum/(n+1);
    s = 0.0;
    for (j=0;j<=n;j++) {
      s += pow((f[j]-favg),2.0)/(n);
    }
    s = sqrt(s);
    if (s < EPSILON) break;
  }
  /* end main loop of the minimization */
  
  /* find the index of the smallest value */
  vs=0;
  for (j=0;j<=n;j++) {
    if (f[j] < f[vs]) {
      vs = j;
    }
  }


  /* printf("The minimum was found at\n");   */
  for (j=0;j<n;j++) {
    start[j] = v[vs][j];
  }
  min=func(n1, p, T, q, v[vs], la, index, X, gamma_k, invPsik, X1);
  k++;
  /*	printf("%d Function Evaluations\n",k);
	printf("%d Iterations through program\n",itr);*/



  /* changed from original */

  FREE_MATRIX(v);
  FREE_VECTOR(f);
  FREE_VECTOR(vr);
  FREE_VECTOR(ve);
  FREE_VECTOR(vc);
  FREE_VECTOR(vm);
  
  /* original memory free
  free(f);
  free(vr);
  free(ve);
  free(vc);
  free(vm);
  for (i=0;i<=n;i++) {
    free (v[i]);
  }  
  free(v);
  */


//printf("Final simplex value la %lf %lf\n",start[0], start[1]);

  return min;
}
