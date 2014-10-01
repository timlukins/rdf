/***************************************************

fitquadric - Fit a quadric to points

....................................................

Requires:
  - gts 
  - gsl
(and of course...) 
  - Matlab   (www.mathworks.com)

....................................................

Version 1.0 - Tim Lukins (July/2006).

****************************************************/

/* Includes */

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multimin.h>
#include <mex.h>

/**
 *
 * Objective function for Euclidean error between points and quadric evaluation.
 *
 */
/*
double
euclidean_fit(const gsl_vector* pp_quadric, void* pp_points)
{
  double rv_mse = 0.0;
  
  Quadric lv_quad;
  lv_quad.a = gsl_vector_get(pp_quadric,0);
  lv_quad.b = gsl_vector_get(pp_quadric,1);
  lv_quad.c = gsl_vector_get(pp_quadric,2);
  lv_quad.h = gsl_vector_get(pp_quadric,3);
  lv_quad.g = gsl_vector_get(pp_quadric,4);
  lv_quad.f = gsl_vector_get(pp_quadric,5);
  lv_quad.u = gsl_vector_get(pp_quadric,6);
  lv_quad.v = gsl_vector_get(pp_quadric,7);
  lv_quad.w = gsl_vector_get(pp_quadric,8);
  lv_quad.d = gsl_vector_get(pp_quadric,9);

  GtsVertexList* lp_points = (GtsVertexList*)pp_points;
  
  for (GtsVertexList::const_iterator p=lp_points->begin(); p!=lp_points->end(); p++)
  {
     GtsVertex* lp_pt = (*p);

     rv_mse += gsl_pow_2(evaluate_quadric(lv_quad,lp_pt->p.x,lp_pt->p.y,lp_pt->p.z));
  }

  rv_mse /= lp_points->size();

  return rv_mse;
}
*/
/* MATLAB entry point */
/* NOTE: actual function name resolved by name of the file! */

void 
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray* lp_pntdata;
  double*  lp_pnts;
  double*  lp_quadric;
  long lv_npts, p;
  gsl_matrix* lp_H      = gsl_matrix_calloc(10,10);
  gsl_matrix* lp_DH     = gsl_matrix_calloc(10,10);
  gsl_vector* lp_evals  = gsl_vector_calloc(10);
  gsl_matrix* lp_evect  = gsl_matrix_calloc(10,10);
  gsl_vector* lp_step   = gsl_vector_calloc(10);
  gsl_eigen_symmv_workspace* lp_espace = gsl_eigen_symmv_alloc(10);
  gsl_multimin_fminimizer* lp_mini     = gsl_multimin_fminimizer_alloc (gsl_multimin_fminimizer_nmsimplex, 10);
  double x,y,z;
  unsigned long i = 0;
  int status      = 0;
  double error    = 0.0;

  /* Test that we have 1 input and 1 output. */

  if (nlhs!=1 && nrhs!=1 )
  {
    printf("quadricfit: Error in number of arguments/assignment.\n");
    return;
  }
  
  /* Initalise return structure */

  plhs[0] = mxCreateDoubleMatrix(10,1, mxREAL); 
  lp_quadric = mxGetPr(plhs[0]);
  
  /* Get the points themselves */
  
  lp_pntdata = prhs[0];
  if (mxGetN(lp_pntdata)!=3)
  {
    printf("quadricfit: Array must be 3xM.\n");
    return;
  }
  lp_pnts = mxGetPr(lp_pntdata);
  lv_npts = mxGetM(lp_pntdata);
 
  /* Perform initial Taubin fitting... */
   
  gsl_matrix_set_all(lp_H,0.0);
  gsl_matrix_set_all(lp_DH,0.0);
  
  for (p=0; p<lv_npts; p++)
  {
    x = lp_pnts[p]; 
    y = lp_pnts[p+lv_npts]; 
    z = lp_pnts[p+2*lv_npts]; 

    double hv[] = {gsl_pow_2(x),
                   gsl_pow_2(y),
                   gsl_pow_2(z),
                   2*x*y, 
                   2*x*z, 
                   2*z*y,
                   2*x,
                   2*y,
                   2*z,
                   1};

    gsl_matrix_view h = gsl_matrix_view_array(hv,10,1);

    gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, &h.matrix, &h.matrix, 1.0, lp_H); 

    double jv[] = {2*x, 0,  0,
                    0, 2*y, 0,
                    0,  0, 2*z,
                   2*y,2*x, 0,
                   2*z, 0, 2*x,
                    0, 2*z,2*y,
                    2,  0,  0,
                    0,  2,  0,
                    0,  0,  2,
                    0,  0,  0};

    gsl_matrix_view j = gsl_matrix_view_array(jv,10,3); /* Jacobian built! */

    gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1.0, &j.matrix, &j.matrix, 1.0, lp_DH);

  }

  /* For generalised eignenvector solution A*v=lambda*B*v -> inv(B)*A*v=lambda*v */
  
  gsl_blas_dtrsm(CblasLeft,CblasUpper,CblasNoTrans,CblasUnit,1.0,lp_DH,lp_H);
  
  gsl_eigen_symmv(lp_H, lp_evals, lp_evect, lp_espace);
  gsl_eigen_symmv_sort (lp_evals, lp_evect, GSL_EIGEN_SORT_VAL_DESC);

  gsl_vector_view lp_best = gsl_matrix_column(lp_evect,0); /* Smallest eigenvalue solution */

  /* Then proceed with some Euclidean minimization... */
/*
  gsl_multimin_function lv_objective;
  lv_objective.f      = &euclidean_fit;
  lv_objective.n      = 10;
  lv_objective.params = (void*)&lp_pts;
  
  gsl_vector_set_all (lp_step, 1.0); // Set inital step size for search

  gsl_multimin_fminimizer_set (lp_mini, &lv_objective, &lp_best.vector, lp_step);

  do
  {
    i++;
    status = gsl_multimin_fminimizer_iterate(lp_mini);

    if (status)
      break;

    error = gsl_multimin_fminimizer_size (lp_mini);
    status = gsl_multimin_test_size (error, *pp_error);
  }
  while (status == GSL_CONTINUE && i < 10000); // Maximum iterations!
 */ 
  /* Update error to final error value... */

  /* *pp_error = error;  TODO return error as well */

  /* Set best solution... */

  /*lp_best = gsl_vector_subvector(lp_mini->x,0,10);*/ /* Acutally, all of vector*/
 
  /* Fill in return result...*/

  lp_quadric[0] = gsl_vector_get(&lp_best.vector,0);
  lp_quadric[1] = gsl_vector_get(&lp_best.vector,1);
  lp_quadric[2] = gsl_vector_get(&lp_best.vector,2);
  lp_quadric[3] = gsl_vector_get(&lp_best.vector,3);
  lp_quadric[4] = gsl_vector_get(&lp_best.vector,4);
  lp_quadric[5] = gsl_vector_get(&lp_best.vector,5);
  lp_quadric[6] = gsl_vector_get(&lp_best.vector,6);
  lp_quadric[7] = gsl_vector_get(&lp_best.vector,7);
  lp_quadric[8] = gsl_vector_get(&lp_best.vector,8);
  lp_quadric[9] = gsl_vector_get(&lp_best.vector,9);

  /* Tidy up all this mess... */

  gsl_matrix_free(lp_H);
  gsl_matrix_free(lp_DH);
  gsl_matrix_free(lp_evect);
  gsl_vector_free(lp_evals);
  gsl_vector_free(lp_step);
  gsl_eigen_symmv_free(lp_espace);
  gsl_multimin_fminimizer_free(lp_mini);

  /* Return... */

  return;
}

