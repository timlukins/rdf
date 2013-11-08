
void
calculate_curvature (GPtrArray* pp_points,double pv_nx,double pv_ny,double pv_nz,double* pp_E,double* pp_F,double* pp_G,double* pp_L,double* pp_M,double* pp_N,double* pp_err, double* pp_a, double* pp_b, double* pp_c)
{
  Matrix* R = matrix_create(3,3);
  Matrix* N = matrix_create(3,3);
  Matrix* r1;
  Matrix* r2;
  Matrix* r3;
  Matrix* I = matrix_create(3,3);
  Matrix* i = matrix_create(3,1);
  Matrix* n = matrix_create(3,1);
  Matrix* v = matrix_create(3,1);
  Matrix* p = matrix_create(3,1);
  Matrix* M = matrix_create(pp_points->len-1,3); // Minus-1 for first
  Matrix* B = matrix_create(pp_points->len-1,1);
  Matrix* t; 
  Matrix* q;
  Matrix* r;
  Matrix* s;
  Matrix* Xu = matrix_create(3,1);
  Matrix* Xv = matrix_create(3,1);
  Matrix* Xuu = matrix_create(3,1);
  Matrix* Xuv = matrix_create(3,1);
  Matrix* Xvv = matrix_create(3,1);
  long j,k;
  double x,y,z;
  double a,b,c;
  GPtrArray* lp_aligned_points;
  GtsPoint* lp_p;
  
  /* Get centre vertex point in question - the first one! */
    
  GtsPoint* lp_v = GTS_POINT(g_ptr_array_index(pp_points,0));
  vector_set(n,pv_nx,pv_ny,pv_nz); 

  /* Build rotation matrix to align to unit normal */
 
  matrix_setidentity(I);
  vector_set(i,0,0,1.0); // Z axis unit straight up  

  vector_normalise(n);
  /* Test if already pointing up/down */
  
  if (matrix_get(n,3,1)==1.0)
  {
    r1 = matrix_create(3,1); // NOTE using LEFT HANDED co-ordinate system - ie. r1==x (not r1=y)
    vector_set(r1,-1.0,0,0);
  } 
  else if (matrix_get(n,3,1)==-1.0) 
  {
    r1 = matrix_create(3,1);
    vector_set(r1,1.0,0,0);
  } 
  else // Work out z rotation 
  {
    N = matrix_multiplyRT(n,n);
    matrix_sub(I,N);
    r1 = matrix_multiply(I,i);
    vector_normalise(r1); 
  } 
  r3 = matrix_copy(n);  
  r2 = vector_cross(n,r1);
  
  matrix_setidentity(R);

/*
  matrix_set(R,1,1,matrix_get(r1,1,1));
  matrix_set(R,1,2,matrix_get(r1,2,1));
  matrix_set(R,1,3,matrix_get(r1,3,1));
  matrix_set(R,2,1,matrix_get(r2,1,1));
  matrix_set(R,2,2,matrix_get(r2,2,1));
  matrix_set(R,2,3,matrix_get(r2,3,1));
  matrix_set(R,3,1,matrix_get(r3,1,1)); 
  matrix_set(R,3,2,matrix_get(r3,2,1));
  matrix_set(R,3,3,matrix_get(r3,3,1));
*/
  //matrix_print(n); // USE TO TEST if R correct rotates normal to z axis!
  //matrix_print(R);

    /* test initial quadric - ACTUALLY, above solution should work!*/

    /*
    lm_R = [pm_P(:,1).^2 pm_P(:,1).*pm_P(:,2) pm_P(:,2).^2 pm_P(:,1) pm_P(:,2)];
    lv_S = pm_P(:,3);
    %Q = inv(lm_R'*lm_R)*lm_R'*lv_S;
    %Q = pinv(lm_R)*lv_S;
    Q = lscov(lm_R,lv_S);
    n = [-Q(4) -Q(5) 1]/(1 + Q(4)^2 + Q(5)^2);
    */
  //%} while (FALSE); // or while n(3)!=1.0 

  /* Go through the points - miss first as it is repeated */

  matrix_setzero(M);
  matrix_setzero(B);

  lp_aligned_points = g_ptr_array_new(); 
  g_ptr_array_add(lp_aligned_points,gts_point_new(gts_point_class(),0,0,0)); // Add centre point at origin
  for (j=1;j<pp_points->len;j++) /* NOTE don't include first point - the origin */ 
  {
    lp_p = GTS_POINT(g_ptr_array_index(pp_points,j));

//printf("%f %f %f\n",lp_p->x,lp_p->y,lp_p->z);

    /* Translate and rotate point relative to centre */   
  
    matrix_set(p,1,1,lp_p->x - lp_v->x);
    matrix_set(p,2,1,lp_p->y - lp_v->y);
    matrix_set(p,3,1,lp_p->z - lp_v->z);

    t = matrix_multiply(R,p);

    /* Push this translated/rotated point into our list for minimisation */ 
    g_ptr_array_add(lp_aligned_points,gts_point_new(gts_point_class(),
      matrix_get(t,1,1), 
      matrix_get(t,2,1), 
      matrix_get(t,3,1)));

    /* Build Least Squares matrices */

    x = matrix_get(t,1,1);
    y = matrix_get(t,2,1);
    z = matrix_get(t,3,1);

  //if (lp_v->x>0.95 && lp_v->y>0.95 && lp_v->x<1.01 && lp_v->y<1.01)  //central patch only
  /* if (lp_v->x>0.80 && lp_v->y>0.80 && lp_v->x<1.20 && lp_v->y<1.20) 
   {
      //printf("%f %f %f\n",lp_p->x,lp_p->y,lp_p->z);
      //matrix_print(R);  
     //printf("%f %f %f\n",x,y,z);
   }*/
    
    matrix_set(B,j,1,z);
    
    matrix_set(M,j,1, x * x);
    matrix_set(M,j,2, x * y);
    matrix_set(M,j,3, y * y);
  }

  /* Perform Least Squares estimation */
  //gw_x = gwave_util::get_inverse((gw_m.Transpose() * gw_m)) * gw_m.Transpose() * gw_b; 
  //r = matrix_multiplyT(M,M);
  //matrix_inverse(r);
  //s = matrix_multiplyT(M,B);
  //q = matrix_multiply(r,s);
  // BIG TODO: add weighting matrix based on distance from centre!
  //q = matrix_create(3,1); 
  q = matrix_solvels(M,B);

  /* Extract quadric co-effs */

  a = matrix_get(q,1,1);
  b = matrix_get(q,2,1);
  c = matrix_get(q,3,1);
 
  *pp_a = a; // Return the quadric
  *pp_b = b;
  *pp_c = c;

  /* NOTE: this is only our INITIAL estimate */
  /* NOW use non-linear estimation for minfit to make even better! */

  double error=0.0;
  for (j=0;j<lp_aligned_points->len;j++)  
  {
    lp_p = GTS_POINT(g_ptr_array_index(lp_aligned_points,j));
    error = ( a*(lp_p->x*lp_p->x) + b*(lp_p->x*lp_p->y) + c*(lp_p->y*lp_p->y) ) ;
    error = (lp_p->z-error)*(lp_p->z-error);
    *pp_err += error;
  }
  *pp_err /= lp_aligned_points->len; // Mean Square Error

  //printf("%d err = %f\n",j,*pp_err);
  if (*pp_err>0.0) // TODO: Tolerance here? 
  {
    *pp_err = matrix_solveMINquadric(&a,&b,&c,lp_aligned_points); // TODO: this could be a switch 
    //printf("Refining fit! Error = %f\n",*pp_err);
  }
  // First and second order partial derivatives of
  //            2              2
  // z(x,y)= a x  + b x y + c y 

  x = 1.0;//nx(2,2)-nx(1,1);  // X & Y need to the be the "new length" of x, y sample points
  y = 1.0;//ny(2,2)-ny(1,1);

  /* Calculate first fund form */

  vector_set(Xu,x, 0, 2*a*x + b*y);
  vector_set(Xv,0, y, b*x + 2*c*y);
  *pp_E           =   vector_dot(Xu,Xu);
  *pp_F           =   vector_dot(Xu,Xv);
  *pp_G           =   vector_dot(Xv,Xv);

  /* Cacluate second fund form */

  vector_set(Xuu,0, 0, 2*a);
  vector_set(Xuv,0, 0, b);
  vector_set(Xvv,0, 0, 2*c);
  *pp_L           =   vector_dot(Xuu,i); // i is straight up Z axis - which is what this quadric is
  *pp_M           =   vector_dot(Xuv,i);
  *pp_N           =   vector_dot(Xvv,i);

  /* Tidy up */

  matrix_delete(q);
  //matrix_delete(r);
  //matrix_delete(s);
  matrix_delete(t);
  matrix_delete(B);
  matrix_delete(M);
  matrix_delete(p);
  matrix_delete(v);
  matrix_delete(n);
  matrix_delete(i);
  matrix_delete(I);
  matrix_delete(r1);
  matrix_delete(r2);
  matrix_delete(r3);
  matrix_delete(R);
  matrix_delete(N);
  matrix_delete(Xu);
  matrix_delete(Xv);
  matrix_delete(Xuu);
  matrix_delete(Xuv);
  matrix_delete(Xvv);
  
  for (k=0;k<lp_aligned_points->len;k++) 
    gts_object_destroy(GTS_OBJECT(g_ptr_array_index(lp_aligned_points,k)));  
  g_ptr_array_free(lp_aligned_points,TRUE); 

}
