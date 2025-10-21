void evaluate_variables_fTTT(double v, double cosp, double y){
  double sinp = sqrt(1-pow(cosp,2));
  pik = pi$plus*v+pi$minus/(2*v)-pTi*cosp;
  pjk = pj$plus*v+pj$minus/(2*v)-pTj*(cosxi*cosp-y*sinxi*sinp);
  p1$minus = (p1[0]-p1[3])/sqrt(2);
  p2$minus = (p2[0]-p2[3])/sqrt(2);
  p1$plus  = (p1[0]+p1[3])/sqrt(2);
  p2$plus  = (p2[0]+p2[3])/sqrt(2);
  p1k = p1$plus*v + p1$minus/(2*v);
  p2k = p2$plus*v + p2$minus/(2*v);
  aij = mi2*pjk/2/pipj/pik;
  aji = mj2*pik/2/pipj/pjk;
  ai1 = mi2*p1k/2/p1pi/pik;
  ai2 = mi2*p2k/2/p2pi/pik;
  aj1 = mj2*p1k/2/p1pj/pjk;
  aj2 = mj2*p2k/2/p2pj/pjk;
  v34 = sqrt(1-4*aij*aji);
  vp = (1+v34)/2;
  vm = (1-v34)/2;
  x = vm/vp;
  dij = 1-2*aij;
  dji = 1-2*aji;
  di1 = 1-2*ai1;
  di2 = 1-2*ai2;
  dj1 = 1-2*aj1;
  dj2 = 1-2*aj2;
}

double fTTT_m2(double v, double cosp, double y){
  double res = 0;
  evaluate_variables_fTTT(v,cosp,y);
  res = ((p1pj*p2k*pik - p1k*p2pj*pik - p1pi*p2k*pjk + p1k*p2pi*pjk)*(1 + v34)*pow(p1k,-1)*pow(p2k,-1)*pow(pik,-1)*pow(pjk,-1)*pow(v,-1)*pow(v34,-1)*pow(1 - pow(cosp,2),-0.5))/4.;
  return res;
}

double fTTT_m1(double v, double cosp, double y){
  double res = 0;
  evaluate_variables_fTTT(v,cosp,y);
  res = (pow(p1k,-1)*pow(p2k,-1)*pow(pik,-1)*pow(pjk,-1)*pow(v,-1)*((p1pj*p2k*pik - p1k*p2pj*pik - p1pi*p2k*pjk + p1k*p2pi*pjk)*
        log((p1p2*pow(p1k,-1)*pow(p2k,-1))/2.) + (-(p1pj*p2k*pik) + p1k*p2k*pipj - p1k*p2pi*pjk + p1p2*pik*pjk)*
        ((1 + di1)*log(2) - 2*ai1*log(ai1) - (1 + di1)*log(p1pi*pow(p1k,-1)*pow(pik,-1)))*pow(1 + di1,-1) + 
       (p1k*p2pj*pik - p1k*p2k*pipj + p1pi*p2k*pjk - p1p2*pik*pjk)*((1 + di2)*log(2) - 2*ai2*log(ai2) - (1 + di2)*log(p2pi*pow(p2k,-1)*pow(pik,-1)))*
        pow(1 + di2,-1) + (p1k*p2pj*pik - p1k*p2k*pipj + p1pi*p2k*pjk - p1p2*pik*pjk)*
        ((1 + dj1)*log(2) - 2*aj1*log(aj1) - (1 + dj1)*log(p1pj*pow(p1k,-1)*pow(pjk,-1)))*pow(1 + dj1,-1) + 
       (-(p1pj*p2k*pik) + p1k*p2k*pipj - p1k*p2pi*pjk + p1p2*pik*pjk)*((1 + dj2)*log(2) - 2*aj2*log(aj2) - (1 + dj2)*log(p2pj*pow(p2k,-1)*pow(pjk,-1)))*
        pow(1 + dj2,-1) + (p1pj*p2k*pik - p1k*p2pj*pik - p1pi*p2k*pjk + p1k*p2pi*pjk)*(1 + v34)*(-2*log(2) - log(1 - pow(cosp,2)))*pow(v34,-1) + 
       (-(p1pj*p2k*pik) + p1k*p2pj*pik + p1pi*p2k*pjk - p1k*p2pi*pjk)*
        (-2*log(v34) + (1 + 2*v34)*log(vp) + (-1 + 2*v34)*log((pipj*pow(pik,-1)*pow(pjk,-1))/2.) + 
          2*((aij - vm)*log(aij*pow(vp,-1)) + (aji - vm)*log(aji*pow(vp,-1)))*pow(dij + dji,-1))*pow(v34,-1))*pow(1 - pow(cosp,2),-0.5))/4.;
  return res;
}

double fTTT_0(double v, double cosp, double y){
  double res = 0;
  evaluate_variables_fTTT(v,cosp,y);
  res = (pow(1 + di1,-1)*pow(1 + di2,-1)*pow(1 + dj1,-1)*pow(1 + dj2,-1)*pow(dij + dji,-1)*pow(p1k,-1)*pow(p2k,-1)*pow(pik,-1)*pow(pjk,-1)*pow(v,-1)*pow(v34,-1)*
     (-12*(-2*log(2) - log(1 - pow(cosp,2)))*((1 + di1)*(1 + di2)*(1 + dj1)*(1 + dj2)*(dij + dji)*
           (-(p1pj*p2k*pik) + p1k*p2pj*pik + p1pi*p2k*pjk - p1k*p2pi*pjk)*v34*log((p1p2*pow(p1k,-1)*pow(p2k,-1))/2.) - 
          (1 + di2)*(1 + dj1)*(1 + dj2)*(dij + dji)*(-(p1pj*p2k*pik) + p1k*p2k*pipj - p1k*p2pi*pjk + p1p2*pik*pjk)*v34*
           ((1 + di1)*log(2) - 2*ai1*log(ai1) - (1 + di1)*log(p1pi*pow(p1k,-1)*pow(pik,-1))) + 
          (1 + di1)*(1 + dj1)*(1 + dj2)*(dij + dji)*(-(p1k*p2pj*pik) + p1k*p2k*pipj - p1pi*p2k*pjk + p1p2*pik*pjk)*v34*
           ((1 + di2)*log(2) - 2*ai2*log(ai2) - (1 + di2)*log(p2pi*pow(p2k,-1)*pow(pik,-1))) + 
          (1 + di1)*(1 + di2)*(1 + dj2)*(dij + dji)*(-(p1k*p2pj*pik) + p1k*p2k*pipj - p1pi*p2k*pjk + p1p2*pik*pjk)*v34*
           ((1 + dj1)*log(2) - 2*aj1*log(aj1) - (1 + dj1)*log(p1pj*pow(p1k,-1)*pow(pjk,-1))) - 
          (1 + di1)*(1 + di2)*(1 + dj1)*(dij + dji)*(-(p1pj*p2k*pik) + p1k*p2k*pipj - p1k*p2pi*pjk + p1p2*pik*pjk)*v34*
           ((1 + dj2)*log(2) - 2*aj2*log(aj2) - (1 + dj2)*log(p2pj*pow(p2k,-1)*pow(pjk,-1))) + 
          (1 + di1)*(1 + di2)*(1 + dj1)*(1 + dj2)*(p1pj*p2k*pik - p1k*p2pj*pik - p1pi*p2k*pjk + p1k*p2pi*pjk)*
           (-2*(dij + dji)*log(v34) + (dij + dji)*(1 + 2*v34)*log(vp) + (dij + dji)*(-1 + 2*v34)*log((pipj*pow(pik,-1)*pow(pjk,-1))/2.) + 
             2*(aij - vm)*log(aij*pow(vp,-1)) + 2*(aji - vm)*log(aji*pow(vp,-1)))) + 
       3*(1 + di1)*(1 + di2)*(1 + dj1)*(1 + dj2)*(dij + dji)*(p1pj*p2k*pik - p1k*p2pj*pik - p1pi*p2k*pjk + p1k*p2pi*pjk)*(1 + v34)*pow(Pi,2) + 
       2*(1 + di1)*(1 + di2)*(1 + dj1)*(1 + dj2)*(dij + dji)*(p1pj*p2k*pik - p1k*p2pj*pik - p1pi*p2k*pjk + p1k*p2pi*pjk)*(1 + v34)*
        (12*log(2)*log(1 - pow(cosp,2)) - pow(Pi,2) + 12*pow(log(2),2) + 3*pow(log(1 - pow(cosp,2)),2)) + 
       (1 + di1)*(1 + di2)*(1 + dj1)*(1 + dj2)*(dij + dji)*(-(p1pj*p2k*pik) + p1k*p2pj*pik + p1pi*p2k*pjk - p1k*p2pi*pjk)*v34*
        (12*log(2)*log(p1p2*pow(p1k,-1)*pow(p2k,-1)) + pow(Pi,2) - 6*pow(log(2),2) - 6*pow(log(p1p2*pow(p1k,-1)*pow(p2k,-1)),2)) + 
       (1 + di2)*(1 + dj1)*(1 + dj2)*(dij + dji)*(-(p1pj*p2k*pik) + p1k*p2k*pipj - p1k*p2pi*pjk + p1p2*pik*pjk)*v34*
        (24*di1*log(1 - ai1)*log(ai1) - 24*ai1*log(ai1)*log((p1pi*pow(p1k,-1)*pow(pik,-1))/2.) + 12*(1 + di1)*log(2)*log(p1pi*pow(p1k,-1)*pow(pik,-1)) + 
          24*di1*Li2(ai1) + 3*pow(Pi,2) - di1*pow(Pi,2) - 6*(1 + di1)*pow(log(2),2) - 12*ai1*pow(log(ai1),2) - 
          6*(1 + di1)*pow(log(p1pi*pow(p1k,-1)*pow(pik,-1)),2)) - 
       (1 + di1)*(1 + dj1)*(1 + dj2)*(dij + dji)*(-(p1k*p2pj*pik) + p1k*p2k*pipj - p1pi*p2k*pjk + p1p2*pik*pjk)*v34*
        (24*di2*log(1 - ai2)*log(ai2) - 24*ai2*log(ai2)*log((p2pi*pow(p2k,-1)*pow(pik,-1))/2.) + 12*(1 + di2)*log(2)*log(p2pi*pow(p2k,-1)*pow(pik,-1)) + 
          24*di2*Li2(ai2) + 3*pow(Pi,2) - di2*pow(Pi,2) - 6*(1 + di2)*pow(log(2),2) - 12*ai2*pow(log(ai2),2) - 
          6*(1 + di2)*pow(log(p2pi*pow(p2k,-1)*pow(pik,-1)),2)) - 
       (1 + di1)*(1 + di2)*(1 + dj2)*(dij + dji)*(-(p1k*p2pj*pik) + p1k*p2k*pipj - p1pi*p2k*pjk + p1p2*pik*pjk)*v34*
        (24*dj1*log(1 - aj1)*log(aj1) - 24*aj1*log(aj1)*log((p1pj*pow(p1k,-1)*pow(pjk,-1))/2.) + 12*(1 + dj1)*log(2)*log(p1pj*pow(p1k,-1)*pow(pjk,-1)) + 
          24*dj1*Li2(aj1) + 3*pow(Pi,2) - dj1*pow(Pi,2) - 6*(1 + dj1)*pow(log(2),2) - 12*aj1*pow(log(aj1),2) - 
          6*(1 + dj1)*pow(log(p1pj*pow(p1k,-1)*pow(pjk,-1)),2)) + 
       (1 + di1)*(1 + di2)*(1 + dj1)*(dij + dji)*(-(p1pj*p2k*pik) + p1k*p2k*pipj - p1k*p2pi*pjk + p1p2*pik*pjk)*v34*
        (24*dj2*log(1 - aj2)*log(aj2) - 24*aj2*log(aj2)*log((p2pj*pow(p2k,-1)*pow(pjk,-1))/2.) + 12*(1 + dj2)*log(2)*log(p2pj*pow(p2k,-1)*pow(pjk,-1)) + 
          24*dj2*Li2(aj2) + 3*pow(Pi,2) - dj2*pow(Pi,2) - 6*(1 + dj2)*pow(log(2),2) - 12*aj2*pow(log(aj2),2) - 
          6*(1 + dj2)*pow(log(p2pj*pow(p2k,-1)*pow(pjk,-1)),2)) - 
       (1 + di1)*(1 + di2)*(1 + dj1)*(1 + dj2)*(p1pj*p2k*pik - p1k*p2pj*pik - p1pi*p2k*pjk + p1k*p2pi*pjk)*
        (-24*log((pipj*pow(pik,-1)*pow(pjk,-1))/2.)*((dij + dji)*log(v34) - (dij + dji)*(0.5 + v34)*log(vp) + (-aij + vm)*log(aij*pow(vp,-1)) + 
             (-aji + vm)*log(aji*pow(vp,-1))) + (dij + dji)*(3 - 2*v34)*pow(Pi,2) + 12*(dij + dji)*v34*pow(log(vp),2) - 
          12*(dij + dji)*(4*vm*Li2(x) + pow(-2*log(v34) + log(vp),2)/2.) + 
          6*(dij + dji)*(-1 + 2*v34)*pow(log((pipj*pow(pik,-1)*pow(pjk,-1))/2.),2) + 
          12*(v34*log(vp)*(log(aij*pow(vp,-1)) + log(aji*pow(vp,-1))) + (-2*log(v34) + log(vp))*(dij*log(aij*pow(vp,-1)) + dji*log(aji*pow(vp,-1))) - 
             2*(-1 + dij + dji)*v34*(log(aij*pow(vp,-1))*log(1 - aij*pow(vp,-1)) + log(aji*pow(vp,-1))*log(1 - aji*pow(vp,-1)) + 
                Li2(aij*pow(vp,-1)) + Li2(aji*pow(vp,-1))) - (v34*(6*Li2(x) + pow(Pi,2)))/3. + (aij - vm)*pow(log(aij*pow(vp,-1)),2) + 
             (aji - vm)*pow(log(aji*pow(vp,-1)),2))))*pow(1 - pow(cosp,2),-0.5))/48.;
  return res;
}

double psi_contribution_fTTT(double v, double cosp, double y){
  double res = 0;
  //res = (-2*fTTT_m1(v,cosp,-y) - 2*fTTT_m1(v,cosp,y) + (1 + y)*(fTTT_m1(v,cosp,-1) + fTTT_m1(v,cosp,1))*(1 + log(2) - y*log(2)))/(2.*Pi*(1 - y)*(1 + y));
  res = (-2*fTTT_m1(v,cosp,-y) - 2*fTTT_m1(v,cosp,y) + (1 + y)*(fTTT_m1(v,cosp,-1) + fTTT_m1(v,cosp,1))*(1 + log(2) - y*log(2)))/(2.*(1 - y)*(1 + y));
  /*res = res +
    ((fTTT_m2(v,cosp,-1) + fTTT_m2(v,cosp,1))*pow(Pi,-1)*pow(-1 + y,-1)*(6*log(2 - 2*y) + (-1 + y)*(pow(Pi,2) - 9*pow(log(2),2))))/12. - 
     (fTTT_m2(v,cosp,-y) + fTTT_m2(v,cosp,y))*(log(1 - y) + log(1 + y))*pow(Pi,-1)*pow(-1 + pow(y,2),-1);*/
  res = res + ((fTTT_m2(v,cosp,-1) + fTTT_m2(v,cosp,1))*pow(-1 + y,-1)*(-log(64) + 6*log(1 - y) + (-1 + y)*pow(Pi,2) - 3*pow(log(2),2) + 3*y*pow(log(2),2)))/12. - 
    (fTTT_m2(v,cosp,-y) + fTTT_m2(v,cosp,y))*log(-0.25*((-1 + y)*(1 + y)))*pow(-1 + pow(y,2),-1);
  return res;
}

double single_integrand_fTTT(double v, double cosp, double y){
  double res = 0;
  res = fTTT_0(v, cosp, 1)/2 + fTTT_0(v, cosp, -1)/2 + psi_contribution_fTTT(v, cosp, y);
  return res;
}

double fTTT_integrand (double *k, size_t dim, void *params)
{
  (void)(dim); /* avoid unused parameter warnings */
  (void)(params);

  // u in (1,inf), v in (0,inf), cosp in (-1,1), y in (0,1)
  double v = -log(k[0]), cosp = 2*k[1]-1, y = k[2];
  double jac = 1/k[0]*2;

  // now I should add here as well the new terms which actually depend on y
  double res = single_integrand_fTTT(v,cosp,y) * jac;

  if (isnan(res) || isinf(res)){
    //cout << "problems" << endl;
    res = 0;
  }
  return res;
}

void numerical_integrals_fTTT (double *result, double *error)
{

  if (acc_goal_fTTT<0) {
    *result = 0.;
    *error  = 0.;
    return;
  }
  
  double acc_goal = acc_goal_fTTT;
  int max_iter = 20;
  
  double res, err, res1, err1;

  double xl_3[3] = { 0, 0, 0 };
  double xu_3[3] = { 1, 1, 1 };

  const gsl_rng_type *T1;
  gsl_rng *r1;

  gsl_monte_function G1 = { &fTTT_integrand, 3, 0 };

  size_t calls = calls_fTTT;

  gsl_rng_env_setup ();

  T1 = gsl_rng_default;
  r1 = gsl_rng_alloc (T1);

  {
    gsl_monte_vegas_state *s1 = gsl_monte_vegas_alloc (3);

    gsl_monte_vegas_params params1;
    
    for(int i=0; i<2; i++){ // 2 warm-up iterations
      int fac[2] = { 3, 1 };
      gsl_monte_vegas_integrate (&G1, xl_3, xu_3, 3, calls/fac[i], r1, s1, &res1, &err1);
    }
    res = res1;
    err = err1;
    
    if(verbose){display_results ("vegas warm-up", res, err);}
    if(verbose){printf ("converging...\n");}

    gsl_monte_vegas_params_get(s1, &params1);

    int j = 0;
    do
      {
	gsl_monte_vegas_integrate (&G1, xl_3, xu_3, 3, calls, r1, s1, &res1, &err1);
	// after first main iteration change stage to 2
	params1.stage = 2;
	gsl_monte_vegas_params_set(s1, &params1);
	res = res1;
	err = err1;
	j = j+1;
	if(verbose){printf ("result = % .6f sigma = % .6f\n", res, err);}
      }
    //while(j<5);
    while ((abs(err/res*100) > acc_goal && j < max_iter ));
    //while (fabs (gsl_monte_vegas_chisq (s4) - 1.0) > 0.5 && j < max_iter);
    
    if(verbose){display_results ("vegas final", res, err);}

    *result = res;
    *error = err;
    
    gsl_monte_vegas_free (s1);
  }

  gsl_rng_free (r1);

}

double Fex2_fTTT (void)
{
  double numerical_integrals, err;    
  numerical_integrals_fTTT(&numerical_integrals, &err);
  // add here analytic piece. I need to evaluate v34 for this
  v34 = sqrt(1-mi2*mj2/pipj/pipj);
  double res = numerical_integrals + (log(pow(p1pi,2)*pow(p1pj,-2)*pow(p2pi,-2)*pow(p2pj,2))*pow(Pi,3)*(1 + pow(v34,-1)))/24.;
  if(print_error){std::cout << "results fTTT: " << res << ", " << err << std::endl;} // error coming mostly from MC integration
  return res;
}
