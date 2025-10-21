
double Fex2_TITj_nf(double B){

  double result = (56. + 3.*pow(Pi,2) 
		   - 36.*Li3(B/(1. + B)) + 60.*log(1. + B) + 
		   18.*pow(log(1. + B),2) + 6.*pow(log(1. + B),3) + 
		   12.*Li2(-B)*(2. + 3.*log(1. + B)))/216.;
  return result;

}

std::complex<double> arccoth(double x){return -mylog_dp(1. - 1./x)/2. + mylog_dp(1. + 1./x)/2.;}


double Li22_integrand(double t, void *params){
  double z = *(double *) params;
  return -1./2.*log(t)*pow(log(1.-z*t),2)/t;
}

double Li22(double z){
  gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
  double result, error;

  gsl_function F;
  F.function = &Li22_integrand;
  F.params = &z;

  gsl_integration_qags (&F, 0, 1, 0, 1e-5, 1000,
                        w, &result, &error);
  
  gsl_integration_workspace_free (w);

  return result;
}


double Fex2_TITj_CA(double B){

  std::complex<double> result;
  std::complex<double> im(0.,1.);

  result = -1.60787849421845561 - 85.2672608708245*im - (51*z3_dp)/8. + (5.*Li2_c(-B))/9. + (3.*Li2_c(-pow(B,-1)))/2. - 
   Li2_c(pow(1 + B,-1)) + (Li2_c(-pow(B,-1))*Li2_c(B*pow(1 + B,-1)))/4. + (91.*Li2_c(pow(1 + B,-0.5)))/12. + (7.*Li2_c((1 - pow(1 + B,0.5))/2.))/3. - 
   (103.*Li2_c(1 - pow(1 + B,0.5)))/12. - (im*M_PI*Li2_c(1 - pow(1 + B,0.5)))/2. - 2.*Li2_c(-pow(1 + B,0.5)) + (11.*im*M_PI*Li2_c(-pow(1 + B,0.5)))/2. - 3.*Li2_c(pow(1 + B,0.5)) - 
   (7.*Li2_c(2*pow(1 + pow(1 + B,0.5),-1)))/3. + Li3_c(-B) - Li3_c(-pow(B,-1))/2. + Li3_c(1 + pow(B,-1))/2. - Li3_c(pow(1 + B,-1)) + (11.*Li3_c(B*pow(1 + B,-1)))/12. + 
   4.*Li3_c(1 - pow(1 + B,-0.5)) - 4.*Li3_c(0.5 - pow(1 + B,-0.5)/2.) + 8.*Li3_c(pow(1 + B,-0.5)) - 4.*Li3_c((1 + pow(1 + B,-0.5))/2.) + 4.*Li3_c((1 - pow(1 + B,0.5))/2.) + 
   3.*Li3_c(1 - pow(1 + B,0.5)) + 5.*Li3_c(-pow(1 + B,0.5)) - 3.*Li3_c(pow(1 + B,0.5)) - 5.*Li3_c(1 + pow(1 + B,0.5)) - 
   3.*Li3_c((1 + pow(1 + B,0.5))*pow(1 - pow(1 + B,0.5),-1)) - 5.*Li3_c((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)) + 8.*Li3_c(pow(1 + pow(1 + B,0.5),-1)) + 
   4.*Li3_c(2*pow(1 + pow(1 + B,0.5),-1)) + 8.*Li3_c((-1 + pow(1 + B,0.5))*pow(1 + pow(1 + B,0.5),-1)) + 4.*Li3_c(pow(1 + B,0.5)*pow(1 + pow(1 + B,0.5),-1)) - 
   (3.*Li4_c(pow(1 + B,-1)))/4. - Li4_c(B*pow(1 + B,-1)) - mylog_dp(B)/8. + (z3_dp*mylog_dp(B))/2. - Li2_c(-B)*mylog_dp(B) + (Li2_c(-pow(B,-1))*mylog_dp(B))/2. + 
   (3.*Li2_c(pow(1 + B,-1))*mylog_dp(B))/2. - (5.*Li2_c(1 - pow(1 + B,0.5))*mylog_dp(B))/2. + (5.*Li2_c(-pow(1 + B,0.5))*mylog_dp(B))/2. - (Li3_c(-pow(B,-1))*mylog_dp(B))/2. - 
   (73.*mylog_dp(1 + B))/36. - (z3_dp*mylog_dp(1 + B))/2. + (7.*Li2_c(-B)*mylog_dp(1 + B))/12. + (Li2_c(-pow(B,-1))*mylog_dp(1 + B))/2. - Li2_c(pow(1 + B,-1))*mylog_dp(1 + B) + 
   2.*Li2_c(pow(1 + B,-0.5))*mylog_dp(1 + B) + (9.*Li2_c(1 - pow(1 + B,0.5))*mylog_dp(1 + B))/8. - (Li2_c(-pow(1 + B,0.5))*mylog_dp(1 + B))/2. + (3.*Li2_c(pow(1 + B,0.5))*mylog_dp(1 + B))/2. + 
   (Li3_c(-pow(B,-1))*mylog_dp(1 + B))/2. - (Li3_c(pow(1 + B,-1))*mylog_dp(1 + B))/4. - (Li3_c(B*pow(1 + B,-1))*mylog_dp(1 + B))/2. - mylog_dp(B)*mylog_dp(1 + B) + im*M_PI*mylog_dp(B)*mylog_dp(1 + B) - 
   (3.*Li2_c(1 + pow(B,-1))*mylog_dp(1 + pow(B,-1)))/2. + (Li2_c(pow(1 + B,-1))*mylog_dp(1 + pow(B,-1)))/2. - 2.*Li2_c(pow(1 + B,-0.5))*mylog_dp(B*pow(1 + B,-1)) - 
   2.*Li2_c(-pow(1 + B,0.5))*mylog_dp(B*pow(1 + B,-1)) - (4.*mylog_dp(1 - pow(1 + B,-0.5)))/3. + 4.*Li2_c(pow(1 + B,-0.5))*mylog_dp(1 - pow(1 + B,-0.5)) + 
   (Li2_c(1 - pow(1 + B,0.5))*mylog_dp(1 - pow(1 + B,-0.5)))/4. + 2.*Li2_c(-pow(1 + B,0.5))*mylog_dp(1 - pow(1 + B,-0.5)) - (mylog_dp(1 + B)*mylog_dp(1 - pow(1 + B,-0.5)))/6. - 
   4.*Li2_c(pow(1 + B,-0.5))*mylog_dp(pow(1 + B,-0.5)/2.) + 4.*Li2_c((1 + pow(1 + B,-0.5))/2.)*mylog_dp((1 + pow(1 + B,-0.5))/2.) + 
   4.*Li2_c(pow(1 + B,0.5)*pow(1 + pow(1 + B,0.5),-1))*mylog_dp((1 + pow(1 + B,-0.5))/2.) + 2.*Li2_c(1 - pow(1 + B,0.5))*mylog_dp(1 + pow(1 + B,-0.5)) + 
   2.*Li2_c(-pow(1 + B,0.5))*mylog_dp(1 + pow(1 + B,-0.5)) + 2.*mylog_dp(B)*mylog_dp(1 + B)*mylog_dp(1 + pow(1 + B,-0.5)) + (Li2_c(-B)*mylog_dp(1 - pow(1 + B,0.5)))/4. + 
   5.*Li2_c(1 + pow(1 + B,0.5))*mylog_dp(1 - pow(1 + B,0.5)) - (3.*mylog_dp(1 + B)*mylog_dp(1 - pow(1 + B,0.5)))/2. - 4.*Li2_c((1 - pow(1 + B,0.5))/2.)*mylog_dp((-1 + pow(1 + B,0.5))/2.) + 
   (35.*mylog_dp(-1 + pow(1 + B,0.5)))/24. + Li2_c(pow(1 + B,-1))*mylog_dp(-1 + pow(1 + B,0.5)) - 6.*Li2_c(pow(1 + B,-0.5))*mylog_dp(-1 + pow(1 + B,0.5)) - 
   4.*Li2_c(1 - pow(1 + B,0.5))*mylog_dp(-1 + pow(1 + B,0.5)) + 3.*Li2_c(-pow(1 + B,0.5))*mylog_dp(-1 + pow(1 + B,0.5)) + 4.*Li2_c(pow(1 + pow(1 + B,0.5),-1))*mylog_dp(-1 + pow(1 + B,0.5)) + 
   (mylog_dp(2)*mylog_dp(-1 + pow(1 + B,0.5)))/3. - mylog_dp(512)*mylog_dp(-1 + pow(1 + B,0.5)) - mylog_dp(B)*mylog_dp(-1 + pow(1 + B,0.5)) + (13.*mylog_dp(1 + B)*mylog_dp(-1 + pow(1 + B,0.5)))/24. - 
   mylog_dp(4)*mylog_dp(1 + B)*mylog_dp(-1 + pow(1 + B,0.5)) + (13.*mylog_dp(B*pow(1 + B,-1))*mylog_dp(-1 + pow(1 + B,0.5)))/6. + (13.*mylog_dp(1 - pow(1 + B,-0.5))*mylog_dp(-1 + pow(1 + B,0.5)))/6. - 
   (13.*mylog_dp(1 + pow(1 + B,-0.5))*mylog_dp(-1 + pow(1 + B,0.5)))/6. + (5.*im*M_PI*mylog_dp(2)*mylog_dp(-pow(1 + B,0.5)))/2. + (mylog_dp(32)*mylog_dp(1 + B)*mylog_dp(-pow(1 + B,0.5)))/4. - 
   4.*Li2_c(-pow(1 + B,0.5))*mylog_dp(2*pow(1 + B,0.5)) + 4.*Li2_c(2*pow(1 + pow(1 + B,0.5),-1))*mylog_dp((1 + pow(1 + B,0.5))/2.) + 2.*mylog_dp(2)*mylog_dp(1 + B)*mylog_dp((1 + pow(1 + B,0.5))/2.) + 
   4.*mylog_dp(2)*mylog_dp((1 + pow(1 + B,-0.5))/2.)*mylog_dp((1 + pow(1 + B,0.5))/2.) - (19.*mylog_dp(-1 + pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))/2.))/3. - 
   2.*mylog_dp(1 + B)*mylog_dp(-1 + pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))/2.) + mylog_dp(1 + pow(1 + B,0.5))/8. - (9.*Li2_c(-B)*mylog_dp(1 + pow(1 + B,0.5)))/4. - 
   Li2_c(pow(1 + B,-1))*mylog_dp(1 + pow(1 + B,0.5)) + 6.*Li2_c(pow(1 + B,-0.5))*mylog_dp(1 + pow(1 + B,0.5)) + (11.*Li2_c(1 - pow(1 + B,0.5))*mylog_dp(1 + pow(1 + B,0.5)))/4. - 
   7.*Li2_c(-pow(1 + B,0.5))*mylog_dp(1 + pow(1 + B,0.5)) + 4.*Li2_c(pow(1 + pow(1 + B,0.5),-1))*mylog_dp(1 + pow(1 + B,0.5)) - mylog_dp(B)*mylog_dp(1 + pow(1 + B,0.5)) + 
   (19.*mylog_dp(1 + B)*mylog_dp(1 + pow(1 + B,0.5)))/6. + (11.*im*M_PI*mylog_dp(1 + B)*mylog_dp(1 + pow(1 + B,0.5)))/4. - (3.*mylog_dp(B)*mylog_dp(1 + B)*mylog_dp(1 + pow(1 + B,0.5)))/4. - 
   (13.*mylog_dp(B*pow(1 + B,-1))*mylog_dp(1 + pow(1 + B,0.5)))/6. + (13.*mylog_dp(1 - pow(1 + B,-0.5))*mylog_dp(1 + pow(1 + B,0.5)))/2. - 
   mylog_dp(1 + B)*mylog_dp(1 - pow(1 + B,-0.5))*mylog_dp(1 + pow(1 + B,0.5)) + (13.*mylog_dp(1 + pow(1 + B,-0.5))*mylog_dp(1 + pow(1 + B,0.5)))/6. + 
   mylog_dp(1 + B)*mylog_dp(1 + pow(1 + B,-0.5))*mylog_dp(1 + pow(1 + B,0.5)) + 5.*im*M_PI*mylog_dp(1 - pow(1 + B,0.5))*mylog_dp(1 + pow(1 + B,0.5)) + 
   (7.*mylog_dp((-1 + pow(1 + B,0.5))/2.)*mylog_dp(1 + pow(1 + B,0.5)))/3. - (7.*mylog_dp(1 + B)*mylog_dp(-1 + pow(1 + B,0.5))*mylog_dp(1 + pow(1 + B,0.5)))/2. - 
   5.*Li2_c(1 - pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))*pow(1 - pow(1 + B,0.5),-1)) + 5.*Li2_c(1 + pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))*pow(1 - pow(1 + B,0.5),-1)) - 
   5.*Li2_c((1 + pow(1 + B,0.5))*pow(1 - pow(1 + B,0.5),-1))*mylog_dp((1 + pow(1 + B,0.5))*pow(1 - pow(1 + B,0.5),-1)) + 
   5.*Li2_c((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1))*mylog_dp((1 + pow(1 + B,0.5))*pow(1 - pow(1 + B,0.5),-1)) - 
   (5.*mylog_dp(1 + B)*mylog_dp(1 + pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))*pow(1 - pow(1 + B,0.5),-1)))/2. + 
   5.*mylog_dp(-pow(1 + B,0.5))*mylog_dp(1 + pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))*pow(1 - pow(1 + B,0.5),-1)) + 
   4.*Li2_c(1 - pow(1 + B,-0.5))*mylog_dp(2*pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1)) - 4.*Li2_c(0.5 - pow(1 + B,-0.5)/2.)*mylog_dp(2*pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1)) - 
   2.*mylog_dp(2)*mylog_dp(1 + B)*mylog_dp(2*pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1)) - 2.*Li2_c(pow(1 + B,-1))*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)) + 
   4.*Li2_c(pow(1 + B,-0.5))*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)) - 8.*Li2_c(1 - pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)) - 
   4.*Li2_c(-pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)) - 
   4.*Li2_c(pow(1 + pow(1 + B,0.5),-1))*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)) - (13.*mylog_dp(B)*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)))/2. + 
   mylog_dp(B)*mylog_dp(1 + B)*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)) + (13.*mylog_dp(-1 + pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)))/2. - 
   mylog_dp(1 + B)*mylog_dp(-1 + pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)) + 
   (13.*mylog_dp(1 + pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)))/2. - 
   mylog_dp(1 + B)*mylog_dp(1 + pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)) + 
   8.*mylog_dp(-1 + pow(1 + B,0.5))*mylog_dp(1 + pow(1 + B,0.5))*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)) - 
   (9.*Li2_c(1 - pow(1 + B,0.5))*mylog_dp((-1 + pow(1 + B,0.5))*pow(1 + pow(1 + B,0.5),-1)))/4. + Li22(-pow(B,-1))/2. - (Li2_c(pow(1 + B,-1))*pow(1 + B,-0.5))/4. + 
   Li2_c(pow(1 + B,-0.5))*pow(1 + B,-0.5) - (2.*mylog_dp(B*pow(1 + B,-1))*pow(1 + B,-0.5))/3. + (2.*mylog_dp(1 - pow(1 + B,-0.5))*pow(1 + B,-0.5))/3. - 
   (mylog_dp(1 + B)*mylog_dp(1 - pow(1 + B,-0.5))*pow(1 + B,-0.5))/4. + (2.*mylog_dp(1 + pow(1 + B,-0.5))*pow(1 + B,-0.5))/3. + (mylog_dp(1 + B)*mylog_dp(1 + pow(1 + B,-0.5))*pow(1 + B,-0.5))/4. + 
    (2.*arccoth(1 + 2*B)*pow(M_PI,2))/3. + (5.*Li2_c(-pow(B,-1))*pow(M_PI,2))/24. + (3.*Li2_c(pow(1 + B,-1))*pow(M_PI,2))/8. + (Li2_c(B*pow(1 + B,-1))*pow(M_PI,2))/12. + 
   (mylog_dp(B)*pow(M_PI,2))/4. + (23.*mylog_dp(1 + B)*pow(M_PI,2))/12. - (5.*mylog_dp(B)*mylog_dp(1 + B)*pow(M_PI,2))/4. - (mylog_dp(1 + B)*mylog_dp(1 + pow(B,-1))*pow(M_PI,2))/6. + 
   (2.*mylog_dp(B*pow(1 + B,-1))*pow(M_PI,2))/3. - (4.*mylog_dp(1 - pow(1 + B,-0.5))*pow(M_PI,2))/3. + (4.*mylog_dp(pow(1 + B,-0.5)/2.)*pow(M_PI,2))/3. - 
   (4.*mylog_dp((1 + pow(1 + B,-0.5))/2.)*pow(M_PI,2))/3. + (4.*mylog_dp(1 + pow(1 + B,-0.5))*pow(M_PI,2))/3. + (17.*mylog_dp(1 - pow(1 + B,0.5))*pow(M_PI,2))/24. - 
   (35.*mylog_dp(-1 + pow(1 + B,0.5))*pow(M_PI,2))/12. - (17.*mylog_dp(1 + pow(1 + B,0.5))*pow(M_PI,2))/24. - (mylog_dp(2*pow(-1 + pow(1 + B,0.5),-1))*pow(M_PI,2))/4. - 
   (11.*mylog_dp(pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1))*pow(M_PI,2))/12. - (13.*mylog_dp(2*pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1))*pow(M_PI,2))/12. - 
   (8.*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1))*pow(M_PI,2))/3. - pow(Li2_c(B*pow(1 + B,-1)),2)/8. + (11.*pow(mylog_dp(B),2))/4. - (im*M_PI*pow(mylog_dp(B),2))/2. + 
   (Li2_c(-pow(B,-1))*pow(mylog_dp(B),2))/8. - (mylog_dp(1 + B)*pow(mylog_dp(B),2))/2. + (25*pow(M_PI,2)*pow(mylog_dp(B),2))/48. + pow(mylog_dp(B),3)/4. - (5.*mylog_dp(1 + B)*pow(mylog_dp(B),3))/12. + 
   (25.*pow(mylog_dp(B),4))/96. - (3.*pow(mylog_dp(1 + B),2))/32. - (im*M_PI*pow(mylog_dp(1 + B),2))/8. - (Li2_c(-pow(B,-1))*pow(mylog_dp(1 + B),2))/2. - 
   (3.*Li2_c(pow(1 + B,-1))*pow(mylog_dp(1 + B),2))/8. + (mylog_dp(2)*pow(mylog_dp(1 + B),2))/2. + (mylog_dp(B)*pow(mylog_dp(1 + B),2))/8. - (mylog_dp(B*pow(1 + B,-1))*pow(mylog_dp(1 + B),2))/4. + 
   (mylog_dp(1 - pow(1 + B,-0.5))*pow(mylog_dp(1 + B),2))/2. - (mylog_dp(pow(1 + B,-0.5)/2.)*pow(mylog_dp(1 + B),2))/2. + 2.*mylog_dp(-1 + pow(1 + B,0.5))*pow(mylog_dp(1 + B),2) + 
   (5.*mylog_dp(1 + pow(1 + B,0.5))*pow(mylog_dp(1 + B),2))/2. - (mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1))*pow(mylog_dp(1 + B),2))/2. + (11*pow(M_PI,2)*pow(mylog_dp(1 + B),2))/16. - 
   (pow(mylog_dp(B),2)*pow(mylog_dp(1 + B),2))/8. - (29.*pow(mylog_dp(1 + B),3))/72. + (11.*mylog_dp(B)*pow(mylog_dp(1 + B),3))/24. - (17.*pow(mylog_dp(1 + B),4))/96. + 
   (Li2_c(-pow(B,-1))*pow(mylog_dp(1 + pow(B,-1)),2))/8. - (Li2_c(B*pow(1 + B,-1))*pow(mylog_dp(1 + pow(B,-1)),2))/8. + (mylog_dp(B)*pow(mylog_dp(1 + pow(B,-1)),2))/4. - 
   (3.*mylog_dp(-pow(B,-1))*pow(mylog_dp(1 + pow(B,-1)),2))/4. - (7*pow(M_PI,2)*pow(mylog_dp(1 + pow(B,-1)),2))/24. + pow(mylog_dp(1 + pow(B,-1)),3)/4. + 
   (mylog_dp(1 + B)*pow(mylog_dp(1 + pow(B,-1)),3))/12. - (7.*pow(mylog_dp(1 + pow(B,-1)),4))/96. - 2.*mylog_dp(2)*pow(mylog_dp((1 + pow(1 + B,-0.5))/2.),2) + 
   2.*mylog_dp((1 + pow(1 + B,-0.5))/2.)*pow(mylog_dp(1 + pow(1 + B,-0.5)),2) - (2.*pow(mylog_dp(1 + pow(1 + B,-0.5)),3))/3. - (16.*pow(mylog_dp(-1 + pow(1 + B,0.5)),2))/3. - 
   4.*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1))*pow(mylog_dp(-1 + pow(1 + B,0.5)),2) - (5.*mylog_dp(2)*pow(mylog_dp(-pow(1 + B,0.5)),2))/2. - 
   (7.*pow(mylog_dp((1 + pow(1 + B,0.5))/2.),2))/6. - 2.*mylog_dp(2)*pow(mylog_dp((1 + pow(1 + B,0.5))/2.),2) + (4.*pow(mylog_dp((1 + pow(1 + B,0.5))/2.),3))/3. - pow(mylog_dp(1 + pow(1 + B,0.5)),2) - 
   (13.*mylog_dp(1 + B)*pow(mylog_dp(1 + pow(1 + B,0.5)),2))/4. + 2.*mylog_dp(-1 + pow(1 + B,0.5))*pow(mylog_dp(1 + pow(1 + B,0.5)),2) - 
   (5.*mylog_dp(-pow(1 + B,0.5))*pow(mylog_dp(1 + pow(1 + B,0.5)),2))/2. - 6.*mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1))*pow(mylog_dp(1 + pow(1 + B,0.5)),2) + 
   (2.*pow(mylog_dp(1 + pow(1 + B,0.5)),3))/3. + (5.*im*M_PI*pow(mylog_dp((1 + pow(1 + B,0.5))*pow(1 - pow(1 + B,0.5),-1)),2))/2. + 
   2.*mylog_dp(2*pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1))*pow(mylog_dp(pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1)),2) - 
   (2.*pow(mylog_dp(pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1)),3))/3. - 2.*mylog_dp(-1 + pow(1 + B,0.5))*pow(mylog_dp(2*pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1)),2) - 
   2.*mylog_dp(pow(-1 + pow(1 + B,0.5),-1)/2.)*pow(mylog_dp(2*pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1)),2) - (4.*pow(mylog_dp(2*pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1)),3))/3. - 
   mylog_dp(1 + B)*pow(mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)),2) + (9.*mylog_dp(-1 + pow(1 + B,0.5))*pow(mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)),2))/4. + 
   (mylog_dp(2*pow(-1 + pow(1 + B,0.5),-1))*pow(mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)),2))/4. + 
   (9.*mylog_dp(pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1))*pow(mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)),2))/4. - 
   (mylog_dp(2*pow(1 + B,0.5)*pow(-1 + pow(1 + B,0.5),-1))*pow(mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)),2))/4. + 
   (8.*pow(mylog_dp((1 + pow(1 + B,0.5))*pow(-1 + pow(1 + B,0.5),-1)),3))/3.;

  return result.real();

}


double Fex2_TITj(double B, double A){
  // Term in <F_{ex,2}>^{(0)} prop to TITj
  // Here B=pT2/m2 for I and A=4*pIpj^2/m2/s
  return CA*Fex2_TITj_CA(B) + nf*Fex2_TITj_nf(B)
         + ((-2.*nf*(56. + 3.*pow(Pi,2))
	 + CA*(808. + 33.*pow(Pi,2) - 756.*z3))*log(A))/432.;
}

double Fex1_1_TITj(double B, double A){
  // Term in <F_{ex,1}>^{(1)} prop to TITj
  // Here B=pT2/m2 for I and A=4*pIpj^2/m2/s
  return ((-pow(M_PI,2) + pow(M_PI,2)*log(A) + 12.*Li2_c(-B) - 12.*Li3_c(-B))/24.).real();
}
