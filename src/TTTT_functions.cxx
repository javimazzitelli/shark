double Fex1_TITj_phi(double phi, void *params){
  // here x = b.pjT/mj = sqrt(B)*cosphi in the ttbar case,
  // that is pTi/mi*cosphi and pTj/mj*cosphi=-pTi/mi*cosphi
  // Now this should change to
  // pTi/mi*cosphi and pTj/mj*(cosxi*cosphi+sinxi*sinphi)
  // I take the expression from eq 25 in 2102.03256
  // Careful, that object (D1) is normalized to zero.
  // I put directly the object normalized to Fex1
  evaluate_external(false);
  int i = *(int *) params;
  double x;
  if (i==3){
    x = pTi/sqrt(mi2)*cos(phi);
  } else if(i==4){
    x = pTj/sqrt(mj2)*(cos(phi)*cosxi+sin(phi)*sinxi);
  }
  double result = pow(asinh(x),2) - x/sqrt(1.+x*x)*asinh(x);
  return result;
}

double Fex1_TITj_phi_im(double phi, void *params){
  evaluate_external(false);
  int i = *(int *) params;
  double x;
  if (i==3){
    x = pTi/sqrt(mi2)*cos(phi);
  } else if(i==4){
    x = pTj/sqrt(mj2)*(cos(phi)*cosxi+sin(phi)*sinxi);
  }
  double result = asinh(x)*M_PI - x/sqrt(1.+x*x)*M_PI/2;
  return result;
}

double Fex1_TITj_phi_avg(int i){
  gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
  double result, error;

  gsl_function F;
  F.function = &Fex1_TITj_phi;
  F.params = &i;

  gsl_integration_qags (&F, 0, 2*M_PI, 0, 1e-5, 1000,
                        w, &result, &error);
  
  gsl_integration_workspace_free (w);

  return result/M_PI/2;
}

double aux_integrand_D34(double x, void *params){
  double phi = *(double *) params;
  // following eq 26
  evaluate_external(false);
  double bwT = x*pTi*cos(phi) + (1-x)*pTj*(cos(phi)*cosxi+sin(phi)*sinxi);
  double w2  = pow(x,2)*mi2 + pow(1-x,2)*mj2 + 2*x*(1-x)*pipj;
  double result = pipj/w2*2*bwT/sqrt(w2+bwT*bwT)*asinh(bwT/sqrt(w2));
  return result;
}

double aux_integral_D34(double phi){
  gsl_integration_workspace * ww  = gsl_integration_workspace_alloc (1000);
  double result, error;

  gsl_function FF;
  FF.function = &aux_integrand_D34;
  FF.params = &phi;

  // This call was generating a call to the gsl error handler.
  // Now I implement a rescue for those situations
  gsl_set_error_handler_off();
  int status;
  
  status = gsl_integration_qags (&FF, 0, 1, 0, 1e-5, 1000,
                        ww, &result, &error);
  
  gsl_integration_workspace_free (ww);

  if (status) {
    // the attempt to integrate failed, I try by splitting the integration region
    double old_result = result;
    if (verbose) {
      std::cout << "shark: aux_integral_D34: first attempt failed" << std::endl;
      std::cout << "shark: aux_integral_D34: old result = " << old_result << std::endl;
    }
    
    gsl_integration_workspace * ww1  = gsl_integration_workspace_alloc (1000);
    gsl_integration_workspace * ww2  = gsl_integration_workspace_alloc (1000);
    gsl_integration_workspace * ww3  = gsl_integration_workspace_alloc (1000);
    double result1, error1,result2, error2,result3, error3;
    int status1,status2,status3;
    
    status1 = gsl_integration_qags (&FF, 0, 0.1, 0, 1e-5, 1000,
			  ww1, &result1, &error1);
    status2 = gsl_integration_qags (&FF, 0.1, 0.9, 0, 1e-5, 1000,
			  ww2, &result2, &error2);
    status3 = gsl_integration_qags (&FF, 0.9, 1, 0, 1e-5, 1000,
			  ww3, &result3, &error3);
    result = result1 + result2 + result3;
    if(verbose) {
      std::cout << "shark: aux_integral_D34: new result = " << result << std::endl;
      std::cout << "shark: aux_integral_D34: ratio = " << result/old_result << std::endl;
    }

    if (status1 || status2 || status3){
      if (verbose){std::cout << "shark: aux_integral_D34: second attempt failed" << std::endl;}
      // if old and new result agree within a given tolerance, accept the point nevertheless
      if (abs(result/old_result - 1) < 0.01){ // 1% tolerance
	if (verbose){std::cout << "shark: aux_integral_D34: point accepted" << std::endl;}
      }
      else {
	if (verbose){std::cout << "shark: aux_integral_D34: point rejected" << std::endl;}
	shark_error_handler("problems in aux_integral_D34","dummy",0,0);
      }
    }
    
    gsl_integration_workspace_free (ww1);
    gsl_integration_workspace_free (ww2);
    gsl_integration_workspace_free (ww3);
  }

  gsl_set_error_handler (&shark_error_handler);
  return result;
}

double aux_integrand_D34_im(double x, void *params){
  double phi = *(double *) params;
  // following eq 26
  evaluate_external(false);
  double bwT = x*pTi*cos(phi) + (1-x)*pTj*(cos(phi)*cosxi+sin(phi)*sinxi);
  double w2  = pow(x,2)*mi2 + pow(1-x,2)*mj2 + 2*x*(1-x)*pipj;
  double result = pipj/w2*2*bwT/sqrt(w2+bwT*bwT)*M_PI/2;
  return result;
}

double aux_integral_D34_im(double phi){
  gsl_integration_workspace * ww  = gsl_integration_workspace_alloc (1000);
  double result, error;

  gsl_function FF;
  FF.function = &aux_integrand_D34_im;
  FF.params = &phi;

  // This call was generating a call to the gsl error handler.
  // Now I implement a rescue for those situations
  gsl_set_error_handler_off();
  int status;
  
  status = gsl_integration_qags (&FF, 0, 1, 0, 1e-5, 1000,
                        ww, &result, &error);
  
  gsl_integration_workspace_free (ww);

  if (status) {
    // the attempt to integrate failed, I try by splitting the integration region
    double old_result = result;
    if (verbose) {
      std::cout << "shark: aux_integral_D34_im: first attempt failed" << std::endl;
      std::cout << "shark: aux_integral_D34_im: old result = " << old_result << std::endl;
    }
    
    gsl_integration_workspace * ww1  = gsl_integration_workspace_alloc (1000);
    gsl_integration_workspace * ww2  = gsl_integration_workspace_alloc (1000);
    gsl_integration_workspace * ww3  = gsl_integration_workspace_alloc (1000);
    double result1, error1,result2, error2,result3, error3;
    int status1,status2,status3;
    
    status1 = gsl_integration_qags (&FF, 0, 0.1, 0, 1e-5, 1000,
			  ww1, &result1, &error1);
    status2 = gsl_integration_qags (&FF, 0.1, 0.9, 0, 1e-5, 1000,
			  ww2, &result2, &error2);
    status3 = gsl_integration_qags (&FF, 0.9, 1, 0, 1e-5, 1000,
			  ww3, &result3, &error3);
    result = result1 + result2 + result3;
    if(verbose) {
      std::cout << "shark: aux_integral_D34_im: new result = " << result << std::endl;
      std::cout << "shark: aux_integral_D34_im: ratio = " << result/old_result << std::endl;
    }

    if (status1 || status2 || status3){
      if (verbose){std::cout << "shark: aux_integral_D34_im: second attempt failed" << std::endl;}
      // if old and new result agree within a given tolerance, accept the point nevertheless
      if (abs(result/old_result - 1) < 0.01){ // 1% tolerance
	if (verbose){std::cout << "shark: aux_integral_D34_im: point accepted" << std::endl;}
      }
      else {
	if (verbose){std::cout << "shark: aux_integral_D34_im: point rejected" << std::endl;}
	//shark_error_handler("problems in aux_integral_D34_im","dummy",0,0);
      }
    }
    
    gsl_integration_workspace_free (ww1);
    gsl_integration_workspace_free (ww2);
    gsl_integration_workspace_free (ww3);
  }

  gsl_set_error_handler (&shark_error_handler);
  return result;
}


double Fex1_T3T4_phi(double phi, void *dummy){
  evaluate_external(false);
  double x3 = pTi/sqrt(mi2)*cos(phi);
  double x4 = pTj/sqrt(mj2)*(cos(phi)*cosxi+sin(phi)*sinxi);  
  return aux_integral_D34(phi) - x3/sqrt(1.+x3*x3)*asinh(x3) - x4/sqrt(1.+x4*x4)*asinh(x4);
}

double Fex1_T3T4_phi_im(double phi, void *dummy){
  evaluate_external(false);
  double x3 = pTi/sqrt(mi2)*cos(phi);
  double x4 = pTj/sqrt(mj2)*(cos(phi)*cosxi+sin(phi)*sinxi);  
  return aux_integral_D34_im(phi) - x3/sqrt(1.+x3*x3)*M_PI/2 - x4/sqrt(1.+x4*x4)*M_PI/2;
}

double Fex1_T3T4_phi_avg(){
  gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
  double result, error;

  gsl_function F;
  F.function = &Fex1_T3T4_phi;

  gsl_integration_qags (&F, 0, 2*M_PI, 0, 1e-5, 1000,
                        w, &result, &error);
  
  gsl_integration_workspace_free (w);

  return result/M_PI/2;
}

double Fex1sq_T3T4T3T4_phi(double phi, void *dummy){
  return pow(Fex1_T3T4_phi(phi,dummy),2) - pow(Fex1_T3T4_phi_im(phi,dummy),2);
}

double Fex1sq_T3T4TITj_phi(double phi, void *params){
  return Fex1_T3T4_phi(phi,params)*Fex1_TITj_phi(phi,params) - Fex1_T3T4_phi_im(phi,params)*Fex1_TITj_phi_im(phi,params);
}

double Fex1sq_TITjTKTl_phi(double phi, void *params){
  int aux = *(int *) params;
  int i,k;
  if(aux==33){i=3;k=3;}
  if(aux==34){i=3;k=4;}
  if(aux==43){i=4;k=3;}
  if(aux==44){i=4;k=4;}
  void *params1 = &i;
  void *params2 = &k;
  return Fex1_TITj_phi(phi,params1)*Fex1_TITj_phi(phi,params2) - Fex1_TITj_phi_im(phi,params1)*Fex1_TITj_phi_im(phi,params2);
}

double Fex1sq_T3T4T3T4_phi_avg(){
  gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
  double result, error;

  gsl_function F;
  F.function = &Fex1sq_T3T4T3T4_phi;

  gsl_integration_qags (&F, 0, 2*M_PI, 0, 1e-5, 1000,
                        w, &result, &error);
  
  gsl_integration_workspace_free (w);

  return result/M_PI/2;
}

double Fex1_T3T4_phi_im_avg(){
  gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
  double result, error;

  gsl_set_error_handler (&shark_error_handler);
  
  gsl_function F;
  F.function = &Fex1_T3T4_phi_im;

  gsl_integration_qags (&F, 0, 2*M_PI, 0, 1e-5, 1000,
                        w, &result, &error);
  
  gsl_integration_workspace_free (w);

  return result/M_PI/2;
}

double Fex1sq_T3T4TITj_phi_avg(int i){
  gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
  double result, error;

  gsl_function F;
  F.function = &Fex1sq_T3T4TITj_phi;
  F.params = &i;

  gsl_integration_qags (&F, 0, 2*M_PI, 0, 1e-5, 1000,
                        w, &result, &error);
  
  gsl_integration_workspace_free (w);

  return result/M_PI/2;
}

double Fex1sq_TITjTKTl_phi_avg(int i, int k){
  gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
  double result, error;

  gsl_function F;
  F.function = &Fex1sq_TITjTKTl_phi;
  int aux = i*10 + k;
  F.params = &aux;

  gsl_integration_qags (&F, 0, 2*M_PI, 0, 1e-5, 1000,
                        w, &result, &error);
  
  gsl_integration_workspace_free (w);

  return result/M_PI/2;
}

double sqFex1_T3T4T3T4(void){
  return pow(Fex1_T3T4(),2);
}

double sqFex1_T3T4TITj(int i){
  evaluate_external(false);
  double B;
  if (i==3){B=pTi*pTi/mi2;}
  if (i==4){B=pTj*pTj/mj2;}
  return Fex1_T3T4()*Fex1_TITj(B);
}

double sqFex1_TITjTKTl(int i, int k){
  evaluate_external(false);
  double Bi,Bk;
  if (i==3){Bi=pTi*pTi/mi2;}
  if (i==4){Bi=pTj*pTj/mj2;}
  if (k==3){Bk=pTi*pTi/mi2;}
  if (k==4){Bk=pTj*pTj/mj2;}
  return Fex1_TITj(Bi)*Fex1_TITj(Bk);
}
