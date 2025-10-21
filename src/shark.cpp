//====== SHARK: Soft function for Heavy quark production in ARbitrary Kinematics ======//
//                                     ,-					       //
//                                   ,'::|					       //
//                                  /::::|					       //
//                                ,'::::o\                                      _..    //
//             ____........-------,..::?88b                                  ,-' /     //
//     _.--"""". . . .      .   .  .  .  ""`-._                           ,-' .;'      //
//    <. - :::::o......  ...   . . .. . .  .  .""--._                  ,-'. .;'	       //
//     `-._  ` `":`:`:`::||||:::::::::::::::::.:. .  ""--._ ,'|     ,-'.  .;'	       //
//         """_=--       //'doo.. ````:`:`::::::::::.:.:.:. .`-`._-'.   .;'	       //
//             ""--.__     P(       \               ` ``:`:``:::: .   .;'	       //
//                    "\""--.:-.     `.                             .:/		       //
//                      \. /    `-._   `.""-----.,-..::(--"".\""`.  `:\		       //
//                       `P         `-._ \          `-:\          `. `:\	       //
//                                       ""            "            `-._) 	       //
//                                                                                     //
//=====================================================================================//

#include "../inc/shark.h"
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "polylogarithm/src/cpp/Li2.hpp"
#include "polylogarithm/src/cpp/Li3.hpp"
#include "polylogarithm/src/cpp/Li4.hpp"
#include <complex>

// using shark namespace
namespace shark
{

double p1[4],p2[4],pi[4],pj[4];
double mi2, mj2, pi$minus, pj$minus, pi$plus, pj$plus, pTi, pTj, pipj, cosxi, sinxi, pTipTj, p1pi, p1pj, p2pi, p2pj, p1p2, p1$minus, p2$minus, p1$plus, p2$plus;
double aux$ni2, aux$nj2, aux$ninj, ni2, nj2, ninj, aux$ni2$u0, aux$nj2$u0, aux$ninj$u0;
double Lu, Bij, Bij$u0, sqrtBij, sqrtBij$u0, LBij$minus, LBij$plus, LBij$u0$minus, LBij$u0$plus, Laux$ni2, Lni2, Lcosp, Li2Bij$plus, Li2Bij$minus, Li2Bij$u0$plus, Li2Bij$u0$minus, Laux$ni2$u0, b;
double Laux1, Laux2, Laux3, Laux4, Laux5, Laux6, Laux7, Laux8, Li3Bij$u0$minus, Li3Bij$u0$plus, Paux1, Paux2, Paux3, Paux4, Paux5, Paux6, Paux7;
double alpha, mu, csi, omega, rho, gammaij, Dij, sqrt$ni2, sqrt$aux;
double pik, pjk, p1k, p2k, aij, aji, ai1, ai2, aj1, aj2, v34, vp, vm, x, dij, dji, di1, di2, dj1, dj2;
double Log2 = log(2), Pi = M_PI, z3 = 1.2020569031595942854;
double CA = 3, nf, beta0, z3_dp = 1.2020569031595942854;
bool verbose = false, print_error = false;
double pipj_save1, p1pj_save1, p2pj_save1, pipj_save2, p1pj_save2, p2pj_save2;
double result_save_T3T4 = -1e10, result_save_fTTT = -1e10;
double acc_goal_T3T4 = 5., acc_goal_fTTT = 5.;
int calls_T3T4 = 5000, calls_fTTT = 5000;

void shark_error_handler (const char * reason,
              const char * file,
              int line,
              int gsl_errno){
  std::cout << "shark: gsl error: " << reason << std::endl;
  std::cout << "shark: problematic phase space point:" << std::endl;
  std::cout << "shark: p1 = (" << p1[0] << ", " << p1[1] << ", " << p1[2] << ", " << p1[3] << ")" << std::endl;
  std::cout << "shark: p2 = (" << p2[0] << ", " << p2[1] << ", " << p2[2] << ", " << p2[3] << ")" << std::endl;
  std::cout << "shark: p3 = (" << pi[0] << ", " << pi[1] << ", " << pi[2] << ", " << pi[3] << ")" << std::endl;
  std::cout << "shark: p4 = (" << pj[0] << ", " << pj[1] << ", " << pj[2] << ", " << pj[3] << ")" << std::endl;
  std::cout << "shark: stop" << std::endl;
  exit(-1);
}
  
void evaluate_external(bool cross){
  if (cross) {
    pTj = sqrt(pi[1]*pi[1]+pi[2]*pi[2]);
    pTi = sqrt(pj[1]*pj[1]+pj[2]*pj[2]);
    pj$minus = (pi[0]-pi[3])/sqrt(2);
    pi$minus = (pj[0]-pj[3])/sqrt(2);
    pj$plus  = (pi[0]+pi[3])/sqrt(2);
    pi$plus  = (pj[0]+pj[3])/sqrt(2);
    mj2 = pi[0]*pi[0] - pTj*pTj - pi[3]*pi[3];
    mi2 = pj[0]*pj[0] - pTi*pTi - pj[3]*pj[3];
    pTipTj = pi[1]*pj[1] + pi[2]*pj[2];
    cosxi = pTipTj/pTi/pTj;
    sinxi = sqrt(1-cosxi*cosxi);
    if(abs(cosxi)>1 && abs(cosxi)<1.000001 && isnan(sinxi)){sinxi=0;} //avoid nan when cosxi=+/-1
    pipj = pi$minus*pj$plus + pj$minus*pi$plus - pTipTj;
    p1pj = p1[0]*pi[0]-p1[1]*pi[1]-p1[2]*pi[2]-p1[3]*pi[3];
    p2pj = p2[0]*pi[0]-p2[1]*pi[1]-p2[2]*pi[2]-p2[3]*pi[3];
    p1pi = p1[0]*pj[0]-p1[1]*pj[1]-p1[2]*pj[2]-p1[3]*pj[3];
    p2pi = p2[0]*pj[0]-p2[1]*pj[1]-p2[2]*pj[2]-p2[3]*pj[3];
    p1p2 = p1[0]*p2[0]-p1[1]*p2[1]-p1[2]*p2[2]-p1[3]*p2[3];
  } else {
    pTi = sqrt(pi[1]*pi[1]+pi[2]*pi[2]);
    pTj = sqrt(pj[1]*pj[1]+pj[2]*pj[2]);
    pi$minus = (pi[0]-pi[3])/sqrt(2);
    pj$minus = (pj[0]-pj[3])/sqrt(2);
    pi$plus  = (pi[0]+pi[3])/sqrt(2);
    pj$plus  = (pj[0]+pj[3])/sqrt(2);
    mi2 = pi[0]*pi[0] - pTi*pTi - pi[3]*pi[3];
    mj2 = pj[0]*pj[0] - pTj*pTj - pj[3]*pj[3];
    pTipTj = pi[1]*pj[1] + pi[2]*pj[2];
    cosxi = pTipTj/pTi/pTj;
    sinxi = sqrt(1-cosxi*cosxi);
    if(abs(cosxi)>1 && abs(cosxi)<1.000001 && isnan(sinxi)){sinxi=0;} //avoid nan when cosxi=+/-1
    pipj = pi$minus*pj$plus + pj$minus*pi$plus - pTipTj;
    p1pi = p1[0]*pi[0]-p1[1]*pi[1]-p1[2]*pi[2]-p1[3]*pi[3];
    p2pi = p2[0]*pi[0]-p2[1]*pi[1]-p2[2]*pi[2]-p2[3]*pi[3];
    p1pj = p1[0]*pj[0]-p1[1]*pj[1]-p1[2]*pj[2]-p1[3]*pj[3];
    p2pj = p2[0]*pj[0]-p2[1]*pj[1]-p2[2]*pj[2]-p2[3]*pj[3];
    p1p2 = p1[0]*p2[0]-p1[1]*p2[1]-p1[2]*p2[2]-p1[3]*p2[3];
  }
}

void set_nf(int nf_in){
  nf = nf_in;
  beta0 = (11*CA - 2*nf)/(12*Pi);
}

double Li2(double x){
  double res = polylogarithm::Li2(x);
  return res;
}

double Li3(double x){
  std::complex<double> xx;
  xx = {x,0.};
  double res = (polylogarithm::Li3(xx)).real();
  return res;
}

double Li4(double x){
  std::complex<double> xx;
  xx = {x,0.};
  double res = (polylogarithm::Li4(xx)).real();
  return res;
}

std::complex<double> Li2_c(double x){
  std::complex<double> xx;
  xx = {x,-1.e-40};
  std::complex<double> res = polylogarithm::Li2(xx);
  return res;
}

std::complex<double> Li3_c(double x){
  std::complex<double> xx;
  xx = {x,-1.e-40};
  std::complex<double> res = polylogarithm::Li3(xx);
  return res;
}

std::complex<double> Li4_c(double x){
  std::complex<double> xx;
  xx = {x,-1.e-40};
  std::complex<double> res = polylogarithm::Li4(xx);
  return res;
}

std::complex<double> mylog(double x){
  if (x>0) {
    return {log(x), 0.};
  } else {
    return {log(-x), Pi};
  }
}

std::complex<double> mylog_dp(double x){
  if (x>0) {
    return {log(x), 0.};
  } else {
    return {log(-x), M_PI};
  }
}

double mypow(double x, double y){
  if (x==0 and y<0){
    throw std::runtime_error("Error: dividing by zero\n");
  }
  return pow(x,y);
}

#include "T3T4_functions.cxx"
#include "TITj_functions.cxx"
#include "H1_functions.cxx"
#include "TTTT_functions.cxx"
#include "fTTT_functions.cxx"

double H1_TT(int i, int j, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int nf_in){
  gsl_set_error_handler (&shark_error_handler);
  set_nf(nf_in);
  for (int mu=0; mu<4; mu++){
    p1[mu] = p1_in[mu];
    p2[mu] = p2_in[mu];
    pi[mu] = p3_in[mu];
    pj[mu] = p4_in[mu];
  }
  evaluate_external(false);

  double result, B;

  if (i==3 && j==4){
    // do nothing
  } else if (i==3 && j==1){
    B = pTi*pTi/mi2;
  } else if (i==4 && j==1){
    B = pTj*pTj/mj2;
  } else if (i==3 && j==2){
    B = pTi*pTi/mi2;
  } else if (i==4 && j==2){
    B = pTj*pTj/mj2;
  } else{ // any other case is zero
    return 0;
  }
  
  if (i==3 && j==4){
    result = Fex1_T3T4();
  }
  else { //TITj case
    result = Fex1_TITj(B);
  }

  return result;
}

double H2_TT(int i, int j, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int nf_in){
  gsl_set_error_handler (&shark_error_handler);
  set_nf(nf_in);
  for (int mu=0; mu<4; mu++){
    p1[mu] = p1_in[mu];
    p2[mu] = p2_in[mu];
    pi[mu] = p3_in[mu];
    pj[mu] = p4_in[mu];
  }
  evaluate_external(false);

  double result, A, B;

  if (i==3 && j==4){
    // do nothing
  } else if (i==3 && j==1){
    B = pTi*pTi/mi2;
    A = 2.*pow(p1pi,2)/mi2/p1p2;
  } else if (i==4 && j==1){
    B = pTj*pTj/mj2;
    A = 2.*pow(p1pj,2)/mj2/p1p2;
  } else if (i==3 && j==2){
    B = pTi*pTi/mi2;
    A = 2.*pow(p2pi,2)/mi2/p1p2;
  } else if (i==4 && j==2){
    B = pTj*pTj/mj2;
    A = 2.*pow(p2pj,2)/mj2/p1p2;
  } else{ // any other case is zero
    return 0;
  }
  
  if (i==3 && j==4){
    if (pipj_save1 == pipj && p1pj_save1 == p1pj && p2pj_save1 == p2pj && result_save_T3T4 != -1e10){
      result = result_save_T3T4;
    }
    else{
      pipj_save1 = pipj, p1pj_save1 = p1pj, p2pj_save1 = p2pj;
      result = Fex2_and_11_T3T4();
      result_save_T3T4 = result;
    }
  }
  else { //TITj case
    result = Fex2_TITj(B,A)-2.*Pi*beta0*Fex1_1_TITj(B,A);
  }

  return result;
}

double H2_TTTT(int i, int j, int k, int l, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int nf_in){
  gsl_set_error_handler (&shark_error_handler);
  set_nf(nf_in);
  // contributions coming from <Fex1(0)^2> -1/2 <Fex,1(0)>^2
  for (int mu=0; mu<4; mu++){
    p1[mu] = p1_in[mu];
    p2[mu] = p2_in[mu];
    pi[mu] = p3_in[mu];
    pj[mu] = p4_in[mu];
  }
  evaluate_external(false);
  double result = 0;

  if(i==3 && j==4 && k==3 && l==4){
    result = Fex1sq_T3T4T3T4_phi_avg() - sqFex1_T3T4T3T4()/2.;
  }
  else if(i==3 && j==4 && (k==3 || k==4) && (l==1 || l==2)){
    result = Fex1sq_T3T4TITj_phi_avg(k) - sqFex1_T3T4TITj(k)/2.;
    result = 2*result; // with factor 2 no need to include ij<->kl
  }
  else if((i==3 || i==4) && (j==1 || j==2) && (k==3 || k==4) && (l==1 || l==2)){
    result = Fex1sq_TITjTKTl_phi_avg(i,k) - sqFex1_TITjTKTl(i,k)/2.;
  }

  return result;
}

double H2_fTTT(int i, int j, int k, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int nf_in){
  gsl_set_error_handler (&shark_error_handler);
  set_nf(nf_in);
  for (int mu=0; mu<4; mu++){
    p1[mu] = p1_in[mu];
    p2[mu] = p2_in[mu];
    pi[mu] = p3_in[mu];
    pj[mu] = p4_in[mu];
  }
  
  double result = 0;

  if(i==1 && j==2 && k==3){
    evaluate_external(false);
    if (pipj_save2 == pipj && p1pj_save2 == p1pj && p2pj_save2 == p2pj && result_save_fTTT != -1e10){
      result = result_save_fTTT;
    }
    else{
      pipj_save2 = pipj, p1pj_save2 = p1pj, p2pj_save2 = p2pj;
      result = Fex2_fTTT();
      result_save_fTTT = result;
    }
  }

  return result;
}

}
// end of shark namespace, now defining global functions for the interface:


double H1_TT(int i, int j, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int nf_in){
  return shark::H1_TT(i,j,p1_in,p2_in,p3_in,p4_in,nf_in);
}
double H2_TT(int i, int j, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int nf_in){
  return shark::H2_TT(i,j,p1_in,p2_in,p3_in,p4_in,nf_in);
}
double H2_TTTT(int i, int j, int k, int l, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int nf_in){
  return shark::H2_TTTT(i,j,k,l,p1_in,p2_in,p3_in,p4_in,nf_in);
}
double H2_fTTT(int i, int j, int k, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int nf_in){
  return shark::H2_fTTT(i,j,k,p1_in,p2_in,p3_in,p4_in,nf_in);
}

// Functions with all pointers in arguments
double H1_TT_shark(int *i, int *j, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int *nf_in){
  return shark::H1_TT(i[0],j[0],p1_in,p2_in,p3_in,p4_in,nf_in[0]);
}
double H2_TT_shark(int *i, int *j, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int *nf_in){
  return shark::H2_TT(i[0],j[0],p1_in,p2_in,p3_in,p4_in,nf_in[0]);
}
double H2_TTTT_shark(int *i, int *j, int *k, int *l, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int *nf_in){
  return shark::H2_TTTT(i[0],j[0],k[0],l[0],p1_in,p2_in,p3_in,p4_in,nf_in[0]);
}
double H2_fTTT_shark(int *i, int *j, int *k, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int *nf_in){
  return shark::H2_fTTT(i[0],j[0],k[0],p1_in,p2_in,p3_in,p4_in,nf_in[0]);
}

// Precision settings
void setprecision_T3T4_shark(double *acc){
  shark::acc_goal_T3T4 = acc[0];
}
void setprecision_T3T4(double acc){
  shark::acc_goal_T3T4 = acc;
}
void setprecision_fTTT_shark(double *acc){
  shark::acc_goal_fTTT = acc[0];
}
void setprecision_fTTT(double acc){
  shark::acc_goal_fTTT = acc;
}
void setcalls_T3T4_shark(int *calls){
  shark::calls_T3T4 = calls[0];
}
void setcalls_T3T4(int calls){
  shark::calls_T3T4 = calls;
}
void setprecision_mode(int mode){
  // Some pre-defined settings
  if (mode==1){
    // default setup:
    shark::acc_goal_T3T4 = 5.;
    shark::acc_goal_fTTT = 5.;
    shark::calls_T3T4 = 5000;
    shark::calls_fTTT = 5000;
  }
  else if (mode==2){
    // higher precision
    shark::acc_goal_T3T4 = 1.;
    shark::acc_goal_fTTT = 1.;
    shark::calls_T3T4 = 20000;
    shark::calls_fTTT = 20000;
  }
  else if (mode==3){
    // even higher precision
    shark::acc_goal_T3T4 = 0.01;
    shark::acc_goal_fTTT = 0.01;
    shark::calls_T3T4 = 1000000;
    shark::calls_fTTT = 1000000;
  }
  else if (mode==0){
    // fast mode: only analytical results are evaluated
    shark::acc_goal_T3T4 = -1;
    shark::acc_goal_fTTT = -1;
    shark::calls_T3T4 = 0;
  }
  else {
    std::cout << "shark: invalid value in precision mode" << std::endl;
    exit(EXIT_FAILURE);
  }
}
void setprecision_mode_shark(int *mode){
  setprecision_mode(mode[0]);
}
