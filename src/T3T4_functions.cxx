
void evaluate_int_variables(double u, double v, double cosp, double y){
  double sinp = sqrt(1 - cosp*cosp);
  aux$ni2  = 4*mi2*pow(v,2)*pow(pi$minus + pi$minus*u - 2*cosp*pTi*v + 2*pi$plus*pow(v,2),-2);
  aux$nj2  = 4*mj2*pow(v,2)*pow(pj$minus*(1 + u) +  2*v*(-(cosp*cosxi*pTj) + pj$plus*v + pTj*sinxi*y*sinp),-2);
  aux$ninj = 4*pipj*pow(v,2)*pow(pj$minus*(1 + u) +  2*v*(-(cosp*cosxi*pTj) + pj$plus*v + pTj*sinxi*y*sinp),-1)*pow(pi$minus + pi$minus*u - 2*cosp*pTi*v + 2*pi$plus*pow(v,2),-1);
  ni2  = 1 - aux$ni2*u;
  nj2  = 1 - aux$nj2*u;
  ninj = 1 - aux$ninj*u;
  aux$ni2$u0  = 4*mi2*pow(v,2)*pow(pi$minus - 2*cosp*pTi*v + 2*pi$plus*pow(v,2),-2);
  aux$nj2$u0  = 4*mj2*pow(v,2)*pow(pj$minus +  2*v*(-(cosp*cosxi*pTj) + pj$plus*v + pTj*sinxi*y*sinp),-2);
  aux$ninj$u0 = 4*pipj*pow(v,2)*pow(pj$minus +  2*v*(-(cosp*cosxi*pTj) + pj$plus*v + pTj*sinxi*y*sinp),-1)*pow(pi$minus - 2*cosp*pTi*v + 2*pi$plus*pow(v,2),-1);
}


void evaluate_extra_variables(double u, double v, double cosp){
  //double sinp = sqrt(1 - cosp*cosp);
  Lu = log(u);
  Bij = pow(aux$ni2 - aux$ninj,2)/(aux$ninj*aux$ninj - aux$ni2*aux$nj2);
  Bij$u0 = pow(aux$ni2$u0 - aux$ninj$u0,2)/(aux$ninj$u0*aux$ninj$u0 - aux$ni2$u0*aux$nj2$u0);
  sqrtBij = sqrt(Bij);
  sqrtBij$u0 = sqrt(Bij$u0);
  LBij$minus = log(1-sqrtBij);
  LBij$plus  = log(1+sqrtBij);
  LBij$u0$minus = log(1-sqrtBij$u0);
  LBij$u0$plus  = log(1+sqrtBij$u0);
  Laux$ni2 = log(aux$ni2);
  Lni2 = log(ni2);
  Lcosp = log(1-cosp*cosp);
  Laux$ni2$u0 = log(aux$ni2$u0);
  b = sqrt(1-mi2*mi2/pipj/pipj);

  Li2Bij$plus = Li2((1+sqrtBij)/2);
  Li2Bij$minus = Li2((1-sqrtBij)/2);
  Li2Bij$u0$plus = Li2((1+sqrtBij$u0)/2);
  Li2Bij$u0$minus = Li2((1-sqrtBij$u0)/2);
    
}


void evaluate_extra_variables_delta(double v, double cosp, double y){
  Laux1 = log(pipj*pow(-(cosp*pTi) + pi$plus*v + (pi$minus*pow(v,-1))/2.,-1)*pow(pj$plus*v + (pj$minus*pow(v,-1))/2. - pTj*(cosp*cosxi - sinxi*y*pow(1 - pow(cosp,2),0.5)),-1));
  Laux2 = log((1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5))/2.);
  Laux3 = log(pow(mi2,2)*pow(pipj,-2)*pow(1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5),-2));
  Laux4 = log(mi2*pow(pipj,-1)*(-(cosp*pTi) + pi$plus*v + (pi$minus*pow(v,-1))/2.)*pow(pj$plus*v + (pj$minus*pow(v,-1))/2. - pTj*(cosp*cosxi - sinxi*y*pow(1 - pow(cosp,2),0.5)),-1)*pow(1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5),-1));
  Laux5 = log(mi2*pow(pipj,-1)*(pj$plus*v + (pj$minus*pow(v,-1))/2. - pTj*(cosp*cosxi - sinxi*y*pow(1 - pow(cosp,2),0.5)))*pow(-(cosp*pTi) + pi$plus*v + (pi$minus*pow(v,-1))/2.,-1)*pow(1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5),-1));
  Laux6 = log(pow(1 - pow(mi2,2)*pow(pipj,-2),0.5));
  Laux7 = log(1 - mi2*pow(pipj,-1)*(-(cosp*pTi) + pi$plus*v + (pi$minus*pow(v,-1))/2.)*pow(pj$plus*v + (pj$minus*pow(v,-1))/2. - pTj*(cosp*cosxi - sinxi*y*pow(1 - pow(cosp,2),0.5)),-1)*pow(1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5),-1));
  Laux8 = log(1 - mi2*pow(pipj,-1)*(pj$plus*v + (pj$minus*pow(v,-1))/2. - pTj*(cosp*cosxi - sinxi*y*pow(1 - pow(cosp,2),0.5)))*pow(-(cosp*pTi) + pi$plus*v + (pi$minus*pow(v,-1))/2.,-1)*pow(1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5),-1));

  Paux1 = Li2((1 - pow(1 - pow(mi2,2)*pow(pipj,-2),0.5))*pow(1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5),-1));
  Paux2 = Li3(1 - (1 - pow(1 - pow(mi2,2)*pow(pipj,-2),0.5))*pow(1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5),-1));
  Paux3 = Li3((1 - pow(1 - pow(mi2,2)*pow(pipj,-2),0.5))*pow(1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5),-1));
  Paux4 = Li2(mi2*pow(pipj,-1)*(-(cosp*pTi) + pi$plus*v + (pi$minus*pow(v,-1))/2.)*pow(pj$plus*v + (pj$minus*pow(v,-1))/2. - pTj*(cosp*cosxi - sinxi*y*pow(1 - pow(cosp,2),0.5)),-1)*pow(1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5),-1));
  Paux5 = Li2(mi2*pow(pipj,-1)*(pj$plus*v + (pj$minus*pow(v,-1))/2. - pTj*(cosp*cosxi - sinxi*y*pow(1 - pow(cosp,2),0.5)))*pow(-(cosp*pTi) + pi$plus*v + (pi$minus*pow(v,-1))/2.,-1)*pow(1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5),-1));
  Paux6 = Li3(mi2*pow(pipj,-1)*(-(cosp*pTi) + pi$plus*v + (pi$minus*pow(v,-1))/2.)*pow(pj$plus*v + (pj$minus*pow(v,-1))/2. - pTj*(cosp*cosxi - sinxi*y*pow(1 - pow(cosp,2),0.5)),-1)*pow(1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5),-1));
  Paux7 = Li3(mi2*pow(pipj,-1)*(pj$plus*v + (pj$minus*pow(v,-1))/2. - pTj*(cosp*cosxi - sinxi*y*pow(1 - pow(cosp,2),0.5)))*pow(-(cosp*pTi) + pi$plus*v + (pi$minus*pow(v,-1))/2.,-1)*pow(1 + pow(1 - pow(mi2,2)*pow(pipj,-2),0.5),-1));

  Li3Bij$u0$plus = Li3((1+sqrtBij$u0)/2);
  Li3Bij$u0$minus = Li3((1-sqrtBij$u0)/2);
}

void evaluate_extra_variables_reg(){
  alpha = -((mypow(ninj,2) - ni2*nj2)*mypow(ni2 - 2*ninj + mypow(ninj,2) + nj2 - ni2*nj2,-1));
  mu = abs(ni2 + ninj)*mypow(sqrt(ni2 + ninj*(2 + ninj) + nj2 - ni2*nj2),-1);
  csi = sqrt(mypow(ni2 - ninj,2)*mypow(ni2 - 2*ninj + mypow(ninj,2) + nj2 - ni2*nj2,-1));
  omega = (sqrt(ni2)*(ni2 - 2*ninj + nj2))*mypow(ni2 - ninj,-2);
  rho = (-mypow(ninj,2) + ni2*nj2)*mypow(ni2 - ninj,-2);
  gammaij = (mypow(ninj,2) - ni2*nj2)*mypow(ni2 + 2*ninj + mypow(ninj,2) + nj2 - ni2*nj2,-1);
  Dij = ((1 + ninj)*(ni2 + ninj))*mypow(ni2 + 2*ninj + mypow(ninj,2) + nj2 - ni2*nj2,-1);
  sqrt$ni2 = sqrt(ni2);
  sqrt$aux = sqrt(ni2+gammaij-ni2*gammaij);
}


double DIST2t(double u, double v, double cosp){
  double sinp = sqrt(1 - cosp*cosp);
  double res = (pow(aux$ni2 - aux$ninj,-1)*(3*aux$ninj*(aux$ni2 + aux$nj2)*(2*Li2Bij$minus - 2*Li2Bij$plus + (LBij$minus - LBij$plus)*(LBij$minus + LBij$plus - 2*Log2))*sqrtBij + 
       12*aux$ninj*(aux$ni2 + aux$nj2)*(LBij$minus - LBij$plus)*(Log2 - Lu)*sqrtBij + 
       6*(Laux$ni2 + Lcosp + 2*Log2)*((aux$ni2 - aux$ninj)*(aux$ni2 - aux$nj2)*(Laux$ni2 - 2*Log2 + Lu) + aux$ninj*(aux$ni2 + aux$nj2)*(-LBij$minus + LBij$plus)*sqrtBij) - 
       3*(aux$ni2 - aux$ninj)*(aux$ni2 - aux$nj2)*pow(Laux$ni2 - Lu,2) + 12*(aux$ni2 - aux$ninj)*(aux$ni2 - aux$nj2)*pow(Log2 - Lu,2) + 
       2*(aux$ni2 - aux$ninj)*(aux$ni2 - aux$nj2)*pow(Pi,2) + 2*aux$ninj*
        (6*aux$ni2*Laux$ni2*Lcosp - 6*aux$ninj*Laux$ni2*Lcosp - 12*aux$ni2*Lcosp*Log2 + 12*aux$ninj*Lcosp*Log2 + 6*aux$ninj*Laux$ni2*(LBij$minus - LBij$plus)*sqrtBij + 
          6*aux$ninj*(LBij$minus - LBij$plus)*Lcosp*sqrtBij - 6*aux$ninj*Li2Bij$minus*sqrtBij + 6*aux$ninj*Li2Bij$plus*sqrtBij + 6*aux$ninj*LBij$minus*Log2*sqrtBij - 
          6*aux$ninj*LBij$plus*Log2*sqrtBij + 6*Lu*((aux$ni2 - aux$ninj)*(2*Laux$ni2 + Lcosp - 2*Log2) + 2*aux$ninj*(LBij$minus - LBij$plus)*sqrtBij) + 
          3*aux$ni2*pow(Laux$ni2,2) - 3*aux$ninj*pow(Laux$ni2,2) - 3*aux$ninj*sqrtBij*pow(LBij$minus,2) + 3*aux$ninj*sqrtBij*pow(LBij$plus,2) - 12*aux$ni2*pow(Log2,2) + 
	 12*aux$ninj*pow(Log2,2) + 9*(aux$ni2 - aux$ninj)*pow(Lu,2) + 2*aux$ni2*pow(Pi,2) - 2*aux$ninj*pow(Pi,2)))*pow(sinp,-1)*pow(u,-1)*pow(v,-1))/6.;
  return res;
}


double DISTt(double u, double v, double cosp){
  double sinp = sqrt(1 - cosp*cosp);  
  double res = ((3*aux$ninj*(aux$ni2 + aux$nj2)*(2*Li2Bij$minus - 2*Li2Bij$plus + (LBij$minus - LBij$plus)*(LBij$minus + LBij$plus - 2*Log2))*sqrtBij*pow(aux$ni2 - aux$ninj,-1) + 
       12*aux$ninj*(aux$ni2 + aux$nj2)*(LBij$minus - LBij$plus)*(Log2 - Lu)*sqrtBij*pow(aux$ni2 - aux$ninj,-1) + 
       6*(Laux$ni2 + Lcosp + 2*Log2)*((aux$ni2 - aux$ninj)*(aux$ni2 - aux$nj2)*(Laux$ni2 - 2*Log2 + Lu) + aux$ninj*(aux$ni2 + aux$nj2)*(-LBij$minus + LBij$plus)*sqrtBij)*
        pow(aux$ni2 - aux$ninj,-1) - 3*aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*
        (2*Li2Bij$u0$minus - 2*Li2Bij$u0$plus + (LBij$u0$minus - LBij$u0$plus)*(LBij$u0$minus + LBij$u0$plus - 2*Log2))*sqrtBij$u0*pow(aux$ni2$u0 - aux$ninj$u0,-1) + 
       12*aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*(-LBij$u0$minus + LBij$u0$plus)*(Log2 - Lu)*sqrtBij$u0*pow(aux$ni2$u0 - aux$ninj$u0,-1) - 
       6*(Laux$ni2$u0 + Lcosp + 2*Log2)*((aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*(Laux$ni2$u0 - 2*Log2 + Lu) + 
          aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*(-LBij$u0$minus + LBij$u0$plus)*sqrtBij$u0)*pow(aux$ni2$u0 - aux$ninj$u0,-1) - 3*(aux$ni2 - aux$nj2)*pow(Laux$ni2 - Lu,2) + 
       3*(aux$ni2$u0 - aux$nj2$u0)*pow(Laux$ni2$u0 - Lu,2) + 12*(aux$ni2 - aux$nj2)*pow(Log2 - Lu,2) - 12*(aux$ni2$u0 - aux$nj2$u0)*pow(Log2 - Lu,2) + 
       2*(aux$ni2 - aux$nj2)*pow(Pi,2) - 2*(aux$ni2$u0 - aux$nj2$u0)*pow(Pi,2) + 
       2*aux$ninj*pow(aux$ni2 - aux$ninj,-1)*(aux$ni2*(-12*Log2*Lu + 6*Lcosp*(-2*Log2 + Lu) + 6*Laux$ni2*(Lcosp + 2*Lu) + 3*pow(Laux$ni2,2) - 12*pow(Log2,2) + 9*pow(Lu,2) + 
             2*pow(Pi,2)) - aux$ninj*(-12*Lcosp*Log2 + 6*Lcosp*Lu - 12*Log2*Lu - 6*LBij$minus*Lcosp*sqrtBij + 6*LBij$plus*Lcosp*sqrtBij + 6*Li2Bij$minus*sqrtBij - 
             6*Li2Bij$plus*sqrtBij - 6*LBij$minus*Log2*sqrtBij + 6*LBij$plus*Log2*sqrtBij - 12*LBij$minus*Lu*sqrtBij + 12*LBij$plus*Lu*sqrtBij + 
             6*Laux$ni2*(Lcosp + 2*Lu + (-LBij$minus + LBij$plus)*sqrtBij) + 3*pow(Laux$ni2,2) + 3*sqrtBij*pow(LBij$minus,2) - 3*sqrtBij*pow(LBij$plus,2) - 12*pow(Log2,2) + 
             9*pow(Lu,2) + 2*pow(Pi,2))) + 2*aux$ninj$u0*pow(aux$ni2$u0 - aux$ninj$u0,-1)*
        (-(aux$ni2$u0*(-12*Log2*Lu + 6*Lcosp*(-2*Log2 + Lu) + 6*Laux$ni2$u0*(Lcosp + 2*Lu) + 3*pow(Laux$ni2$u0,2) - 12*pow(Log2,2) + 9*pow(Lu,2) + 2*pow(Pi,2))) + 
          aux$ninj$u0*(-12*Lcosp*Log2 + 6*Lcosp*Lu - 12*Log2*Lu - 6*LBij$u0$minus*Lcosp*sqrtBij$u0 + 6*LBij$u0$plus*Lcosp*sqrtBij$u0 + 6*Li2Bij$u0$minus*sqrtBij$u0 - 
             6*Li2Bij$u0$plus*sqrtBij$u0 - 6*LBij$u0$minus*Log2*sqrtBij$u0 + 6*LBij$u0$plus*Log2*sqrtBij$u0 - 12*LBij$u0$minus*Lu*sqrtBij$u0 + 
             12*LBij$u0$plus*Lu*sqrtBij$u0 + 6*Laux$ni2$u0*(Lcosp + 2*Lu + (-LBij$u0$minus + LBij$u0$plus)*sqrtBij$u0) + 3*pow(Laux$ni2$u0,2) + 
		       3*sqrtBij$u0*pow(LBij$u0$minus,2) - 3*sqrtBij$u0*pow(LBij$u0$plus,2) - 12*pow(Log2,2) + 9*pow(Lu,2) + 2*pow(Pi,2))))*pow(sinp,-1)*pow(u,-1)*pow(v,-1))/6.;
  return res;
}

double DELTAt(double v, double cosp, double y){
  if(abs(mi2-mj2) > 0.001){
    std::cout << "mi different from mj. STOP" << std::endl;
    exit(1);
  }
  double sinp = sqrt(1 - cosp*cosp);
  double res = (pow(aux$ni2$u0 - aux$ninj$u0,-1)*(-24*aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*
         (2*Li2Bij$u0$minus - 2*Li2Bij$u0$plus + (LBij$u0$minus - LBij$u0$plus)*(LBij$u0$minus + LBij$u0$plus - 2*Log2))*Log2*sqrtBij$u0 - 
        48*(aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*z3 + 8*(aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*pow(Laux$ni2$u0,3) + 
        48*aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*(-LBij$u0$minus + LBij$u0$plus)*sqrtBij$u0*pow(Log2,2) - 
        12*(Laux$ni2$u0 + Lcosp + 2*Log2)*(-(aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*
              (2*Li2Bij$u0$minus - 2*Li2Bij$u0$plus + (LBij$u0$minus - LBij$u0$plus)*(LBij$u0$minus + LBij$u0$plus - 2*Log2))*sqrtBij$u0) + 
           4*aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*(-LBij$u0$minus + LBij$u0$plus)*Log2*sqrtBij$u0 + 
           2*(aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*pow(Laux$ni2$u0,2) - 4*(aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*pow(Log2,2)) - 
        32*(aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*pow(Log2,3) + 
        8*((aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*(Laux$ni2$u0 + Lcosp) + aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*(-LBij$u0$minus + LBij$u0$plus)*sqrtBij$u0)*
         pow(Pi,2) - 4*(2*(aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*(Laux$ni2$u0 - Log2) + 
           aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*(-LBij$u0$minus + LBij$u0$plus)*sqrtBij$u0)*
         (-12*Lcosp*Log2 - 6*Laux$ni2$u0*(Lcosp + 2*Log2) - 3*pow(Laux$ni2$u0,2) - 3*pow(Lcosp,2) - 12*pow(Log2,2) + pow(Pi,2)) - 
        2*aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*sqrtBij$u0*(-18*LBij$u0$minus*Li2Bij$u0$minus + 12*LBij$u0$plus*Li2Bij$u0$minus - 18*LBij$u0$minus*Li2Bij$u0$plus + 
           12*LBij$u0$plus*Li2Bij$u0$plus + 24*Li3Bij$u0$minus - 24*Li3Bij$u0$plus + 6*LBij$u0$minus*LBij$u0$plus*Log2 - 12*LBij$u0$plus*pow(LBij$u0$minus,2) + 
           6*Log2*pow(LBij$u0$minus,2) + 2*pow(LBij$u0$minus,3) + 6*LBij$u0$minus*pow(LBij$u0$plus,2) - 2*pow(LBij$u0$plus,3) - 6*LBij$u0$minus*pow(Log2,2) + 
           LBij$u0$minus*pow(Pi,2)) - 4*(aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*
         (12*z3 + 3*(Lcosp + 2*Log2)*pow(Laux$ni2$u0,2) + pow(Laux$ni2$u0,3) + 6*Log2*pow(Lcosp,2) + pow(Lcosp,3) + 12*Lcosp*pow(Log2,2) + 8*pow(Log2,3) + 
           Laux$ni2$u0*(12*Lcosp*Log2 + 3*pow(Lcosp,2) + 12*pow(Log2,2) - pow(Pi,2)) - Lcosp*pow(Pi,2) - 2*Log2*pow(Pi,2)))*pow(sinp,-1)*pow(v,-1))/48. + 
   (aux$ninj$u0*pow(aux$ni2$u0 - aux$ninj$u0,-1)*(aux$ninj$u0*(-18*LBij$u0$minus*Li2Bij$u0$minus*sqrtBij$u0 + 12*LBij$u0$plus*Li2Bij$u0$minus*sqrtBij$u0 - 
           18*LBij$u0$minus*Li2Bij$u0$plus*sqrtBij$u0 + 12*LBij$u0$plus*Li2Bij$u0$plus*sqrtBij$u0 + 24*Li3Bij$u0$minus*sqrtBij$u0 - 24*Li3Bij$u0$plus*sqrtBij$u0 + 
           6*LBij$u0$minus*LBij$u0$plus*Log2*sqrtBij$u0 + 48*z3 - 6*(Lcosp + (-LBij$u0$minus + LBij$u0$plus)*sqrtBij$u0)*pow(Laux$ni2$u0,2) - 2*pow(Laux$ni2$u0,3) - 
           12*LBij$u0$plus*sqrtBij$u0*pow(LBij$u0$minus,2) + 6*Log2*sqrtBij$u0*pow(LBij$u0$minus,2) + 2*sqrtBij$u0*pow(LBij$u0$minus,3) + 
           6*LBij$u0$minus*sqrtBij$u0*pow(LBij$u0$plus,2) - 2*sqrtBij$u0*pow(LBij$u0$plus,3) + 6*(4*Log2 + (LBij$u0$minus - LBij$u0$plus)*sqrtBij$u0)*pow(Lcosp,2) + 
           2*pow(Lcosp,3) - 6*LBij$u0$minus*sqrtBij$u0*pow(Log2,2) + 32*pow(Log2,3) - 8*Log2*pow(Pi,2) + 3*LBij$u0$minus*sqrtBij$u0*pow(Pi,2) - 
           2*LBij$u0$plus*sqrtBij$u0*pow(Pi,2) - 2*Laux$ni2$u0*(6*(-LBij$u0$minus + LBij$u0$plus)*Lcosp*sqrtBij$u0 + 
              3*sqrtBij$u0*(2*Li2Bij$u0$minus - 2*Li2Bij$u0$plus - 2*LBij$u0$minus*Log2 + 2*LBij$u0$plus*Log2 + pow(LBij$u0$minus,2) - pow(LBij$u0$plus,2)) + 
              3*pow(Lcosp,2) + pow(Pi,2)) - 6*Lcosp*(2*(-LBij$u0$minus + LBij$u0$plus)*Log2*sqrtBij$u0 + 
              sqrtBij$u0*(2*Li2Bij$u0$minus - 2*Li2Bij$u0$plus + pow(LBij$u0$minus,2) - pow(LBij$u0$plus,2)) - 8*pow(Log2,2) + pow(Pi,2))) + 
        2*aux$ni2$u0*(3*Lcosp*pow(Laux$ni2$u0,2) + pow(Laux$ni2$u0,3) - 12*Log2*pow(Lcosp,2) - pow(Lcosp,3) + Laux$ni2$u0*(3*pow(Lcosp,2) + pow(Pi,2)) + 
           3*Lcosp*(-8*pow(Log2,2) + pow(Pi,2)) + 4*(-6*z3 - 4*pow(Log2,3) + Log2*pow(Pi,2))))*pow(sinp,-1)*pow(v,-1))/2/6 - 
   Pi*pow(v,-1)*(2 - mi2*(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v))*pow(pipj,-1)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1) - 
      mi2*(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v))*pow(pipj,-1)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1))*
    (-(pipj*pow(Laux1 - Lcosp - Log2,3)*pow(Pi,-1)*pow(sinp,-1)*pow(v,2)*pow(-pi$minus + 2*v*(cosp*pTi - pi$plus*v),-1)*
          pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1))/3. + 
      2*Log2*pipj*pow(-Laux1 + Lcosp + Log2,2)*pow(Pi,-1)*pow(sinp,-1)*pow(v,2)*pow(-pi$minus + 2*v*(cosp*pTi - pi$plus*v),-1)*
       pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1) + 
      ((-Laux1 + Lcosp + Log2)*pipj*pow(Pi,-1)*(12*pow(Log2,2) - pow(Pi,2))*pow(sinp,-1)*pow(v,2)*pow(-pi$minus + 2*v*(cosp*pTi - pi$plus*v),-1)*
         pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1))/3. + 
      ((2*b*Laux2 + (-1 + b)*Laux3)*pipj*pow(b,-1)*pow(Pi,-1)*(-18*Lcosp*Log2 + 6*Laux1*(Lcosp + 3*Log2) - 3*pow(Laux1,2) - 3*pow(Lcosp,2) - 27*pow(Log2,2) + pow(Pi,2))*
         pow(sinp,-1)*pow(v,2)*pow(-pi$minus + 2*v*(cosp*pTi - pi$plus*v),-1)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1))/6. + 
      (2*pipj*pow(Pi,-1)*(6*z3 + 4*pow(Log2,3) - Log2*pow(Pi,2))*pow(sinp,-1)*pow(v,2)*pow(-pi$minus + 2*v*(cosp*pTi - pi$plus*v),-1)*
         pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1))/3. + 
      2*(-Laux1 + Lcosp + 3*Log2)*pipj*pow(Pi,-1)*pow(sinp,-1)*pow(v,2)*pow(-pi$minus + 2*v*(cosp*pTi - pi$plus*v),-1)*
       pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1)*
       ((((1 + b)*Laux2*Laux3 - 2*(Laux3*Laux6 + Paux1))*pow(b,-1))/2. + pow(Laux2,2)/2. + ((-9.5 + 7*pow(b,-1))*pow(Pi,2))/6. + 
         (pow(b,-1)*pow(pipj,-1)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1)*
            (mi2*pow(Laux4,2)*((1 + b)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),2) + (-1 + b)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),2)) - 
              mi2*pow(Laux5,2)*((1 - b)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),2) - (1 + b)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),2)))*
            pow(2 - mi2*(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v))*pow(pipj,-1)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1) - 
              mi2*(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v))*pow(pipj,-1)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1),-1))/4.) + 
      2*pipj*pow(Pi,-1)*pow(sinp,-1)*pow(v,2)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1)*
       (-(z3*(2.3333333333333335 + pow(b,-1))) + pow(Laux2,3)/6. + pow(b,-1)*
          (2*Paux2 + Paux3 + ((-1 + b)*(pow(Laux4,3) + pow(Laux5,3)))/2/6 + 2*Laux6*pow(Pi,2) - Laux2*(Paux1 + ((10 + 19*b)*pow(Pi,2))/12.) + 
            (-8*Laux2*Laux3*Laux6 + (6 + 2*b)*Laux3*pow(Laux2,2) + 3*(-1 + b)*Laux3*pow(Pi,2))/8.) + 
         ((-7*Laux3*pow(Pi,2))/6. + (-(Laux3*Paux1) + 2*Paux3 + 2*Laux4*Paux4 + 2*Laux5*Paux5 - 2*Paux6 - 2*Paux7 + 2*z3 + Laux7*pow(Laux4,2) + Laux8*pow(Laux5,2))*
             (-1 + mi2*(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v))*pow(pipj,-1)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1) + 
               mi2*(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v))*pow(pipj,-1)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1)) + 
            (pow(b,-1)*pow(pipj,-1)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1)*
               (-4*(Laux5*(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v))*
                     (pipj*(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v)) - mi2*(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v))) - 
                    Laux4*(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v))*(mi2*(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v)) - 
                       pipj*(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v))))*(Paux1 - (7*pow(Pi,2))/6.) + 
                 2*Laux6*mi2*(pow(Laux4,2) - pow(Laux5,2))*(pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),2) - 
                    pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),2)) + 
                 Laux2*(mi2*pow(Laux5,2)*((1 + b)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),2) + 
                       (-1 + b)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),2)) - 
                    mi2*pow(Laux4,2)*((1 - b)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),2) - 
                       (1 + b)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),2)))))/4.)*
          pow(2 - mi2*(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v))*pow(pipj,-1)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1) - 
	      mi2*(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v))*pow(pipj,-1)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1),-1)));
  // Adding corrections due to the fact that what I had was R$old
  res = res + (v*pow(b,-1)*(3*(-1 + b)*Laux3*pow(Pi,2) + 2*b*(28*z3 + 3*Laux1*pow(Pi,2) + 3*Laux2*pow(Pi,2) - 3*Lcosp*pow(Pi,2) - 9*Log2*pow(Pi,2)))*pow(sinp,-1)*
     (-2*pipj*(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v))*(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v)) + 
       mi2*(4*pi$minus*v*(-(cosp*pTi) + pi$plus*v) + 4*pj$minus*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v) + pow(pi$minus,2) + pow(pj$minus,2) + 
          4*pow(v,2)*(2*pj$plus*pTj*sinp*sinxi*y*v - 2*cosp*(pi$plus*pTi*v + cosxi*pTj*(pTj*sinp*sinxi*y + pj$plus*v)) + 
             pow(cosp,2)*(pow(pTi,2) + pow(pTj,2)*(pow(cosxi,2) - pow(sinxi*y,2))) + pow(pTj,2)*pow(sinxi*y,2) + (pow(pi$plus,2) + pow(pj$plus,2))*pow(v,2))))*
	       pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-2)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-2))/12.;
  return res;
}


double AUX2(double u, double v, double cosp){
  // kept only the r0_epsilon and R0_epsilon pieces
  double sinp = sqrt(1 - cosp*cosp);
  std::complex<double> ii(0.,1.);
  std::complex<double> res;
  res = mypow(sinp,-1)*mypow(u,-2)*mypow(v,-1)*(-6*Lni2*(-1 + ni2) - 2*Lni2*(-1 + ninj) +
      (mypow(ni2,-1.5)*(2*sqrt$ni2*((-1 + ni2)*(ni2 - ninj) + ni2*(Lni2*(-6 + 7*ni2 + nj2) + 2*ni2*log(-4*mypow(-1 + ni2,-1)))) +
      2*(-1 + ni2)*(5*ni2 - ninj)*log(2 - 2*mypow(1 + sqrt$ni2,-1)) -
      log(-4*ni2*mypow(-1 + ni2,-1))*(ninj - ni2*(5 + ninj) + 2*nj2*mypow(ni2,1.5) + 5*mypow(ni2,2) + 2*mypow(ni2,2.5))))/2. +
      Dij*(-2 + ni2 + nj2)*(-log(1 - sqrt$aux) + log(1 + sqrt$aux))*mypow(sqrt$aux,-1) +
      (-1 + ninj)*(4*Log2 - 2*log(-1 + mypow(ni2,-1)) + (1 + Dij)*(log(-1 + mypow(ni2,-1)) + 2*log(1 - mypow(1 + sqrt$ni2,-1)))*mypow(sqrt$ni2,-1) -
      Dij*gammaij*(sqrt$ni2*(log(1 - mu) - log(1 + mu)) - mu*log(1 - ni2) + 2*mu*log(1 + sqrt$ni2))*mypow(-1 + gammaij,-1)*mypow(mu,-1)*(-1 + mypow(mu,2))*mypow(sqrt$ni2,-1)*
      mypow(ni2 - mypow(mu,2),-1)));
  res = - res * (Lcosp + 2*Log2 + Lu);
  
  if (isnan(res.real()) || isinf(res.real())) {
    if(u<0.0001){ //if u is small and original result is NaN or Inf, try using approximate expression
      res = ((-((aux$ninj - aux$nj2)*(2*aux$ninj - aux$nj2)) + aux$ni2*aux$ninj*(-5 + log(256)) - 2*aux$ni2*(aux$ni2 + 2*aux$ninj)*log(aux$ni2*u) + 
       (1 + log(16))*mypow(aux$ni2,2))*mypow(sinp,-1)*mypow(v,-1))/4.;
      res = - res * (Lcosp + 2*Log2 + Lu);
    }
    else{
      return 0.;
    }
  }
  return res.real();
}


double DIST2t_0(double u, double v, double cosp, double y){
  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  double sinp = sqrt(1 - cosp*cosp);
  double res = (-((aux$ni2 - aux$ninj)*(aux$ni2 + 2*aux$ninj - aux$nj2)*(Laux$ni2 - 2*Log2 + Lu))
	     - aux$ninj*(aux$ni2 - 2*aux$ninj + aux$nj2)*(-LBij$minus + LBij$plus)*sqrtBij)*
               pow(aux$ni2 - aux$ninj,-1)*pow(sinp,-1)*pow(u,-1)*pow(v,-1);
  return res;
}


double DISTt_0(double u, double v, double cosp, double y){
  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  double sinp = sqrt(1 - cosp*cosp);
  double res = ((-aux$ni2 - 2*aux$ninj + aux$nj2)*Laux$ni2 + aux$ni2$u0*Laux$ni2$u0 + 2*aux$ninj$u0*Laux$ni2$u0 - aux$nj2$u0*Laux$ni2$u0 + 2*aux$ni2*Log2 - 2*aux$ni2$u0*Log2 +
	    4*aux$ninj*Log2 - 4*aux$ninj$u0*Log2 - 2*aux$nj2*Log2 + 2*aux$nj2$u0*Log2 - aux$ni2*Lu + aux$ni2$u0*Lu - 2*aux$ninj*Lu + 2*aux$ninj$u0*Lu + aux$nj2*Lu -
	    aux$nj2$u0*Lu + aux$ninj*(aux$ni2 - 2*aux$ninj + aux$nj2)*(LBij$minus - LBij$plus)*sqrtBij*pow(aux$ni2 - aux$ninj,-1) +
	    aux$ninj$u0*(aux$ni2$u0 - 2*aux$ninj$u0 + aux$nj2$u0)*(-LBij$u0$minus + LBij$u0$plus)*sqrtBij$u0*pow(aux$ni2$u0 - aux$ninj$u0,-1))*pow(sinp,-1)*pow(u,-1)*pow(v,-1);
  return res;
}


double AUX2_0(double u, double v, double cosp, double y){
  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  evaluate_extra_variables_reg();
  double sinp = sqrt(1 - cosp*cosp);
  double res = mypow(sinp,-1)*mypow(u,-2)*mypow(v,-1)*(-6*Lni2*(-1 + ni2) - 2*Lni2*(-1 + ninj) +
      (mypow(ni2,-1.5)*(2*sqrt$ni2*((-1 + ni2)*(ni2 - ninj) + ni2*(Lni2*(-6 + 7*ni2 + nj2) + 2*ni2*log(-4*mypow(-1 + ni2,-1)))) +
      2*(-1 + ni2)*(5*ni2 - ninj)*log(2 - 2*mypow(1 + sqrt$ni2,-1)) -
      log(-4*ni2*mypow(-1 + ni2,-1))*(ninj - ni2*(5 + ninj) + 2*nj2*mypow(ni2,1.5) + 5*mypow(ni2,2) + 2*mypow(ni2,2.5))))/2. +
      Dij*(-2 + ni2 + nj2)*(-log(1 - sqrt$aux) + log(1 + sqrt$aux))*mypow(sqrt$aux,-1) +
      (-1 + ninj)*(4*Log2 - 2*log(-1 + mypow(ni2,-1)) + (1 + Dij)*(log(-1 + mypow(ni2,-1)) + 2*log(1 - mypow(1 + sqrt$ni2,-1)))*mypow(sqrt$ni2,-1) -
      Dij*gammaij*(sqrt$ni2*(log(1 - mu) - log(1 + mu)) - mu*log(1 - ni2) + 2*mu*log(1 + sqrt$ni2))*mypow(-1 + gammaij,-1)*mypow(mu,-1)*(-1 + mypow(mu,2))*mypow(sqrt$ni2,-1)*
      mypow(ni2 - mypow(mu,2),-1)));
  return res;
}

double DELTAt_m2(double v, double cosp, double y){
  if(abs(mi2-mj2) > 0.001){
    std::cout << "mi different from mj. STOP" << std::endl;
    exit(1);
  }
  double u=0.5; //dummy
  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  evaluate_extra_variables_delta(v,cosp,y);
  double sinp = sqrt(1 - cosp*cosp);
  double res = (pow(sinp,-1)*pow(v,-1)*(aux$ni2$u0 + 2*aux$ninj$u0 - aux$nj2$u0 + 4*mi2*pow(v,2)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-2) +
	     4*mi2*pow(v,2)*pow(pi$minus - 2*cosp*pTi*v + 2*pi$plus*pow(v,2),-2) -
	     8*pipj*pow(v,2)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1)*pow(pi$minus - 2*cosp*pTi*v + 2*pi$plus*pow(v,2),-1)))/2.;
  return res;
}

double DELTAt_m1(double v, double cosp, double y){
  if(abs(mi2-mj2) > 0.001){
    std::cout << "mi different from mj. STOP" << std::endl;
    exit(1);
  }
  double u=0.5; //dummy
  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  evaluate_extra_variables_delta(v,cosp,y);
  double sinp = sqrt(1 - cosp*cosp);
  double res = ((aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*(-LBij$u0$minus + LBij$u0$plus)*sqrtBij$u0 +
	 (aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*(Laux$ni2$u0 - log(-16*(-1 + pow(cosp,2)))))*pow(aux$ni2$u0 - aux$ninj$u0,-1)*pow(sinp,-1)*pow(v,-1))/2. +
    aux$ninj$u0*((-aux$ni2$u0 + aux$ninj$u0)*Laux$ni2$u0 + aux$ni2$u0*Lcosp - aux$ninj$u0*Lcosp + 4*aux$ni2$u0*Log2 - 4*aux$ninj$u0*Log2 +
	 aux$ninj$u0*(-LBij$u0$minus + LBij$u0$plus)*sqrtBij$u0)*pow(-aux$ni2$u0 + aux$ninj$u0,-1)*pow(sinp,-1)*pow(v,-1) +
    ((-1 + b)*Laux3 + 2*b*(Laux1 + Laux2 - Lcosp - 3*Log2))*v*pow(b,-1)*pow(sinp,-1)*
    (-2*pipj*(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v))*(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v)) +
     mi2*(4*pi$minus*v*(-(cosp*pTi) + pi$plus*v) + 4*pj$minus*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v) + pow(pi$minus,2) + pow(pj$minus,2) +
	  4*pow(v,2)*(2*pj$plus*pTj*sinp*sinxi*y*v - 2*cosp*(pi$plus*pTi*v + cosxi*pTj*(pTj*sinp*sinxi*y + pj$plus*v)) +
	  pow(cosp,2)*(pow(pTi,2) + pow(pTj,2)*(pow(cosxi,2) - pow(sinxi*y,2))) + pow(pTj,2)*pow(sinxi*y,2) + (pow(pi$plus,2) + pow(pj$plus,2))*pow(v,2))))*
    pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-2)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-2);
  return res;
}

double DELTAt_0(double v, double cosp, double y){
  if(abs(mi2-mj2) > 0.001){
    std::cout << "mi different from mj. STOP" << std::endl;
    exit(1);
  }
  double u=0.5; //dummy
  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  evaluate_extra_variables_delta(v,cosp,y);
  double sinp = sqrt(1 - cosp*cosp);
  double res = (pow(aux$ni2$u0 - aux$ninj$u0,-1)*(-3*aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*
         (2*Li2Bij$u0$minus - 2*Li2Bij$u0$plus + (LBij$u0$minus - LBij$u0$plus)*(LBij$u0$minus + LBij$u0$plus - 2*Log2))*sqrtBij$u0 + 
        12*aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*(-LBij$u0$minus + LBij$u0$plus)*Log2*sqrtBij$u0 - 
        12*(Laux$ni2$u0 + Lcosp + 2*Log2)*((aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*(Laux$ni2$u0 - Log2) + 
           (aux$ninj$u0*(aux$ni2$u0 + aux$nj2$u0)*(-LBij$u0$minus + LBij$u0$plus)*sqrtBij$u0)/2.) + 
        6*(aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*pow(Laux$ni2$u0,2) - 12*(aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*pow(Log2,2) + 
        (aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*(12*Lcosp*Log2 + 6*Laux$ni2$u0*(Lcosp + 2*Log2) + 3*pow(Laux$ni2$u0,2) + 3*pow(Lcosp,2) + 12*pow(Log2,2) - 
           pow(Pi,2)) - 2*(aux$ni2$u0 - aux$ninj$u0)*(aux$ni2$u0 - aux$nj2$u0)*pow(Pi,2))*pow(sinp,-1)*pow(v,-1))/12. + 
   (aux$ninj$u0*pow(-aux$ni2$u0 + aux$ninj$u0,-1)*(aux$ni2$u0*(2*Laux$ni2$u0*Lcosp - 8*Lcosp*Log2 + pow(Laux$ni2$u0,2) - pow(Lcosp,2) - 8*pow(Log2,2) + pow(Pi,2)) - 
        aux$ninj$u0*(2*Li2Bij$u0$minus*sqrtBij$u0 - 2*Li2Bij$u0$plus*sqrtBij$u0 - 2*LBij$u0$minus*Log2*sqrtBij$u0 + 2*LBij$u0$plus*Log2*sqrtBij$u0 - 
           2*Lcosp*(4*Log2 + (LBij$u0$minus - LBij$u0$plus)*sqrtBij$u0) + 2*Laux$ni2$u0*(Lcosp + (-LBij$u0$minus + LBij$u0$plus)*sqrtBij$u0) + pow(Laux$ni2$u0,2) + 
           sqrtBij$u0*pow(LBij$u0$minus,2) - sqrtBij$u0*pow(LBij$u0$plus,2) - pow(Lcosp,2) - 8*pow(Log2,2) + pow(Pi,2)))*pow(sinp,-1)*pow(v,-1))/2. + 
   (v*pow(Pi,2)*pow(sinp,-1)*(mi2*(pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-2) + pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-2)) - 
        2*pipj*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1)))/2. - 
   (pipj*v*pow(sinp,-1)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1)*
      (2 - mi2*(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v))*pow(pipj,-1)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1) - 
        mi2*(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v))*pow(pipj,-1)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1))*
      (12*Log2*(-Laux1 + Lcosp + Log2)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1) + 
        3*(2*b*Laux2 + (-1 + b)*Laux3)*(Laux1 - Lcosp - 3*Log2)*pow(b,-1)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1) + 
        3*pow(-Laux1 + Lcosp + Log2,2)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1) + 
        (12*pow(Log2,2) - pow(Pi,2))*pow(pi$minus - 2*cosp*pTi*v + 2*pi$plus*pow(v,2),-1) + 
        (pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1)*(((6 + 6*b)*Laux2*Laux3 - 12*(Laux3*Laux6 + Paux1))*pow(b,-1) + 6*pow(Laux2,2) + 
             2*(-9.5 + 7*pow(b,-1))*pow(Pi,2) + 3*pow(b,-1)*pow(pipj,-1)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1)*
              pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1)*
              ((-1 + b)*mi2*pow(Laux5,2)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),2) + 
              (1 + b)*mi2*pow(Laux5,2)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),2) + 
              mi2*pow(Laux4,2)*((1 + b)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),2) + (-1 + b)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),2)))
              *pow(2 - mi2*(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v))*pow(pipj,-1)*pow(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v),-1) - 
              mi2*(pi$minus + 2*v*(-(cosp*pTi) + pi$plus*v))*pow(pipj,-1)*pow(pj$minus + 2*v*(-(cosp*cosxi*pTj) + pTj*sinp*sinxi*y + pj$plus*v),-1),-1)))/2.))/3.;
  return res;
}


double psi_contribution_DIST2t(double u, double v, double cosp, double y){
  double res = 0;
  res = ((-2*DIST2t_0(u,v,cosp,-y) - 2*DIST2t_0(u,v,cosp,y) + (1 + y)*(DIST2t_0(u,v,cosp,-1)
      + DIST2t_0(u,v,cosp,1))*(1 + Log2 - y*Log2))*pow(1 - y,-1)*pow(1 + y,-1))/2.;
  return res;
}

double psi_contribution_DISTt(double u, double v, double cosp, double y){
  double res = 0;
  res = ((-2*DISTt_0(u,v,cosp,-y) - 2*DISTt_0(u,v,cosp,y) + (1 + y)*(DISTt_0(u,v,cosp,-1)
      + DISTt_0(u,v,cosp,1))*(1 + Log2 - y*Log2))*pow(1 - y,-1)*pow(1 + y,-1))/2.;
  return res;
}

double psi_contribution_AUX2(double u, double v, double cosp, double y){
  double res = 0;
  res = ((-2*AUX2_0(u,v,cosp,-y) - 2*AUX2_0(u,v,cosp,y) + (1 + y)*(AUX2_0(u,v,cosp,-1)
      + AUX2_0(u,v,cosp,1))*(1 + Log2 - y*Log2))*pow(1 - y,-1)*pow(1 + y,-1))/2.;
  return res;
}

double psi_contribution_DELTAt_0(double v, double cosp, double y){
  double res = 0;
  res = ((-2*DELTAt_0(v,cosp,-y) - 2*DELTAt_0(v,cosp,y) + (1 + y)*(DELTAt_0(v,cosp,-1)
      + DELTAt_0(v,cosp,1))*(1 + Log2 - y*Log2))*pow(1 - y,-1)*pow(1 + y,-1))/2.;
  return res;
}

double psi_contribution_DELTAt_m1(double v, double cosp, double y){
  double res = 0;
  res = (DELTAt_m1(v,cosp,-1)*(3*Log2*(-2 - Log2 + y*Log2) + 6*log(1 - y) + (-1 + y)*pow(Pi,2))*pow(-1 + y,-1))/12. + 
   (DELTAt_m1(v,cosp,1)*(3*Log2*(-2 - Log2 + y*Log2) + 6*log(1 - y) + (-1 + y)*pow(Pi,2))*pow(-1 + y,-1))/12. + 
   DELTAt_m1(v,cosp,-y)*(log(0.25 - y/4.) + log(1 + y))*pow(1 - pow(y,2),-1) + DELTAt_m1(v,cosp,y)*(log(0.25 - y/4.) + log(1 + y))*pow(1 - pow(y,2),-1);
  return res;
}

double psi_contribution_DELTAt_m2(double v, double cosp, double y){
  double res = 0;
  res = (DELTAt_m2(v,cosp,-y)*(24*Log2*(log(1 - y) + log(1 + y)) - 2*(pow(Pi,2) + 12*pow(Log2,2)) - 6*pow(log(1 - y) + log(1 + y),2))*pow(1 - pow(y,2),-1))/12. + 
   (DELTAt_m2(v,cosp,y)*(24*Log2*(log(1 - y) + log(1 + y)) - 2*(pow(Pi,2) + 12*pow(Log2,2)) - 6*pow(log(1 - y) + log(1 + y),2))*pow(1 - pow(y,2),-1))/12. - 
   (DELTAt_m2(v,cosp,-1)*pow(-1 + y,-1)*(-12*Log2*log(2 - 2*y) + 3*log(4 - 4*y)*log(1 - y) + pow(Pi,2) - log(4)*pow(Pi,2) + y*log(4)*pow(Pi,2) + log(8)*pow(Pi,2) - 
        y*log(8)*pow(Pi,2) + 15*pow(Log2,2) + 19*pow(Log2,3) - 19*y*pow(Log2,3) - log(4)*pow(log(8),2) + y*log(4)*pow(log(8),2) + 12*z3 - 12*y*z3))/12. - 
   (DELTAt_m2(v,cosp,1)*pow(-1 + y,-1)*(-12*Log2*log(2 - 2*y) + 3*log(4 - 4*y)*log(1 - y) + pow(Pi,2) - log(4)*pow(Pi,2) + y*log(4)*pow(Pi,2) + log(8)*pow(Pi,2) - 
        y*log(8)*pow(Pi,2) + 15*pow(Log2,2) + 19*pow(Log2,3) - 19*y*pow(Log2,3) - log(4)*pow(log(8),2) + y*log(4)*pow(log(8),2) + 12*z3 - 12*y*z3))/12.;
  return res;
}

double psi_contribution_DELTAt(double v, double cosp, double y){
  double res = psi_contribution_DELTAt_0(v,cosp,y) + psi_contribution_DELTAt_m1(v,cosp,y) + psi_contribution_DELTAt_m2(v,cosp,y);
  return res;
}



double integrand(double u, double v, double cosp, double y){

  double res=0;
  
  evaluate_external(false);

  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  res = res + DIST2t(u,v,cosp);

  evaluate_extra_variables_reg();
  res = res + AUX2(u,v,cosp);
  
  evaluate_external(true);

  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  res = res + DIST2t(u,v,cosp);

  evaluate_extra_variables_reg();
  res = res + AUX2(u,v,cosp);

  return res/2;

}

double integrand_u_lt_1(double u, double v, double cosp, double y){

  double res=0;
  
  evaluate_external(false);

  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  res = res + DISTt(u,v,cosp);
    
  evaluate_external(true);

  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  res = res + DISTt(u,v,cosp);

  return res/2;

}

double integrand_u_lt_1_AUX2(double u, double v, double cosp, double y){

  double res=0;
  
  evaluate_external(false);
  
  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  
  evaluate_extra_variables_reg();
  res = res + AUX2(u,v,cosp);
    
  evaluate_external(true);
  
  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  
  evaluate_extra_variables_reg();
  res = res + AUX2(u,v,cosp);

  return res/2;

}

double integrand_psi(double u, double v, double cosp, double y){

  double res=0;

  evaluate_external(false);

  res = res + psi_contribution_DIST2t(u,v,cosp,y);
  res = res + psi_contribution_AUX2(u,v,cosp,y);
  res = res + psi_contribution_DISTt(1/u,v,cosp,y)/u/u ;
  res = res + psi_contribution_AUX2(1/u,v,cosp,y)/u/u;
  res = res + psi_contribution_DELTAt(v,cosp,y)/u/u ;
  
  evaluate_external(true);

  res = res + psi_contribution_DIST2t(u,v,cosp,y);
  res = res + psi_contribution_AUX2(u,v,cosp,y);
  res = res + psi_contribution_DISTt(1/u,v,cosp,y)/u/u ;
  res = res + psi_contribution_AUX2(1/u,v,cosp,y)/u/u;
  res = res + psi_contribution_DELTAt(v,cosp,y)/u/u ;

  return res/2;

}

double integrand_delta_u(double v, double cosp, double y){

  double res=0;
  double u=0.5; //dummy
  
  evaluate_external(false);
  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  evaluate_extra_variables_reg();
  evaluate_extra_variables_delta(v,cosp,y);
  
  res = DELTAt(v,cosp,y);
  
  evaluate_external(true);
  evaluate_int_variables(u,v,cosp,y);
  evaluate_extra_variables(u,v,cosp);
  evaluate_extra_variables_reg();
  evaluate_extra_variables_delta(v,cosp,y);

  res = res + DELTAt(v,cosp,y);

  return res/2;

}

double delta_psi (double *k, size_t dim, void *params)
{
  (void)(dim); /* avoid unused parameter warnings */
  (void)(params);

  // u in (1,inf), v in (0,inf), cosp in (-1,1)
  double u = 1/k[0], v = -log(k[1]), cosp = 2*k[2]-1;
  double jac = u*u/k[1]*2;

  // this would be the old piece, actually independent from y
  // since it comes from the delta terms in the y distribution
  double res;
  try{
    res = (integrand(u,v,cosp,1) + integrand(u,v,cosp,-1))/2. * jac;
  }
  catch (...){
    return 0.;
  }
  
  if (isnan(res) || isinf(res)){
    res = 0;
  }
  return res;
}

double distr_psi (double *k, size_t dim, void *params)
{
  (void)(dim); /* avoid unused parameter warnings */
  (void)(params);

  // u in (1,inf), v in (0,inf), cosp in (-1,1), y in (0,1)
  double u = 1/k[0], v = -log(k[1]), cosp = 2*k[2]-1, y = k[3];
  double jac = u*u/k[1]*2;

  // now I should add here as well the new terms which actually depend on y
  double res;
  try{
    res = integrand_psi(u,v,cosp,y) * jac;
  }
  catch (...){
    return 0.;
  }
  
  if (isnan(res) || isinf(res)){
    res = 0;
  }
  return res;
}

double u_lt_1 (double *k, size_t dim, void *params)
{
  (void)(dim); /* avoid unused parameter warnings */
  (void)(params);

  // u in (0,1), v in (0,inf), cosp in (-1,1)
  double u = pow(k[0],3), v = -log(k[1]), cosp = 2*k[2]-1;
  double jac = 1./k[1]*2 * 3*pow(k[0],2);

  // this would be the old piece, actually independent from y
  // since it comes from the delta terms in the y distribution
  double res;
  try{
    res = (integrand_u_lt_1(u,v,cosp,1) + integrand_u_lt_1(u,v,cosp,-1))/2. * jac;
  }
  catch (...){
    return 0.;
  }
  
  if (isnan(res) || isinf(res)){
    res = 0;
  }
  return res;
}

double u_lt_1_AUX2 (double *k, size_t dim, void *params)
{
  (void)(dim); /* avoid unused parameter warnings */
  (void)(params);

  // u in (0,1), v in (0,inf), cosp in (-1,1)
  double u = k[0], v = -log(k[1]), cosp = 2*k[2]-1;
  double jac = 1./k[1]*2;

  // this would be the old piece, actually independent from y
  // since it comes from the delta terms in the y distribution
  double res;
  try{
    res = (integrand_u_lt_1_AUX2(u,v,cosp,1) + integrand_u_lt_1_AUX2(u,v,cosp,-1))/2. * jac;
  }
  catch (...){
    return 0.;
  }
  
  if (isnan(res) || isinf(res)){
    res = 0;
  }
  return res;
}


double delta_u (double *k, size_t dim, void *params)
{
  (void)(dim); /* avoid unused parameter warnings */
  (void)(params);

  // v in (0,inf), cosp in (-1,1)
  double v = -log(k[0]), cosp = 2*k[1]-1;
  double jac = 1/k[0]*2;

  double res;
  try{
    res = (integrand_delta_u(v,cosp,1) + integrand_delta_u(v,cosp,-1))/2. * jac;
  }
  catch (...){
    return 0.;
  }
  
  if (isnan(res) || isinf(res)){
    res = 0;
  }
  return res;
}

double sum_all (double *k, size_t dim, void *params)
{
  (void)(dim); /* avoid unused parameter warnings */
  (void)(params);

  double res;
  try{
    res = delta_psi(k,dim,params)+distr_psi(k,dim,params)+delta_u(k,dim,params)+u_lt_1(k,dim,params)+u_lt_1_AUX2(k,dim,params);
  }
  catch (...){
    return 0.;
  }
  
  if (isnan(res) || isinf(res)){
    res = 0;
  }
  return res;
}

void
display_results (const char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
}

double coth(double x){
  return cosh(x)/sinh(x);
}

double regular_fitted(double v34){
  double a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14;
  a1 = -1049.587264678528;
  a2 = -84.11162191300883;
  a3 = -323.037342978441;
  a4 = -14.371579647998743;
  a5 = -33.069454531192505;
  a6 = 0.597806952805064;
  a7 = -2.7828058225701717;
  a8 = -1.601175437181741;
  a9 = -137.97192760395922;
  a10 = 273.57438294065133;
  a11 = -895.3590089478016;
  a12 = -24.871468909240022;
  a13 = 154.06904375688052;
  a14 = 0.0914833296862021;
  double result = (12*a1 + 12*a11*coth(log((1 - v34)*pow(1 + v34,-1))/2.) + 6*a12*coth(log((1 - v34)*pow(1 + v34,-1))/2.)*log((1 - v34)*pow(1 + v34,-1)) + 
     2*a9*coth(log((1 - v34)*pow(1 + v34,-1))/2.)*(-6*Li2(2*v34*pow(-1 + v34,-1)) + 6*log((1 + v34)*pow(1 - v34,-1))*log((-1 + pow(v34,-1))/2.) + 
        pow(Pi,2)) + 6*a10*coth(log((1 - v34)*pow(1 + v34,-1))/2.)*(log(-2*v34*pow(-1 + v34,-1))*log(-1 + 2*pow(1 + v34,-1)) + 2*pow(Pi,2)) + 
     3*a14*coth(log((1 - v34)*pow(1 + v34,-1))/2.)*log(pow(-1 + v34,2)*pow(1 + v34,-2))*(-3*pow(Pi,2) + pow(atanh(v34),2)) + 
     12*a3*pow(coth(log((1 - v34)*pow(1 + v34,-1))/2.),2) + 6*a4*log((1 - v34)*pow(1 + v34,-1))*pow(coth(log((1 - v34)*pow(1 + v34,-1))/2.),2) + 
     a8*log(-1 + 2*pow(1 + v34,-1))*(-6*Li2(2*v34*pow(-1 + v34,-1)) - 6*log(2*v34)*log(1 + v34) + 6*log(1 - v34)*log(-2*v34*(1 + v34)*pow(-1 + v34,-1)) + 
        13*pow(Pi,2))*pow(coth(log((1 - v34)*pow(1 + v34,-1))/2.),2) + 
     3*a6*log(pow(-1 + v34,2)*pow(1 + v34,-2))*(-3*pow(Pi,2) + pow(atanh(v34),2))*pow(coth(log((1 - v34)*pow(1 + v34,-1))/2.),2) - 
     12*(pow(Pi,2) - pow(atanh(v34),2))*(a2 + a13*coth(log((1 - v34)*pow(1 + v34,-1))/2.) + a5*pow(coth(log((1 - v34)*pow(1 + v34,-1))/2.),2)) + 
     2*a7*pow(coth(log((1 - v34)*pow(1 + v34,-1))/2.),2)*(6*Li3((1 - v34)*pow(1 + v34,-1)) - 2*log((1 - v34)*pow(1 + v34,-1))*pow(Pi,2) + 
        pow(log((1 - v34)*pow(1 + v34,-1)),3)))/12.;
  return -2*result;
}

void numerical_integrals_T3T4 (double *result, double *error, double results_rest)
{

  if (acc_goal_T3T4<0) {
    *result = 0.;
    *error  = 0.;
    return;
  }
  
  double v34 = sqrt(1-mi2*mj2/pipj/pipj);
  
  double acc_goal = acc_goal_T3T4;
  int max_iter = 10;
  
  double res, err, res_tot;
  double res1, err1, res2, err2, res3, err3, res4, err4, res5, err5, res6, err6;

  double xl_4[4] = { 0, 0, 0, 0 };
  double xu_4[4] = { 1, 1, 1, 1 };
  double xl_3[3] = { 0, 0, 0 };
  double xu_3[3] = { 1, 1, 1 };
  double xl_2[2] = { 0, 0 };
  double xu_2[2] = { 1, 1 };

  const gsl_rng_type *T1,*T2,*T3,*T4,*T5;
  gsl_rng *r1,*r2,*r3,*r4,*r5;

  gsl_monte_function G1 = { &delta_psi, 3, 0 };
  gsl_monte_function G2 = { &distr_psi, 4, 0 };
  gsl_monte_function G3 = { &delta_u, 2, 0 };
  gsl_monte_function G4 = { &u_lt_1, 3, 0 };
  gsl_monte_function G5 = { &u_lt_1_AUX2, 3, 0 };

  size_t calls = calls_T3T4;

  gsl_rng_env_setup ();

  T1 = gsl_rng_default;
  r1 = gsl_rng_alloc (T1);
  T2 = gsl_rng_default;
  r2 = gsl_rng_alloc (T2);
  T3 = gsl_rng_default;
  r3 = gsl_rng_alloc (T3);
  T4 = gsl_rng_default;
  r4 = gsl_rng_alloc (T4);
  T5 = gsl_rng_default;
  r5 = gsl_rng_alloc (T5);

  int fac1=2, fac2=5;

  bool compute_G1 = true, compute_G2 = true, compute_G3 = true, compute_G4 = true, compute_G5 = true, compute_G6 = false;
  res6 = regular_fitted(v34);
  err6 = 0.;

  {
    gsl_monte_vegas_state *s1 = gsl_monte_vegas_alloc (3);
    gsl_monte_vegas_state *s2 = gsl_monte_vegas_alloc (4);
    gsl_monte_vegas_state *s3 = gsl_monte_vegas_alloc (2);
    gsl_monte_vegas_state *s4 = gsl_monte_vegas_alloc (3);
    gsl_monte_vegas_state *s5 = gsl_monte_vegas_alloc (3);

    gsl_monte_vegas_params params1,params2,params3,params4,params5,params6;
    
    for(int i=0; i<2; i++){ // 2 warm-up iterations
      int fac[2] = { 3, 1 };
      gsl_monte_vegas_integrate (&G1, xl_3, xu_3, 3, calls/fac[i], r1, s1, &res1, &err1);
      gsl_monte_vegas_integrate (&G2, xl_4, xu_4, 4, calls/fac[i], r2, s2, &res2, &err2);
      gsl_monte_vegas_integrate (&G3, xl_2, xu_2, 2, calls*fac1/fac[i], r3, s3, &res3, &err3);
      gsl_monte_vegas_integrate (&G4, xl_3, xu_3, 3, calls*fac2/fac[i], r4, s4, &res4, &err4);
      gsl_monte_vegas_integrate (&G5, xl_3, xu_3, 3, calls*fac2/fac[i], r5, s5, &res5, &err5);
      if(false){std::cout << std::endl << res1 << " " << err1 << "   " << res2 << " " << err2 << "   " << res3 << " " << err3 << "   " << res4 << " " << err4 << "   " << res5 << " " << err5 << "   " << res6 << " " << err6 << std::endl;}
    }
    res = res1 + res2 + res3 + res4 + res5 + res6;
    err = sqrt(err1*err1 + err2*err2 + err3*err3 + err4*err4 + err5*err5);
    
    if(verbose){display_results ("vegas warm-up", res, err);}
    if(verbose){printf ("converging...\n");}

    gsl_monte_vegas_params_get(s1, &params1);
    gsl_monte_vegas_params_get(s2, &params2);
    gsl_monte_vegas_params_get(s3, &params3);
    gsl_monte_vegas_params_get(s4, &params4);
    gsl_monte_vegas_params_get(s5, &params5);

    int j = 0;
    do
      {
	if(compute_G1){gsl_monte_vegas_integrate (&G1, xl_3, xu_3, 3, calls, r1, s1, &res1, &err1);}
	if(compute_G2){gsl_monte_vegas_integrate (&G2, xl_4, xu_4, 4, calls, r2, s2, &res2, &err2);}
	if(compute_G3){gsl_monte_vegas_integrate (&G3, xl_2, xu_2, 2, calls*fac1, r3, s3, &res3, &err3);}
	if(compute_G4){gsl_monte_vegas_integrate (&G4, xl_3, xu_3, 3, calls*fac2, r4, s4, &res4, &err4);}
	if(compute_G5){gsl_monte_vegas_integrate (&G5, xl_3, xu_3, 3, calls*fac2, r5, s5, &res5, &err5);}
	// after first main iteration change stage to 2
	params1.stage = 2;
	params2.stage = 2;
	params3.stage = 2;
	params4.stage = 2;
	params5.stage = 2;
	gsl_monte_vegas_params_set(s1, &params1);
	gsl_monte_vegas_params_set(s2, &params2);
	gsl_monte_vegas_params_set(s3, &params3);
	gsl_monte_vegas_params_set(s4, &params4);
	gsl_monte_vegas_params_set(s5, &params5);
	if(abs(err1/res*100) < acc_goal/2){compute_G1 = false;}
	if(abs(err2/res*100) < acc_goal/2){compute_G2 = false;}
	if(abs(err3/res*100) < acc_goal/2){compute_G3 = false;}
	if(abs(err4/res*100) < acc_goal/2){compute_G4 = false;}
	if(abs(err5/res*100) < acc_goal/2){compute_G5 = false;}
	if(false){std::cout << res1 << " " << err1 << "   " << res2 << " " << err2 << "   " << res3 << " " << err3 << "   " << res4 << " " << err4 << "   " << res5 << " " << err5 << "   " << res6 << " " << err6 << std::endl;}
	res = res1 + res2 + res3 + res4 + res5 + res6;
	res_tot = results_rest + res*2/Pi*CA/8;
	err = sqrt(err1*err1 + err2*err2 + err3*err3 + err4*err4 + err5*err5);
	j = j+1;
	if(verbose){printf ("result = % .6f sigma = % .6f\n", res, err);}
      }
    while ((abs(err/res*100) > acc_goal && j < max_iter && abs((err*2/Pi*CA/8)/res_tot*100) > acc_goal/2.));
    
    if(verbose){display_results ("vegas final", res, err);}

    *result = res;
    *error = err;
    
    gsl_monte_vegas_free (s1);
    gsl_monte_vegas_free (s2);
    gsl_monte_vegas_free (s3);
    gsl_monte_vegas_free (s4);
    gsl_monte_vegas_free (s5);
  }

  gsl_rng_free (r1);
  gsl_rng_free (r2);
  gsl_rng_free (r3);
  gsl_rng_free (r4);
  gsl_rng_free (r5);

}


double L_integrand(double x, void *params){
  int n = *(int *) params;
  evaluate_external(false);
  double p2_x = pow(x,2)*mi2 +  pow(1-x,2)*mj2 + 2*x*(1-x)*pipj;
  double pT2_x = pow(x,2)*pTi*pTi +  pow(1-x,2)*pTj*pTj + 2*x*(1-x)*pTipTj;  
  double result = pow(log(1+pT2_x/p2_x),n)*pipj/p2_x;
  if(verbose){
    std::cout << "no crossing" << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "mi2 = " << mi2 << std::endl;
    std::cout << "mj2 = " << mj2 << std::endl;
    std::cout << "pipj = " << pipj << std::endl;
    std::cout << "p2_x   = " << p2_x << std::endl;
    std::cout << "pT2_x  = " << pT2_x << std::endl;
    std::cout << "result = " << result << std::endl;
  }

  evaluate_external(true);
  p2_x = pow(x,2)*mi2 +  pow(1-x,2)*mj2 + 2*x*(1-x)*pipj;
  pT2_x = pow(x,2)*pTi*pTi +  pow(1-x,2)*pTj*pTj + 2*x*(1-x)*pTipTj;  
  result = result + pow(log(1+pT2_x/p2_x),n)*pipj/p2_x;
  if(verbose){
    std::cout << "with crossing" << std::endl;
    std::cout << "x = " << x << std::endl;
    std::cout << "mi2 = " << mi2 << std::endl;
    std::cout << "mj2 = " << mj2 << std::endl;
    std::cout << "pipj = " << pipj << std::endl;
    std::cout << "p2_x   = " << p2_x << std::endl;
    std::cout << "pT2_x  = " << pT2_x << std::endl;
    std::cout << "result = " << result << std::endl;
  }
  // Divide by 2 to average
  if(verbose){std::cout << "final result = " << result/2 << std::endl;}
  return result/2;
}

double P_integrand(double x, void *params){
  int n = *(int *) params;
  evaluate_external(false);
  double p2_x = pow(x,2)*mi2 +  pow(1-x,2)*mj2 + 2*x*(1-x)*pipj;
  double pT2_x = pow(x,2)*pTi*pTi +  pow(1-x,2)*pTj*pTj + 2*x*(1-x)*pTipTj;
  double result;
  if(n==2){
    result = Li2(-pT2_x/p2_x)*pipj/p2_x;
  } else if (n==3){
    result = Li3(-pT2_x/p2_x)*pipj/p2_x;
  } else {
    std::cout << "Pn not supported" << std::endl;
    exit(1);
  }

  evaluate_external(true);
  p2_x = pow(x,2)*mi2 +  pow(1-x,2)*mj2 + 2*x*(1-x)*pipj;
  pT2_x = pow(x,2)*pTi*pTi +  pow(1-x,2)*pTj*pTj + 2*x*(1-x)*pTipTj;
  if(n==2){
    result = result + Li2(-pT2_x/p2_x)*pipj/p2_x;
  } else if (n==3){
    result = result + Li3(-pT2_x/p2_x)*pipj/p2_x;
  } else {
    std::cout << "Pn not supported" << std::endl;
    exit(1);
  }
  // Divide by 2 to average
  return result/2;
}

double integrate_1D(char flag, int n){
  
  gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);

  double result, error;

  gsl_function F;
  if(flag == 'L'){
    F.function = &L_integrand;
  } else if (flag == 'P'){
    F.function = &P_integrand;
  } else {
    exit(1);
  }
  F.params = &n;

  gsl_integration_qags (&F, 0, 1, 0, 1e-5, 1000,
                        w, &result, &error);

  gsl_integration_workspace_free (w);

  return result;

}

double S_tilde_jj(int n){
  double result = 0;
  double B,t;
  evaluate_external(false);
  for (int i=0; i<2; i++){
    if(i==0){
      B = pTi*pTi/mi2;
    } else {
      B = pTj*pTj/mj2;
    }
    t = sqrt(1+B);

    if(n == -2){
      result = result + 4;
    } else if (n==-1){
      result = result + 8*(log(1+B)-1);
    } else if (n==0){
      result = result + (2*(30 - 24*Li2(-B) - 24*log(1 + B) - 5*pow(Pi,2)))/3.;
    } else if (n==1){
      result = result + 2*((2*(-24 - 18*z3 - 12*Li3(-B) + 12*Li3(pow(1 + B,-1)) - 24*log(1 + B) - 3*pow(Pi,2) + log(1 + B)*pow(Pi,2) - 12*pow(log(1 + B),2) + 6*log(B)*pow(log(1 + B),2) - 
          4*pow(log(1 + B),3)))/3. + (-12 + 312*z3 + 254*Li2(1 - t) - 72*Li2(-t) + 24*Li2(pow(t,-2)) - 182*Li2(pow(t,-1)) + 56*Li2(2*pow(1 + t,-1)) - 96*Li3((1 - t)/2.) - 
          192*Li3(1 - t) - 192*Li3(-t) + 192*Li3((1 + t)*pow(1 - t,-1)) - 192*Li3(pow(t,-1)) + 96*Li3(((-1 + t)*pow(t,-1))/2.) - 96*Li3((-1 + t)*pow(t,-1)) + 
          96*Li3((1 + pow(t,-1))/2.) - 192*Li3(pow(1 + t,-1)) - 96*Li3(2*pow(1 + t,-1)) - 192*Li3((-1 + t)*pow(1 + t,-1)) - 96*Li3(t*pow(1 + t,-1)) + 
          8*Li2((1 - t)/2.)*(-7 + 12*log((-1 + t)/2.)) - 35*log(-1 + t) + 96*Li2(1 - t)*log(-1 + t) + 48*Li2(-t)*log(-1 + t) - 24*Li2(pow(t,-2))*log(-1 + t) + 
          144*Li2(pow(t,-1))*log(-1 + t) - 96*Li2(pow(1 + t,-1))*log(-1 + t) - 8*log(2)*log(-1 + t) + 108*log(4)*log(-1 + t) + 248*log(t) - 54*Li2(1 - t)*log(t) + 
          96*Li2(-t)*log(t) + 48*Li2(pow(t,-2))*log(t) - 96*Li2(pow(t,-1))*log(t) - 26*log(-1 + t)*log(t) + 96*log(2)*log(-1 + t)*log(t) + 96*Li2(-t)*log(2*t) - 
          96*Li2(2*pow(1 + t,-1))*log((1 + t)/2.) + 152*log(-1 + t)*log((1 + t)/2.) - 96*log(2)*log(t)*log((1 + t)/2.) + 96*log(-1 + t)*log(t)*log((1 + t)/2.) - 
          3*log(1 + t) + 54*Li2(1 - t)*log(1 + t) + 48*Li2(-t)*log(1 + t) + 24*Li2(pow(t,-2))*log(1 + t) - 144*Li2(pow(t,-1))*log(1 + t) - 96*Li2(pow(1 + t,-1))*log(1 + t) - 
          56*log((-1 + t)/2.)*log(1 + t) - 272*log(t)*log(1 + t) + 288*log(-1 + t)*log(t)*log(1 + t) + 96*Li2(((-1 + t)*pow(t,-1))/2.)*log(2*t*pow(-1 + t,-1)) - 
          96*Li2((-1 + t)*pow(t,-1))*log(2*t*pow(-1 + t,-1)) + 96*log(2)*log(t)*log(2*t*pow(-1 + t,-1)) + 192*Li2(1 - t)*log((1 + t)*pow(-1 + t,-1)) + 
          96*Li2(-t)*log((1 + t)*pow(-1 + t,-1)) + 48*Li2(pow(t,-2))*log((1 + t)*pow(-1 + t,-1)) - 96*Li2(pow(t,-1))*log((1 + t)*pow(-1 + t,-1)) + 
          96*Li2(pow(1 + t,-1))*log((1 + t)*pow(-1 + t,-1)) - 156*log(-1 + t)*log((1 + t)*pow(-1 + t,-1)) + 48*log(-1 + t)*log(t)*log((1 + t)*pow(-1 + t,-1)) - 
          156*log(1 + t)*log((1 + t)*pow(-1 + t,-1)) - 192*log(-1 + t)*log(1 + t)*log((1 + t)*pow(-1 + t,-1)) + 48*log(t)*log(1 + t)*log((1 + t)*pow(-1 + t,-1)) + 
          48*Li2(-t)*log(1 - pow(t,-2)) + 48*Li2(pow(t,-1))*log(1 - pow(t,-2)) - 52*log(-1 + t)*log(1 - pow(t,-2)) + 52*log(1 + t)*log(1 - pow(t,-2)) + 
          96*Li2(pow(t,-1))*log(pow(t,-1)/2.) + 32*log((-1 + t)*pow(t,-1)) - 6*Li2(1 - t)*log((-1 + t)*pow(t,-1)) - 48*Li2(-t)*log((-1 + t)*pow(t,-1)) - 
          96*Li2(pow(t,-1))*log((-1 + t)*pow(t,-1)) - 52*log(-1 + t)*log((-1 + t)*pow(t,-1)) + 8*log(t)*log((-1 + t)*pow(t,-1)) - 156*log(1 + t)*log((-1 + t)*pow(t,-1)) + 
          48*log(t)*log(1 + t)*log((-1 + t)*pow(t,-1)) - 96*Li2((1 + pow(t,-1))/2.)*log((1 + pow(t,-1))/2.) - 96*Li2(t*pow(1 + t,-1))*log((1 + pow(t,-1))/2.) - 
          96*log(2)*log((1 + t)/2.)*log((1 + pow(t,-1))/2.) - 48*Li2(1 - t)*log(1 + pow(t,-1)) - 48*Li2(-t)*log(1 + pow(t,-1)) + 52*log(-1 + t)*log(1 + pow(t,-1)) - 
          52*log(1 + t)*log(1 + pow(t,-1)) - 48*log(t)*log(1 + t)*log(1 + pow(t,-1)) + 3*log(-1 + pow(t,2)) + 48*Li2(1 - t)*log(-1 + pow(t,2)) - 
          48*Li2(-t)*log(-1 + pow(t,2)) - 24*Li2(pow(t,-2))*log(-1 + pow(t,2)) + 24*log(-1 + t)*log(-1 + pow(t,2)) + 48*log(t)*log(-1 + pow(t,2)) + 
          24*log(1 + t)*log(-1 + pow(t,2)) + 48*log(t)*log(1 + t)*log(-1 + pow(t,2)) + 156*log((1 + t)*pow(-1 + t,-1))*log(-1 + pow(t,2)) - 
          48*log(t)*log((1 + t)*pow(-1 + t,-1))*log(-1 + pow(t,2)) - 96*log(t)*log(1 + pow(t,-1))*log(-1 + pow(t,2)) + 54*Li2(1 - t)*log((-1 + t)*pow(1 + t,-1)) + 
          33*pow(Pi,2) - 16*log(2)*pow(Pi,2) + 70*log(-1 + t)*pow(Pi,2) - 80*log(t)*pow(Pi,2) - 16*log(1 + t)*pow(Pi,2) + 6*log(2*pow(-1 + t,-1))*pow(Pi,2) + 
          22*log(t*pow(-1 + t,-1))*pow(Pi,2) + 26*log(2*t*pow(-1 + t,-1))*pow(Pi,2) + 64*log((1 + t)*pow(-1 + t,-1))*pow(Pi,2) - 16*log(1 - pow(t,-2))*pow(Pi,2) - 
          32*log(pow(t,-1)/2.)*pow(Pi,2) + 32*log((-1 + t)*pow(t,-1))*pow(Pi,2) + 32*log((1 + pow(t,-1))/2.)*pow(Pi,2) - 32*log(1 + pow(t,-1))*pow(Pi,2) - 
          8*log(-1 + pow(t,2))*pow(Pi,2) + 6*Li2(pow(t,-2))*pow(t,-1) - 24*Li2(pow(t,-1))*pow(t,-1) + 16*log(1 - pow(t,-2))*pow(t,-1) - 16*log((-1 + t)*pow(t,-1))*pow(t,-1) + 
          12*log(t)*log((-1 + t)*pow(t,-1))*pow(t,-1) - 16*log(1 + pow(t,-1))*pow(t,-1) - 12*log(t)*log(1 + pow(t,-1))*pow(t,-1) - 56*pow(log(2),2) + 8*log(4)*pow(log(2),2) + 
          128*pow(log(-1 + t),2) + 96*log((1 + t)*pow(-1 + t,-1))*pow(log(-1 + t),2) + 85*pow(log(t),2) - 48*log(2)*pow(log(t),2) - 120*log(-1 + t)*pow(log(t),2) - 
          168*log(1 + t)*pow(log(t),2) + 48*log((1 + t)*pow(-1 + t,-1))*pow(log(t),2) + 24*log(1 - pow(t,-2))*pow(log(t),2) + 48*log(pow(t,-1)/2.)*pow(log(t),2) - 
          48*log((-1 + t)*pow(t,-1))*pow(log(t),2) - 48*log(-1 + pow(t,2))*pow(log(t),2) + 112*pow(log(t),3) + 28*pow(log((1 + t)/2.),2) + 48*log(2)*pow(log((1 + t)/2.),2) - 
          32*pow(log((1 + t)/2.),3) + 24*pow(log(1 + t),2) - 48*log(-1 + t)*pow(log(1 + t),2) + 96*log(t)*pow(log(1 + t),2) + 
          144*log((1 + t)*pow(-1 + t,-1))*pow(log(1 + t),2) - 16*pow(log(1 + t),3) - 48*log(2*t*pow(-1 + t,-1))*pow(log(t*pow(-1 + t,-1)),2) + 
          16*pow(log(t*pow(-1 + t,-1)),3) + 48*log(-1 + t)*pow(log(2*t*pow(-1 + t,-1)),2) + 48*log(pow(-1 + t,-1)/2.)*pow(log(2*t*pow(-1 + t,-1)),2) + 
          32*pow(log(2*t*pow(-1 + t,-1)),3) - 54*log(-1 + t)*pow(log((1 + t)*pow(-1 + t,-1)),2) + 48*log(t)*pow(log((1 + t)*pow(-1 + t,-1)),2) - 
          6*log(2*pow(-1 + t,-1))*pow(log((1 + t)*pow(-1 + t,-1)),2) - 54*log(t*pow(-1 + t,-1))*pow(log((1 + t)*pow(-1 + t,-1)),2) + 
          6*log(2*t*pow(-1 + t,-1))*pow(log((1 + t)*pow(-1 + t,-1)),2) - 64*pow(log((1 + t)*pow(-1 + t,-1)),3) + 48*log(2)*pow(log((1 + pow(t,-1))/2.),2) - 
          48*log((1 + pow(t,-1))/2.)*pow(log(1 + pow(t,-1)),2) + 16*pow(log(1 + pow(t,-1)),3) - 48*pow(log(-1 + pow(t,2)),2))/3.);
    } else {
      std::cout << "Unsupported S_tilde_jj" << std::endl;
      exit(1);
    }
  }
  return result/2;
}


double Fex2_and_11_T3T4 (void)
{
  // I'm computing here the result corresponding to
  // CA*Fex2_T3T4_CA + nf*Fex2_T3T4_nf - 2*Pi*beta0*Fex11_T3T4;
  // where Fex2_T3T4_CA, Fex2_T3T4_nf and Fex1_T3T4 are the
  // results in the grids in the old implementation
  evaluate_external(false);

  if(verbose){std::cout << "computing L integrals" << std::endl;}
  if(verbose){std::cout << "computing L0 integral" << std::endl;}
  double L0 = integrate_1D('L',0);
  if(verbose){std::cout << "computing L1 integral" << std::endl;}
  double L1 = integrate_1D('L',1);
  if(verbose){std::cout << "computing L2 integral" << std::endl;}
  double L2 = integrate_1D('L',2);
  //double L3 = integrate_1D('L',3);

  if(verbose){std::cout << "computing P integrals" << std::endl;}
  double P2 = integrate_1D('P',2);
  //double P3 = integrate_1D('P',3);

  if(verbose){std::cout << "computing S_tilde_jj" << std::endl;}
  double S_tilde_jj_1 = S_tilde_jj(1);

  double Bi, Bj;
  Bi = pTi*pTi/mi2;
  Bj = pTj*pTj/mj2;
  
  double Fex11_piece = -P2 + Pi*Pi/12*(L0-1) + Li2(-Bi)/2 + Li2(-Bj)/2;
  
  double result;
  result = (CA*(-808 + 808*L0 + 804*L1 + 198*L2 - 396*P2 - 27*S_tilde_jj_1
		  + (396*Li2(-Bi) - 804*log(1+Bi) - 198*pow(log(1+Bi),2))/2
		  + (396*Li2(-Bj) - 804*log(1+Bj) - 198*pow(log(1+Bj),2))/2
		  - 33*pow(Pi,2) + 33*L0*pow(Pi,2) ))/216.
    + (nf*(56 - 56*L0 - 60*L1 - 18*L2 + 36*P2
	     + (- 36*Li2(-Bi) + 60*log(1+Bi) + 18*pow(log(1+Bi),2))/2
	     + (- 36*Li2(-Bj) + 60*log(1+Bj) + 18*pow(log(1+Bj),2))/2
	     + 3*pow(Pi,2) - 3*L0*pow(Pi,2) ))/108.;
  result = result - 2*Pi*beta0*Fex11_piece;

  if(verbose){std::cout << "Computing reg34m[1]+reg34m0[1]+(1/2)R$[34,1]+sing34m[1]+sing34m0[1]" << std::endl;}
  double Fex2_integrals, err;    
  numerical_integrals_T3T4(&Fex2_integrals, &err, result);
  Fex2_integrals = Fex2_integrals*2/Pi*CA/8;
  err = err*2/Pi*CA/8;

  result = result + Fex2_integrals;

  if(false){std::cout << "numerical piece result: " << Fex2_integrals << " +/- " << err << std::endl;}
  if(print_error){std::cout << "results T3T4: " << result << ", " << err << std::endl;} // dominant error coming from MC integration
           
  return result;
}
