
double L34(void){
  // following eq 23 from 2102.03256
  // Different normalization: additional 1/2 here
  // I have i=3 and j=4 in the code
  // careful, in the paper they use jk and I'm using ij

  evaluate_external(false);
  
  double vij = sqrt(1-mi2*mj2/pipj/pipj); //eq 20
  double rij1 = sqrt(mj2/mi2)*p1pi/p1pj; //eq 24
  double rij2 = sqrt(mj2/mi2)*p2pi/p2pj; //eq 24

  double result = 1./2.*log((1+vij)/(1-vij))*log((mi2+pTi*pTi)*(mj2+pTj*pTj)/mi2/mj2)
    - 2.*Li2(2*vij/(1+vij)) - 1./4.*pow(log((1+vij)/(1-vij)),2)
    + Li2(1-sqrt((1-vij)/(1+vij))*rij1) + Li2(1-sqrt((1-vij)/(1+vij))/rij1) + 1./2.*pow(log(rij1),2)
    + Li2(1-sqrt((1-vij)/(1+vij))*rij2) + Li2(1-sqrt((1-vij)/(1+vij))/rij2) + 1./2.*pow(log(rij2),2); // eq 23

  return result/2;
}

double Fex1_T3T4(void){
  evaluate_external(false);
  double vij = sqrt(1-mi2*mj2/pipj/pipj);
  double Bi, Bj;
  Bi = pTi*pTi/mi2;
  Bj = pTj*pTj/mj2;
  double result = L34()/vij - (log(1+Bi)+log(1+Bj))/2.;
  return result;
}

double Fex1_TITj(double B){
  double result = -(log(1+B)+Li2(-B))/2.;
  return result;
}


