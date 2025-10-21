#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "shark.h"

using namespace std;

void print_momenta(double *p1, double *p2, double *p3, double *p4){
  cout << "p1 = (" << p1[0] << ", " << p1[1] << ", " << p1[2] << ", " << p1[3] << ")" << endl;
  cout << "p2 = (" << p2[0] << ", " << p2[1] << ", " << p2[2] << ", " << p2[3] << ")" << endl;
  cout << "p3 = (" << p3[0] << ", " << p3[1] << ", " << p3[2] << ", " << p3[3] << ")" << endl;
  cout << "p4 = (" << p4[0] << ", " << p4[1] << ", " << p4[2] << ", " << p4[3] << ")" << endl;
}

int main (void)
{
  // setting precision mode to 1 (default)
  setprecision_mode(1);
  // number of light flavours
  int nf = 5;
  bool verbose = false;
  // read the kinematics from file events.dat
  int max_events = 100;
  double p1_list[max_events][4],p2_list[max_events][4],p3_list[max_events][4],p4_list[max_events][4];
  int X1 = 0;
  
  ifstream inFile ("events.dat");
  
  while ( !inFile.eof () ) {
    inFile >> p1_list [X1][0];
    inFile >> p1_list [X1][1];
    inFile >> p1_list [X1][2];
    inFile >> p1_list [X1][3];
    inFile >> p2_list [X1][0];
    inFile >> p2_list [X1][1];
    inFile >> p2_list [X1][2];
    inFile >> p2_list [X1][3];
    inFile >> p3_list [X1][0];
    inFile >> p3_list [X1][1];
    inFile >> p3_list [X1][2];
    inFile >> p3_list [X1][3];
    inFile >> p4_list [X1][0];
    inFile >> p4_list [X1][1];
    inFile >> p4_list [X1][2];
    inFile >> p4_list [X1][3];
    X1++;
    }
  
  X1--;
  int total_points = X1;

  for (int i=0; i<total_points; i++){

    double p1_in[4],p2_in[4],p3_in[4],p4_in[4];
    p1_in[0] = p1_list[i][0];
    p1_in[1] = p1_list[i][1];
    p1_in[2] = p1_list[i][2];
    p1_in[3] = p1_list[i][3];
    p2_in[0] = p2_list[i][0];
    p2_in[1] = p2_list[i][1];
    p2_in[2] = p2_list[i][2];
    p2_in[3] = p2_list[i][3];
    p3_in[0] = p3_list[i][0];
    p3_in[1] = p3_list[i][1];
    p3_in[2] = p3_list[i][2];
    p3_in[3] = p3_list[i][3];
    p4_in[0] = p4_list[i][0];
    p4_in[1] = p4_list[i][1];
    p4_in[2] = p4_list[i][2];
    p4_in[3] = p4_list[i][3];

    print_momenta(p1_in,p2_in,p3_in,p4_in);

    // evaluate the different colour structures of H1:
    cout << "H1_TT(3,4) = " << H1_TT(3,4,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H1_TT(3,1) = " << H1_TT(3,1,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H1_TT(4,2) = " << H1_TT(4,2,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H1_TT(4,1) = " << H1_TT(4,1,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H1_TT(3,2) = " << H1_TT(3,2,p1_in,p2_in,p3_in,p4_in,nf) << endl;

    // evaluate the different colour structures of H2:
    // Two-parton correlators
    cout << "H2_TT(3,4) = " << H2_TT(3,4,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TT(3,1) = " << H2_TT(3,1,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TT(4,2) = " << H2_TT(4,2,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TT(4,1) = " << H2_TT(4,1,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TT(3,2) = " << H2_TT(3,2,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    // Three-parton correlator (only one structure)
    cout << "H2_fTTT(1,2,3) = " << H2_fTTT(1,2,3,p1_in,p2_in,p3_in,p4_in,nf) << endl;

    // Four-parton correlators
    cout << "H2_TTTT(3,4,3,4) = " << H2_TTTT(3,4,3,4,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    
    cout << "H2_TTTT(3,4,3,1) = " << H2_TTTT(3,4,3,1,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TTTT(3,4,3,2) = " << H2_TTTT(3,4,3,2,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TTTT(3,4,4,1) = " << H2_TTTT(3,4,4,1,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TTTT(3,4,4,2) = " << H2_TTTT(3,4,4,2,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    
    cout << "H2_TTTT(3,1,3,1) = " << H2_TTTT(3,1,3,1,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TTTT(3,1,3,2) = " << H2_TTTT(3,1,3,2,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TTTT(3,1,4,1) = " << H2_TTTT(3,1,4,1,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TTTT(3,1,4,2) = " << H2_TTTT(3,1,4,2,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TTTT(3,2,3,2) = " << H2_TTTT(3,2,3,2,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TTTT(3,2,4,1) = " << H2_TTTT(3,2,4,1,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TTTT(3,2,4,2) = " << H2_TTTT(3,2,4,2,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TTTT(4,1,4,1) = " << H2_TTTT(4,1,4,1,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TTTT(4,1,4,2) = " << H2_TTTT(4,1,4,2,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    cout << "H2_TTTT(4,2,4,2) = " << H2_TTTT(4,2,4,2,p1_in,p2_in,p3_in,p4_in,nf) << endl;
    
    cout << endl;
             
  }
  return 0;
}
