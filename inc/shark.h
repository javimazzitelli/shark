#ifdef __cplusplus
extern "C" {
#endif
  double H1_TT_shark(int *i, int *j, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int *nf_in);
  double H2_TT_shark(int *i, int *j, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int *nf_in);
  double H2_TTTT_shark(int *i, int *j, int *k, int *l, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int *nf_in);
  double H2_fTTT_shark(int *i, int *j, int *k, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int *nf_in);
  void setprecision_T3T4_shark(double *acc);
  void setprecision_fTTT_shark(double *acc);
  void setcalls_T3T4_shark(int *calls);
  void setprecision_mode_shark(int *mode);
#ifdef __cplusplus
}
#endif

double H1_TT(int i, int j, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int nf_in);
double H2_TT(int i, int j, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int nf_in);
double H2_TTTT(int i, int j, int k, int l, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int nf_in);
double H2_fTTT(int i, int j, int k, double *p1_in, double *p2_in, double *p3_in, double *p4_in, int nf_in);
void setprecision_T3T4(double acc);
void setprecision_fTTT(double acc);
void setcalls_T3T4(int calls);
void setprecision_mode(int mode);
