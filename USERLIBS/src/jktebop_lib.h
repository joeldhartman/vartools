/* This is a c header for the jktebop library. The jktebop program was developed by:

!=======================================================================
!     PROGRAM JKTEBOP             John Taylor   j.k.taylor@warwick.ac.uk
!                                 Astrophysics Group  Warwick University
!-----------------------------------------------------------------------

*/

void task1_(void);
void task2_(double V[], int LDTYPE[]);
void task34_(int *TASK, double V[], int VARY[], int LDTYPE[], double *DATA, int DTYPE[], int *NDATA, int *NLR, int *NMIN, double *SIGMA);
void task6_(double V[], int VARY[], int LDTYPE[], double *DATA, int DTYPE[], int *NDATA, int *NLR, int *NMIN);
void task578_(int *TASK, double V[], int VARY[], int LDTYPE[], double *DATA, int DTYPE[], int *NDATA, int *NLR, int *NMIN, int *NSIM);
void stopcheck_(int *STATUS);
void outputsim_(int *TASK, double V[], int VARY[], double VEXTRA[], double *VALL, int *NSIM);
void output_(double *DATA, int DTYPE[], int *NDATA, int *NLR, int *NMIN, double V[], int VARY[], int LDTYPE[], int *ITER, double *CHISQ, double VERR[], int *WHAT);
double getphase_(double *HJD, double *PERIOD, double *TZERO);
double getmin_(double *TZERO, double *PERIOD, double *ECC, double *OMEGA, double *CYCLE);
double getmodel_(double V[], int LDTYPE[], double *TIME, int *TYPE, double *LA, double *LB);
void fitebop_(double *DATA, int DTYPE[], int *NDATA, int *NLR, int *NMIN, double V[], int VARY[], int LDTYPE[], int *ITER, double *CHISQ, double VERR[], int *ERROR);
void ebop_(int *INDEX, double *X, double V[], double *Y, double DYDA[], int *NCOEFFS, int VARY[]);
void biax_(double *R, double *Q, double *A, double *B, double *EPS);
void getew_(double *ECOSW, double *ESINW, double *E, double *W);
void light_(double V[], int LDTYPE[], double *PHASE, double *FMAG, float *LP, float *LS);
void mrqmin_(double X[], double Y[], double SIG[], int *ndata, double A[], int IA[], int *ma, double *covar, double *alpha, int *nca, double *chisq, double *alamda, int *ifail);
void mrqcof_(double x[], double y[], double sig[], int *ndata, double a[], double *alpha, double beta[], int *nalp, double *chisq);
void gaussj_(double *a, int *n, int *np, double *b, int *m, int *mp, int *ifail);
void covsrt_(double *covar, int *npc, int *ma, int ia[], int *mfit);
int seedstart_(void);
double random_(int *SEED);
double randomg_(int *SEED, double *MEAN, double *SD);
double select_(double ARRAY[], int *NUM, int *K);
double sigma_(double ARRAY[], int *NUM);
