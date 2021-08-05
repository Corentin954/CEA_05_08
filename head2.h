// allocation dyn.
void freetab(void *ptr);
int **alloctab(int dim1, int dim2);
double **alloctabd(int dim1, int dim2);

// LOI D'ETAT
double epsilonEOS(int tst, double tau, double p);
void EGS(int tst, double tau, double epsilon, double* Pres, double* Cres);

// Montée en odre spatiale
double phipw(int so, int ind, double* Z, double* X, double** rk);
double fxc(int so,double* X, double** rk);
double phibar(int so, int ind, double* Z, double* X, double** rk);
double delta(int so,int ind, double* Z, double** dk);
double YdZ(int so, int cas, double* P, double* Q, double* u, double* X, double** rk);
double phiQK(int so, double* DK, double* X,double** rk);

void initCoef(double** Ck, double** Ckbar, double** dk, double** Qbar, double** rk);

//
void u0(int tst, double a, double b,double x, double* W);

// Table RUNGE - KUTTA
void initButcher(int sch, double** A, double* THETA, double* ALPHA);

void condlim(int tst, int jdeb, int jfin, int nbghosts, double* Z, int iu);

// The Minmod limiter
double MinMod(double r);

// The Superbee limiter
double Superbee(double r);

// The Christensen limiter
double Christensen(double r, double arg);

// pseusdo-viscosity
double qvis(int q, double Cq, double Cl, double* U, double tau, double c, double* DM);