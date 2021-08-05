// allocation dyn.
void freetab(void *ptr);
int **alloctab(int dim1, int dim2);
double **alloctabd(int dim1, int dim2);

// LOI D'ETAT
double epsilonEOS(int tst, double tau, double p);
void EGS(int tst, double tau, double epsilon, double* Pres, double* Cres);

// Montée en odre spatiale
double phi(int so, double* Z, double** Coef);
double phidiv(int so, double* RZ, double* RHO, double** Coef);
double fxc(int so,double* X, double** rk);
double delta(int so,int ind, double* Z, double** dk);
double YdZ(int so, int cas, double* P, double* Q, double* u, double** dk, double** Ckbar);
double phiQK(int so, double* DK, double** Q);
void initCoef(double** Ck, double** Ckbar, double** dk, double** Qbar, double** rk);

//
int u0(int tst, double a, double b,double x, double* W);

// Table RUNGE - KUTTA
int initButcher(int sch, double** A, double* THETA, double* ALPHA);

int condlim(int ind_cond, int jdeb, int jfin, int nbghosts, double* Z, int iu);

// The Minmod limiter
double MinMod(double r);

// The Superbee limiter
double Superbee(double r);

// The Christensen limiter
double Christensen(double r, double arg);

// pseusdo-viscosity
double qvis(int q, double Cq, double Cl, double* U, double tau, double c, double* DM);

int print_sol(int ind_u, int ideb, int ifin, double* X, double* Xc, double* TAU, double* U, double* E, double* P, double* EPS);

int print_nrj(int n, double* Etot, double* IMPUL, double* Mtot, double* VOL);

void print_err(int err, int Rsch, int sch);



double pastemps(double T, double t, double dt, double dt_old);



int funcGtot(int sch, int tst, double T, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu);

int funcRKint(int sch, int tst, double T, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu, int q , double Cq, double Cl, int z, int dpi, int so);

int funcBBC_JCP2009(int tst, double T, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu, int q , double Cq, double Cl);

int funcBBC_PRED_CORR(int tst, double T, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu, int q , double Cq, double Cl);

int funcBBC_RK2av(int tst, double T, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu, int q , double Cq, double Cl);

int funcvNR(int tst, double T, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu, int q , double Cq, double Cl);

int funcKOVA(int tst, double T, double a, double b, int nx, int Nmax, double CFL, int aff, int ind_cond, int iu, int z);