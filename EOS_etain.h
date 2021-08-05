// allocation dyn.
void freetab(void *ptr);
int **alloctab(int dim1, int dim2);
double **alloctabd(int dim1, int dim2);

// LOI D'ETAT
// x
double fx(double V, int phase);
// THETA
double THETA(double V, int phase);
// Es
double fEs(double V, int phase);
// Ps
double fPs(double V, int phase);
// u
double u(double V, double E, int phase);
// S
double fS(double V, double E, int phase);
// P
double fP(double V, double E, int phase);
// T
double fT(double V, double E, int phase);
// Es'
double Es_prime(double V, int phase);
// Ps'
double Ps_prime(double V, int phase);
// dP/dE
double dPdE(int phase);
// dPdV
double dPdV(double V, int phase);
// dT/dE
double dTdE(int phase);
// dTdV
double dTdV(double V, int phase);

void coeff(int phase, double* K0, double* N0, double* gamma0, double* Crv, double* theta0, double* T0, double* P0, double* rho0, double* v0, double* E0, double* Sr);

double fE_VP(double V, double P, int phase);