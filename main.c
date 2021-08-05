#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head.h"




/////////////////////  MAIN  ////////////////////
int main(int argc, char* argv[]){
  double T;
	int nx=200;
	double CFL=0.3; 
	double a, b;
	int Nmax=2e9;
	int tst=0; 
	int sch=0;	// choix du calcul de Ustar et Pstar (0:Solveur acoustique  1:Despres)
    int aff=1;
    int z=0;
    int q=0;
    int so=2; // spatial order
    int dpi=0; //choix du dPI pour l'intégration de rho.u  
    double Cq=1.5; //--> artificial vicosity q
    double Cl=0.15;
    int iu=0;
    int err=0;
    int ind_cond;


	for (int i = 0; i < argc; i++){
	  if (argv[i][0] == '-'){
	  	switch(argv[i][1]){
	  	  case 'h' : printf("Usage : \n"); 
    	  	  			 printf(" -t test problem \n");
                   printf("                    0 : Sod\n");
                   printf("                    1 : LeBlanc\n");
                   printf("                    2 : Bizarrium\n"); 
                   printf("                    3 : Onde acoustique\n");
                   printf("                    4 : Sod symétrisé\n");
                   printf("                    5 : Woodward (3 états)\n");
                   printf("                    6 : Shu-Oscher (2 états avec sinusoide)\n");
    	  	  			 printf(" -s hydro scheme\n"); 
    	  	  			 printf("    	Godunov energie totale\n");
    	  	  			 printf("                    000 : Després\n");
    	  	  			 printf("                    001 : Jaouen\n");
    	  	  			 printf("                    002 : Solveur acoustique ordre 1\n");
    	  	  			 printf("                    Solveur acoustique ordre 2\n");
    	  	  			 printf("                    010 : sans limiteur\n");
    	  	  			 printf("                    011 : limiteur MinMod\n");
    	  	  			 printf("                    012 : limiteur Superbee\n");
    	  	  			 printf("       Runge-Kutta energie interne\n");
    	  	  			 printf("                    100: RK1 (Matsuno forward-backward)\n");
                   printf("                    101: RK2 (Heun's method)\n");
                   printf("                    102: RK3 SSP\n");
                   printf("                    103: RK4 Classique\n");
                   printf("                    104: RK5 Cash-Karp\n");
                   printf("                    105: RK5 Dormand-Prince\n");
                   printf("                    106: RK3 HEUN\n");
                   printf("                    107: RK3 RALSTON\n");
                   printf("                    108: RK3 Bogaki-Shampine\n");
                   printf("                    109: RH2 Bogaki-Shampine\n");
                   printf("                    110: KUTTA ordre 3\n");
                   printf("                    111: RK2 Explicit Midpoint Method\n");
                   printf("                    112: RK2 Ralston\n");
                   printf("                    113: EULER forward\n");
                   printf("                    114: SSPRK4(5,4) Spiteri-Ruuth\n");
                   printf("                    115: SSPRK3(4,3) Spiteri-Ruuth\n");
                   printf("                    116: SPPRK3(5,3) Spiteri-Ruuth\n");
                   printf("                    117: SPPRK1(3,1) Spiteri-Ruuth\n");
                   printf("                    118: SPPRK2(3,2) Spiteri-Ruuth\n");
                   printf("                    119: SPPRK2(4,2) Spiteri-Ruuth\n");
                   printf("                    120: SPPRK1(2,1) Spiteri-Ruuth\n");
                   printf("                    121: SPPRK3(4,3) Spiteri-Ruuth\n");   
                   printf("    	BBC  \n");
                   printf("                    200: BBC JCP 2009\n");
                   printf("                    201: BBC Predictor-Corrector\n");
                   printf("                    202: BBC RK2 average\n");
                   printf("    	von Neumann-Ritchmyer  \n");
                   printf("                    300: vNR\n");
                   printf("    	Cauchy-Kovalevskaya  \n");
                   printf("                    400: Cauchy-Kovalevskaya\n");


    	  	  			 printf(" -q pseudo viscosity \n");
                   printf("      Ordre 1 : \n");
                   printf("        00 : von Neumann-Ritchmyer\n");
                   printf("        01 : Rosenbluth\n");
                   printf("        02 : Landschoff\n");
                   printf("        03 : Magical q combination\n");
                   printf("        04 : off\n");
                   printf("      Ordre 2 : \n");
                   printf("        10 : von Neumann-Ritchmyer\n");
                   printf("        11 : Rosenbluth\n");
                   printf("        12 : Landschoff\n");
                   printf("        13 : Magical q combination\n");
                   printf("        14 : off\n");

      	  		     printf(" -o spatial order \n");
                   printf("        2 : 2nd order\n");
                   printf("        3 : 3rd order\n");
                   printf("        4 : 4th and 5th order\n");
                   printf("        5 : 6th and 7th order\n");
                   printf("        6 : 8th and 9th order\n");
                 
                   printf(" -n number cells \n");
    	  	  	    		 
          //       printf(" -f final time \n");
           			 
                   printf(" -c final cycle \n");

                   printf(" -k CFL \n");

                   printf(" -a affichage des cycles \n"); 
                   printf("           0 : off\n");
                   printf("           1 : on\n");

                   printf(" -z kinétic energy fix\n");
                   printf("        0 : off\n");
                   printf("        1 : CRAS 2016\n");
                   printf("        2 : JCP Dakin et al. 2019\n");
                   printf("      A chaque sous-pas (Runge-Kutta) :\n");
                   printf("        11 : CRAS 2016\n");
                   printf("        21 : JCP Dakin et al. 2019\n");
                 
                   printf(" -p options\n");
                   printf("        0 : dPI\n");
                   printf("        1 : AV(delta PI)\n");
                 
                   printf(" -d options\n");
                   printf("        valeur du coefficient pseudo-quadratique\n");
                 
                   printf(" -l options\n");
                   printf("        valeur du coefficent pseudo-lineaire\n");
    	  	  			 return 0;
	  	  case 't' : tst=atoi(argv[++i]);
	  	  			     break;
	  	  case 's' : sch=atoi(argv[++i]);
	  	  			     break;
        case 'o' : so=atoi(argv[++i]);
                   break;
        case 'k' : CFL=strtod(argv[++i],NULL);
                   break;
	  	  case 'n' : nx=atoi(argv[++i]);
	  	  			     break;
        case 'a' : aff=atoi(argv[++i]);
                   break;
	  	  case 'c' : Nmax=atoi(argv[++i]);
	  	  			     break;
        case 'z' : z=atoi(argv[++i]);
                   break;
        case 'p' : dpi=atoi(argv[++i]);
                   break;
        case 'q' : q=atoi(argv[++i]);
                   break;
        case 'd' : Cq=strtod(argv[++i],NULL);
                   break;
        case 'l' : Cl=strtod(argv[++i],NULL);
                   break;
	  	  default :  printf("Not recognized option : %c\n",argv[i][1] );
	  	  			     return 1;
	  	}
	  }
	}



	printf("sch= %d\n",sch );
  printf("tst= %d\n",tst );
  printf("CFL= %lf\n",CFL );
  printf("aff= %d\n",aff );
  printf("Nmax= %d\n",Nmax );

  if (tst==0){
    a=0.;  b=1.0;
    T=0.2;
  }
  else if (tst==2){
    a=0.;  b=1.0;  
    T=80e-6; 
  }
  else if(tst==3){
    a=0.;  b=1.0;  
    T=1.0;
  }
  else if(tst==4){
    a=-1.0;  b=1.0;  
    T=0.2;
  }
  else if(tst==1){
    a=0.0;  b=9.0;  
    T=6;
  }
  else if(tst==5){
    a=0.0;  b=1.0;  
    T=0.038;
  }
  else if(tst==6){
    a=-5.0;  b=5.0;  
    T=1.8;
  }
  else {
    printf("Erreur choix tst (main)\n");
    return 1;
  }

  // params pour les cond aux limites
  // changement de signe pour u si iu==1
  
  // Condition de type mur
  if (tst==0|| tst==4 || tst==1 || tst==5){  
    ind_cond=0;
    iu=1;   // indice pour la vitesse
  }
  // Condition de type ? (libre)
  else if (tst==2|| tst==6){
    ind_cond=0;
    iu=0; 
  }
  // Condition de type périodique
  else if (tst==3){
    ind_cond=1;
    iu=0; 
  }
  else {
    printf("Erreur dans le choix de tst (main)\n");
    return 1;
  }


  int Rsch=floor(sch/100);
  sch=sch-100*Rsch;

  if(Rsch==0){
    //int funcGtot(sch, tst, T, a, b, nx, Nmax, CFL, q , Cq, Cl);
    err=funcGtot(sch, tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu);
    if(err){print_err(err,Rsch, sch); return 1;}
  }

  else if(Rsch==1){
    err=funcRKint(sch, tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, q , Cq, Cl, z, dpi, so);
    if(err){print_err(err,Rsch, sch); return 1;}
  }

  else if(Rsch==2){
    if(sch==0){
      err=funcBBC_JCP2009(tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, q , Cq, Cl);
      if(err){print_err(err,Rsch, sch); return 1;}
    }

    else if(sch==1){
      err=funcBBC_PRED_CORR(tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, q , Cq, Cl);
      if(err){print_err(err,Rsch, sch); return 1;}
    }

    else if(sch==2){
      err=funcBBC_RK2av(tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, q , Cq, Cl);
      if(err){print_err(err,Rsch, sch); return 1;}
    }
    else{printf("Erreur sch (main)\n");}
  }

  else if(Rsch==3){
    err=funcvNR(tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, q , Cq, Cl);
    if(err){print_err(err,Rsch, sch); return 1;}
  }

  else if(Rsch==4){
    err=funcKOVA(tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, z);
    if(err){print_err(err,Rsch, sch); return 1;}
  }
  else{printf("erreur choix sch (main)\n"); return 1;}
  
  return 0;
}
