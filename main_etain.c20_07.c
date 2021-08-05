#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "head_etain.h"
#include "head.h"


/*  A FAIRE
  Premièrement lire le fichier qui donnes les tailles des fichiers pour allocation
  (tester les valeurs retours des allocations)

  Creer un fonction qui va lire les fichiers contenant
  les discretisations des lignes de chgmt de phase
  et les sommets des polygones
  Ainsi que la tabulation

  (le calcul se fera dans une autre fonction, cela permet de ne calculer qu'une fois les frontières
   et la tabulation pour plusieurs lancements du main : 
   int ecr_bench...()      )
  

  int lec_tabulation(double)
  Reprendre la fonction use_bench.. sauf qu'il faut allouer la mémoire de la tabulation et des frontières dans le main
  ( -->  Implique de faire un fichier qui contient la taille des fichiers)
  puis on transmet leur adresse à lec_tab

  Ensuite on transmet ces adresses au fonctions_schemas qui vont lancer les focntion P(V,E) et C(V,E)  (T ???)
   
  
  +

  Creation d'une fonction qui renvoie P et C à parir d'un couple(V,E) et des adresses des tableaux contenant la tabulation
  (modif de lec_bench_PT)

  ---------

  Imaginer une méthode qui grace a un argument en ligne de commande du main permettent de choisir 
  l'implementation de P et C que l'on souhaite




  IMPLEMENTATION DE LA METHODE DU CALCUL ENTIER A CHAQUE COUPLE (V,E) :

  On garde la lecture des frontières 
  On transmet les frontières au main_schemas 
  Puis on cree une focntion qui fait le test sur chaque zone et renvoie une erreur si on est hors domaine
  int calcul_PT



  IMPLEMENTATION DE LA METHODE DE LA TABULATION CONTENANT LES COUPLES (P,T) :

  On a besoin uniquement de l'adresse de la tabulation de P et de T
  Puis in cree une fonction qui va chercher la bonne valeur selon le couple (V,E)





  DANS LA FOCNTION ecr_bench_PT :   ---> MODIFIER LA FOCNTION graine... 
                                         pour qu'elle prenne en argument un entier corespondant à la ligne de ghmt de phase
                                         et les sommets des polygones (et non pas segVa12 etc ...
                                         Permet de limité les adresses de tableaux a transmettre au sommets de polygones


*/





/////////////////////  MAIN  ////////////////////
int main(int argc, char* argv[]){
  double T;
	int nx=200;
	double CFL=0.20; 
	double a, b;
	int Nmax=2e6;
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
                   printf("                    00 : Sod\n");
                   printf("                    01 : LeBlanc\n");
                   printf("                    02 : Bizarrium\n"); 
                   printf("                    03 : Onde acoustique\n");
                   printf("                    04 : Sod symétrisé\n");
                   printf("                    05 : Woodward (3 états)\n");
                   printf("                    06 : Shu-Oscher (2 états avec sinusoide)\n");
                   printf("                Etain:\n");
                   printf("                    10 : Projection d'une plaque d'un milimetre contre un mur a gauche\n");
                   printf("                    11 : Choc symetrisé\n");

    	  	  			 printf(" -s hydro scheme\n"); 
    	  	  			 printf("    	 Godunov energie totale\n");
    	  	  			 printf("                    000 : Després\n");
    	  	  			 printf("                    001 : Jaouen\n");
    	  	  			 printf("                    002 : Solveur acoustique ordre 1\n");
    	  	  			 printf("                    Solveur acoustique ordre 2\n");
    	  	  			 printf("                    010 : sans limiteur\n");
    	  	  			 printf("                    011 : limiteur MinMod\n");
    	  	  			 printf("                    012 : limiteur Superbee\n");
    	  	  			 printf("      Runge-Kutta energie interne\n");
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
                   printf("    	 BBC  \n");
                   printf("                    200: BBC JCP 2009\n");
                   printf("                    201: BBC Predictor-Corrector\n");
                   printf("                    202: BBC RK2 average\n");
                   printf("    	 von Neumann-Ritchmyer  \n");
                   printf("                    300: vNR\n");
                   printf("    	 Cauchy-Kovalevskaya  \n");
                   printf("                    400: Cauchy-Kovalevskaya\n");

    	  	  			 printf(" -q pseudo viscosity \n");
                   printf("      Ordre 1 : \n");
                   printf("        00 : von Neumann-Ritchmyer\n");
                   printf("        01 : Rosenbluth\n");
                   printf("        02 : Landschoff\n");
                   printf("        03 : Magical q combination\n");
                   printf("        04 : Quadratique + lineaire\n");
                   printf("        05 : off\n");
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
  
  int etain=tst/10;
  printf("etain= %d\n",etain );

// ########   Lecture des fichiers binaires contenants les discretisations des frontières  ############

  // Discrétisation du benchmark/tabulation
  int nb_points;
  int Nv_dis,Ne_dis;
  double Vmin_dis, Vmax_dis;
  double Emin_dis, Emax_dis;

  // Ecriture dans un fichier texte des paramètres de la discretisation
  /* pour pouvoir utiliser la discretisations juste en lisant les fichiers
  Organisation du fichier :
    nb_points
    Nv_dis , Ne_dis
    Vmin_dis , Vmax_dis
    Emin_dis , Emax_dis
  */

  FILE *fparams;
  if((fparams = fopen("fichiers/PARAMS_dis.txt", "r")) == NULL){printf("erreur ouverture fichier PARAMS_dis.txt \n"); return 1;}
  fscanf(fparams, "%d \n",&nb_points );
  fscanf(fparams, "%d %d\n",&Nv_dis,&Ne_dis );
  fscanf(fparams, "%lf %lf\n",&Vmin_dis,&Vmax_dis );
  fscanf(fparams, "%lf %lf\n",&Emin_dis,&Emax_dis );
  fclose(fparams);
  

  int** tabZONE=alloctab(5,(Nv_dis-1)*(Ne_dis-1));
  int* tab_Nb_ZONE=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(int)); 
  double* tabVx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEx=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));
  double* tabVy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double)); double* tabEy=malloc((Nv_dis-1)*(Ne_dis-1)*sizeof(double));

  // Nombre de points sur les frontières
  int NTG1, NPB1, NPB2, NTD2, NTD3, NPH3, NTG3;
  int N12, N13, N23, NA, NB, NC, NAB, NAC, NBC;

  fnb_pts(nb_points, &NTG1, &NPB1, &NPB2, &NTD2 , &NTD3, &NPH3, &NTG3, &N12, &N13, &N23, &NA, &NB, &NC, &NAB, &NAC, &NBC );

  double** SA=alloctabd(NA+1,2);
  double** SB=alloctabd(NB+1,2);
  double** SC=alloctabd(NC+1,2); // on ajoute S[N]=S[0]  --> utile pour l'algo
  double** SAB=alloctabd(NAB+1,2);
  double** SAC=alloctabd(NAC+1,2);
  double** SBC=alloctabd(NBC+1,2);
  
  if(etain){
    printf("Entree lec_tabulation\n");
    err=lec_tabulation(nb_points, Nv_dis, Ne_dis, tabZONE, tab_Nb_ZONE, tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
    if(err){printf("Erreur d'ouvertures des fichiers binaires (lec_tabulation)\n"); return 1;}
  }
    

// ##############################################################################

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
  else if (tst==10){
    a=0.;  b=1.0e-3;
    T=2e-7;
  }
  else if (tst==11){
    a=0.;  b=2.0e-3;
    T=2e-7;
  }
  else {
    printf("Erreur choix tst (main)\n");
    return 1;
  }

  // params pour les cond aux limites
  // changement de signe pour u si iu==1
  // Condition de type mur
  //  ind_cond=0;
  //  iu=1;   // indice pour la vitesse
  // Condition de type ? (libre)
  //  ind_cond=0;
  //  iu=0; 
  // Condition de type périodique
  //  ind_cond=1;
  //  iu=0; 


  // params pour les cond aux limites
  // changement de signe pour u si iu==1
  // Condition de type mur à gauche et libre à droite

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
  else if (tst==10){  
    ind_cond=2;
    iu=1;   // indice pour la vitesse
  }
  else if (tst==11){  
    ind_cond=0;
    iu=0;   // indice pour la vitesse
  }
  else {
    printf("Erreur dans le choix de tst (main)\n");
    return 1;
  }


  int Rsch=floor(sch/100);
  sch=sch-100*Rsch;

  if(Rsch==0){
    //int funcGtot(sch, tst, T, a, b, nx, Nmax, CFL, q , Cq, Cl);
    printf("Entree funcGtot_etain\n");
    err=funcGtot_etain(sch, tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu,
                 nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                 tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
    if(err){print_err(err,Rsch, sch); return 1;}
  }
  else if(Rsch==1){
    err=funcRKint_etain(sch, tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, q , Cq, Cl, z, dpi, so,
                  nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                  tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
    if(err){print_err(err,Rsch, sch); return 1;}
  }
  else if(Rsch==2){
    if(sch==0){
      err=funcBBC_JCP2009_etain(tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, q , Cq, Cl,
                          nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis,  tabZONE, tab_Nb_ZONE, 
                        
                          tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
      if(err){print_err(err,Rsch, sch); return 1;}
    }

    else if(sch==1){
      err=funcBBC_PRED_CORR_etain(tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, q , Cq, Cl,
                            nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                            tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
      if(err){print_err(err,Rsch, sch); return 1;}
    }

    else if(sch==2){
      err=funcBBC_RK2av_etain(tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, q , Cq, Cl,
                        nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                        tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
      if(err){print_err(err,Rsch, sch); return 1;}
    }
    else{printf("Erreur sch (main)\n");}
  }

  else if(Rsch==3){
    err=funcvNR_etain(tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, q , Cq, Cl,
                nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
    if(err){print_err(err,Rsch, sch); return 1;}
  }

  else if(Rsch==4){
    err=funcKOVA_etain(tst, T, a, b, nx, Nmax, CFL, aff, ind_cond, iu, z,
                 nb_points, Nv_dis, Ne_dis, Vmin_dis, Vmax_dis, Emin_dis, Emax_dis, tabZONE, tab_Nb_ZONE, 
                 tabVx, tabEx, tabVy, tabEy, SA, SB, SC, SAB, SAC, SBC);
    if(err){print_err(err,Rsch, sch); return 1;}
  }
  else{printf("erreur choix sch (main)\n"); return 1;}
  
  return 0;
}
