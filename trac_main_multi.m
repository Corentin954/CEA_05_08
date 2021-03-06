clear all




%///////////////             RUNGE-KUTTA               ///////////////////
%///////////////    Plot  schema + sol analytique        ///////////////////


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choix de la loi d'etat et du cas test
%  tst : 0 (gaz parfait)   2
tst = 200; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deux_plaques=floor(tst/100);

if tst==0  % Sod
 T=0.2; % pas en ligne de commande
 nx=200;
elseif tst==1  % LeBlanc
  %T=6;
  nx=400;
  %a=0; b=9;
elseif tst==2  % Bizarrium
  T=80e-6; % pas en ligne de commande
  nx=1000;
elseif tst==3  % Onde acoustique
  T=1/2; % pas en ligne de commande
  nx=200;
elseif tst==4  % Sod sym?tris?
  %T=0.2;
  nx=400;
elseif tst==5  % Woodward
  %T=0.038;
  nx=300;
  %a=0; b=1;
elseif tst==6  % Shu-Osher
  %T=1.8;
  nx=300;
  %a=-5; b=5;
elseif tst==100   % etain
  T=2e-7;
  nx=200;
elseif deux_plaques  % etain
  T=1e-7;
  nx=300;
elseif tst==200  % etain
  T=1e-7;
  nx=300;
else
 disp(["erreur valeur tst= (Script Octave)",num2str(tst)]);
end




% =======================================================================
% Creation de l'executable
system("sh compil_main_multi.sh");


##%%%%%%%%%%%%%%%%%%%%%%%%  Presentation mi-stage
##Lname=['Godunov acous o1'; 
##       'RK1, q o1 q+l, so 2, kin fix'; 
##       'BBC JCP 2009, q o1 q+l';
##       'SSP RK3, q o1 q+l, so 2, kin fix'; ];
##
##Lsch=[002,100,200,102];  % scheme
##Lso=[-1,2,-1,2];    % spatial order
##Lz=[-1,1,-1,1];     % kinetic energy fix
##Lq=[4,4,4,4];      % artificial viscosity
##Ldpi=[-1,0,-1,0];   % dPI vs delta PI bar

##%%%%%%%%%%%%%%%%%%%%%%%%  Presentation mi-stage
##Lname=['Godunov acous o1';
##       'Godunov acous o2, GAD MinMod'; 
##       "vNR, q o1 Q+L";
##       'BBC JCP 2009, q o1 Q+L'; 
##       'BBC Pred-Corr, q o1 Q+L'; 
##       'BBC RK2av, q o1 Q+L';
##       'RK1, q o1 Q+L, so 3, kin fix'; 
##       'RK5 Dormand-Prince, q o1 Q+L, so 2, kin fix';
##       'RK5 Dormand-Prince, q o1 Q+L, so 3, kin fix';
##       'RK5 Dormand-Prince, q o1 Q+L, so 5, kin fix';
##       'RK5 Dormand-Prince, q o1 Q+L, so 7, kin fix'; 
##        ];
##
##Lsch=[002,011,300,200,201,202,100,105,105,105,105];  % scheme
##Lso=[-1,-1,-1,-1,-1,-1,3,2,3,5,7];    % spatial order
##Lz=[-1,-1,-1,-1,-1,-1,1,1,1,1,1];     % kinetic energy fix
##Lq=  (004)*[-1,-1,1,1,1,1,1,1,1,1,1];      % artificial viscosity
##Ldpi=[-1,-1,-1,-1,-1,-1,0,0,0,0,0];     % dPI vs delta PI bar
##LTtau=1e-8*[1,1,1,1,1,1,1,1,1,1,1];
##Lcin_phase=0*[1,1,1,1,1,1,1,1,1,1,1];

%%%%%%%%%%%%%%%%%%%%%%%%  Presentation mi-stage
##Lname=['Godunov acous o2, GAD MinMod, equi thermo'; 
##       "vNR, q o1 Q+L, Ttau=1e-8";
##       'BBC JCP 2009, q o1 Q+L, Ttau=5e-9'; 
##       'BBC Pred-Corr, q o1 Q+L, Ttau=5e-9'; 
##       'BBC RK2av, q o1 Q+L, Ttau=5e-9';
##       'RK3, q o1 Q+L, so 2, kin fix, Ttau=5e-9'; 
##       'RK3, q o1 Q+L, so 3, kin fix, Ttau=5e-9'; 
##       'RK3, q o1 Q+L, so 5, kin fix, Ttau=5e-9';
##       'RK3, q o1 Q+L, so 7, kin fix, Ttau=5e-9';
##        ];
##
##Lsch=[011,300,200,201,202,102,102,102,102];  % scheme
##Lso=[-1,-1,-1,-1,-1,2,3,5,7];    % spatial order
##Lz=[-1,-1,-1,-1,-1,1,1,1,1];     % kinetic energy fix
##Lq=  (004)*[-1,1,1,1,1,1,1,1,1];      % artificial viscosity
##Ldpi=[-1,-1,-1,-1,-1,0,0,0,0,0];     % dPI vs delta PI bar
##LTtau=5e-9*[0,1,1,1,1,1,1,1,1];
##Lcin_phase=0*[0,1,1,1,1,1,1,1,1];

%%%%%%%%%%%%%%%%%%%%%%  Godunov energie totale
##Lname=["vNR, q o1 Q+L";
##        'Godunov acous o1'; 
##        'Godunov Despres'; 
##        'Godunov Jaouen'; 
##        'Godunov acous o2, without lim'; 
##        'Godunov acous o2, lim MinMod'; 
##        'Godunov acous o2, lim Superbee'; ] ;
##
##Lsch=[002,001,000,010,011,012];  % scheme
##Lso= [-1,-1,-1,-1,-1,-1];    % spatial order
##Lz=  [-1,-1,-1,-1,-1,-1];     % kinetic energy fix
##Lq=  [-1,-1,-1,-1,-1,-1];      % artificial viscosity
##Ldpi=[-1,-1,-1,-1,-1,-1];   % dPI vs delta PI bar
##LTtau=1e-8*[1,1,1,1,1,1,1,1,1,1];
##Lcin_phase=[1,1,1,1,1,1,1,1,1,1];


%%%%%%%%%%%%%%%%%%%%%%  Spiteri-Ruuth
##Lname=[ 'Godunov acous o1, sans cin';  
##        'Godunov acous o1, cin phase Ttau=1e-8';
##        'Godunov acous o1, cin phase Ttau=1e-9';
##        'Godunov acous o1, cin phase Ttau=1e-10';
##        'Godunov acous o1, cin phase Ttau=1e-11';  ] 
##
##Lsch=[002,002,002,002,002];  % scheme
##Lso= [-1,-1,-1,-1,-1];    % spatial order
##Lz=  [-1,-1,-1,-1,-1];     % kinetic energy fix
##Lq=  -1*[-1,1,1,1,1];      % artificial viscosity
##Ldpi=[-1,-1,-1,-1,-1];   % dPI vs delta PI bar
##LTtau=[-1,1e-8,1e-9,1e-10,1e-11];
##Lcin_phase=[0,1,1,1,1];

##%%%%%%%%%%%%%%%%%%%%%%  vNR
####Lname=['vNR, q o1 Q+L, equi thermo';
####       'vNR, q o1 Q+L, cin phase Ttau=1e-9';  ];
####
####Lsch=[300,300];   % scheme
####Lso= [-1,-1];     % spatial order
####Lz=  [-1,-1];     % kinetic energy fix
####Lq=  4*[1,1];     % artificial viscosity
####Ldpi=[-1,-1];   % dPI vs delta PI bar
####LTtau=1e-9*[-1,1];
####Lcin_phase=[0,1];

%%%%%%%%%%%%%%%%%%%%%%  Godunov energie totale
##Lname=['Godunov acous o1, ';   
##       'Godunov acous o2 lim Superbee, ';   
##       'Godunov acous o2 lim MinMod, ';   ];
##       
##       
##Lsch=[002,012,011];   % scheme
##Lso= [-1,-1,-1];     % spatial order
##Lz=  [-1,-1,-1];     % kinetic energy fix
##Lq=  [-1,-1,-1];     % artificial viscosity
##Ldpi=[-1,-1,-1];     % dPI vs delta PI bar
##LTtau=5e-9*[1,1,1];
##Lcin_phase=1*[1,1,1];


##%%%%%%%%%%%%%%%%%%%%%%%%  RUNGE-KUTTA
##Lname=['RK1, q o1 Q+L, so 3, with out kin fix, tau=5e-9';];
##
##Lsch=[100];  % scheme
##Lso= [3];    % spatial order
##Lz=  [0];     % kinetic energy fix
##Lq=  (004)*[1];      % artificial viscosity
##Ldpi=[0];     % dPI vs delta PI bar
##LTtau=5e-9*[1];
##Lcin_phase=0*[1];

%%%%%%%%%%%%%%%%%%%%%%%%  BBC
##Lname=['Godunov acous o1, equi thermo'; 
##       'Godunov acous o2 DAD MinMod, equi thermo'; 
##       'BBC JCP 2009, q o1 Q+L, equi thermo';
##       'BBC PredCorr, q o1 Q+L,equi thermo';
##       'BBC RK2av, q o1 Q+L, equi thermo'; ];
##
##Lsch=[002,011,200,201,202];  % scheme
##Lso= [-1,-1,-1,-1,-1];    % spatial order
##Lz=  [-1,-1,-1,-1,-1];     % kinetic energy fix
##Lq=  (004)*[-1,-1,1,1,1];      % artificial viscosity
##Ldpi=[-1,-1,-1,-1,-1];     % dPI vs delta PI bar
##LTtau=0*[1,1,1,1,1];
##Lcin_phase=0*[1,1,1,1,1];

Lname=[ 'Godunov acous o1, sans cin';  ] 

Lsch=[002];  % scheme
Lso= [-1];    % spatial order
Lz=  [-1,];     % kinetic energy fix
Lq=  -1*[-1];      % artificial viscosity
Ldpi=[-1,];   % dPI vs delta PI bar
LTtau=[-1];
Lcin_phase=[0];

ls=length(Lsch);


% CFL
numCFL=0.2;
 
Ntau=15;

fac_dt=1e-2;

% Coeff de pseudo
Cq=1.5;
Cl=0.15;
%Cl=0.29;

% affichage
aff=1;

sch_cin_phase=0;

% Multi fluide
tst_multi=4;

% Max iter
   Nmax=1e5;

% symoble pour la legende
Symleg=['.-'; '-.'; '+-';'--', '*-'; 'x-'; 'o-'; 'd-'; 's-';'h-';'v-';':' ];


Etot=zeros(ls, 1e5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:ls
  disp([' ']);
  Lname(j,:)
  ["a.exe  -a ",num2str(aff)," -q ",num2str(Lq(j))," -d ",num2str(Cq)," -l ",num2str(Cl)," -p ",num2str(Ldpi(j))," -o ",num2str(Lso(j))," -z ",num2str(Lz(j))," -c ",num2str(Nmax)," -n ",num2str(nx)," -t ",num2str(tst)," -s ",num2str(Lsch(j))," -m ",num2str(numCFL)," -i ",num2str(Ntau)," -u ",num2str(LTtau(j))," -w ",num2str(Lcin_phase(j))," -x ",num2str(fac_dt)," -y ",num2str(sch_cin_phase)]
  system(["a.exe  -a ",num2str(aff)," -q ",num2str(Lq(j))," -d ",num2str(Cq)," -l ",num2str(Cl)," -p ",num2str(Ldpi(j))," -o ",num2str(Lso(j))," -z ",num2str(Lz(j))," -c ",num2str(Nmax)," -n ",num2str(nx)," -t ",num2str(tst)," -s ",num2str(Lsch(j))," -m ",num2str(numCFL)," -i ",num2str(Ntau)," -u ",num2str(LTtau(j))," -w ",num2str(Lcin_phase(j))," -x ",num2str(fac_dt)," -y ",num2str(sch_cin_phase)," -r ",num2str(tst_multi)])
  
  resu=load("sol.txt");
  [~,n]=size(resu);
  %size(resu)
  Xd(j,:) = resu(1,:); % X d?centr?e (fronti?res)
  Xc(j,:) = resu(2,1:n-1);
  tau(j,:) = resu(3,1:n-1);
  u(j,:) = resu(4,:);
  e(j,:) = resu(5,1:n-1);
  p(j,:) = resu(6,1:n-1);
  epsilon(j,:) = resu(7,1:n-1);
  
  % energie du syt?me
##  NRJ=load("NRJ.txt");
##  Nbiter(j)=length(NRJ(1,:));
##  Etot(j,1:Nbiter(j))=NRJ(1,:);
##  IMPUL(j,1:Nbiter(j))=NRJ(2,:);
##  Mtot(j,1:Nbiter(j))=NRJ(3,:);
##  VOL(j,1:Nbiter(j))=NRJ(4,:);
  
  resu=load("frac_mass.txt");
  [~,n]=size(resu);
  fmA(j,:)=resu(1,1:n-1);
  fmB(j,:)=resu(2,1:n-1);
  fmC(j,:)=resu(3,1:n-1);
  
endfor

% Diagramme V,E
SA=load("fichiers/pSA.txt");
SB=load("fichiers/pSB.txt");
SC=load("fichiers/pSC.txt");
SAB=load("fichiers/pSAB.txt");
SAC=load("fichiers/pSAC.txt");
SBC=load("fichiers/pSBC.txt");


[na,~]=size(SA); % de taille NA+1 car SA[NA]=S[0] (code C)
[nb,~]=size(SB);
[nc,~]=size(SC);
[nab,~]=size(SAB); % de taille NBA+1 car SA[NA]=S[0] (code C)
[nac,~]=size(SAC);
[nbc,~]=size(SBC);

% Traingle triple
Va=129.905671875844291e-6;  Ea=72.168190127265447e3;
Vb=127.666112185063994e-6;  Eb=99.480256119236515e3;
Vc=132.538478805442338e-6;  Ec=136.154346228414312e3;



%%  
if 1 %0
  
  figure
  hold on
  for j=1:ls
    plot(tau(j,:)*1e6,epsilon(j,:),Symleg(j,:));
  endfor
  % Triangle
  plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'r');
  % Ploynomes
  S=SA;
  plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
  S=SB;
  plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
  S=SC;
  plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne sp?cfique ,  t=',num2str(T),' s , nx=',num2str(nx)]);
  xlabel('V (cm^3/kg')
  ylabel("E  (J/Kg)");
  if ls==1
    legend(Lname(1,:),'location','northeast');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','northeast');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','northeast');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','northeast');
  elseif ls==5
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','northeast');
  elseif ls==6
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','northeast');
  elseif ls==7
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
  elseif ls==8
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
  elseif ls==9
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
  elseif ls==10
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
  elseif ls==11
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
  endif  
  hold off

  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),epsilon(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne sp?cfique ,  t=',num2str(T),' s , nx=',num2str(nx)]);
  xlabel('m')
  ylabel("E  (J/Kg)");
  if ls==1
    legend(Lname(1,:),'location','northeast');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','northeast');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','northeast');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','northeast');
  elseif ls==5
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','northeast');
  elseif ls==6
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','northeast');
  elseif ls==7
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
  elseif ls==8
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
  elseif ls==9
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
  elseif ls==10
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
  elseif ls==11
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
  endif  
  hold off

  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),tau(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\tau nx=',num2str(nx)]);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  elseif ls==7
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
  elseif ls==8
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
  elseif ls==9
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
  elseif ls==10
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
  elseif ls==11
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
  endif  
  hold off
  
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),1./tau(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\rho nx=',num2str(nx)]);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  elseif ls==7
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
  elseif ls==8
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
  elseif ls==9
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
  elseif ls==10
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
  elseif ls==11
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
  endif  
  hold off

  figure
  hold on
  for j=1:ls
    plot(Xd(j,:),u(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['u nx=',num2str(nx)]);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  elseif ls==7
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
  elseif ls==8
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
  elseif ls==9
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
  elseif ls==10
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
  elseif ls==11
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
  endif  
  hold off

  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),p(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['p nx=',num2str(nx)]);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  elseif ls==7
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
  elseif ls==8
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
  elseif ls==9
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
  elseif ls==10
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
  elseif ls==11
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
  endif  
  hold off
  
  
  % fraction massique
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),fmA(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase A nx=',num2str(nx)]);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  elseif ls==7
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
  elseif ls==8
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
  elseif ls==9
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
  elseif ls==10
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
  elseif ls==11
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
  endif  
  hold off
  
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),fmB(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase B nx=',num2str(nx)]);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  elseif ls==7
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
  elseif ls==8
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
  elseif ls==9
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
  elseif ls==10
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
  elseif ls==11
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
  endif  
  hold off
  
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),fmC(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['phase C nx=',num2str(nx)]);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  elseif ls==7
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
  elseif ls==8
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
  elseif ls==9
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
  elseif ls==10
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
  elseif ls==11
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
  endif  
  hold off
  
  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),fmA(j,:)+fmB(j,:)+fmC(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Somme des fractions massiques nx=',num2str(nx)]);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  elseif ls==7
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
  elseif ls==8
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
  elseif ls==9
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
  elseif ls==10
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
  elseif ls==11
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
  endif  
  hold off
  
##
##  % Plot de l'energie du syst?me au cours du temps
##  figure
##  hold on
##  for j=1:ls
##    t=0:T/(Nbiter(j)-1):T;
##    plot(t,Etot(j,1:Nbiter(j)),'.-');
##  end
##  grid();
##  set(gca, "fontsize", 40);
##  title(['energie totale Bizarrium']);
##  if ls==1
##    legend(Lname(1,:),'location','southwest');
##  elseif ls==2
##    legend(Lname(1,:),Lname(2,:),'location','southwest');
##  elseif ls==3
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
##  elseif ls==4
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
##  elseif ls==5
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
##  elseif ls==6
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
##  elseif ls==7
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
##  elseif ls==8
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
##  elseif ls==9
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
##  elseif ls==10
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
##   elseif ls==11
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
##   endif  
##  hold off
##  
##  
##  % Plot de l'impulsion au cours du temps
##  figure
##  hold on
##  for j=1:ls
##    t=0:T/(Nbiter(j)-1):T;
##    plot(t,IMPUL(j,1:Nbiter(j)),'.-');
##  end
##  grid();
##  set(gca, "fontsize", 40);
##  title(['impulsion']);
##  if ls==1
##    legend(Lname(1,:),'location','southwest');
##  elseif ls==2
##    legend(Lname(1,:),Lname(2,:),'location','southwest');
##  elseif ls==3
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
##  elseif ls==4
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
##  elseif ls==5
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
##  elseif ls==6
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
##  elseif ls==7
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
##  elseif ls==8
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
##  elseif ls==9
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
##  elseif ls==10
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
##   elseif ls==11
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
##   endif  
##  hold off
##  
##  % Plot de la masse totale au cours du temps
##  figure
##  hold on
##  for j=1:ls
##    t=0:T/(Nbiter(j)-1):T;
##    plot(t,Mtot(j,1:Nbiter(j)),'.-');
##  end
##  grid();
##  set(gca, "fontsize", 40);
##  title(['masse totale']);
##  if ls==1
##    legend(Lname(1,:),'location','southwest');
##  elseif ls==2
##    legend(Lname(1,:),Lname(2,:),'location','southwest');
##  elseif ls==3
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
##  elseif ls==4
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
##  elseif ls==5
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
##  elseif ls==6
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
##  elseif ls==7
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
##  elseif ls==8
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
##  elseif ls==9
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
##  elseif ls==10
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
##   elseif ls==11
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
##   endif  
##  hold off
##  
##  % Plot du volume total au cours du temps
##  figure
##  hold on
##  for j=1:ls
##    t=0:T/(Nbiter(j)-1):T;
##    plot(t,VOL(j,1:Nbiter(j)),'.-');
##  end
##  grid();
##  set(gca, "fontsize", 40);
##  title(['volume total']);
##  if ls==1
##    legend(Lname(1,:),'location','southwest');
##  elseif ls==2
##    legend(Lname(1,:),Lname(2,:),'location','southwest');
##  elseif ls==3
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
##  elseif ls==4
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
##  elseif ls==5
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
##  elseif ls==6
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
##  elseif ls==7
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','southwest');
##  elseif ls==8
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),'location','southwest');
##  elseif ls==9
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),'location','southwest');
##  elseif ls==10
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),'location','southwest');
##   elseif ls==11
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),Lname(8,:),Lname(9,:),Lname(10,:),Lname(11,:),'location','southwest');
##   endif  
##  hold off
##  
  
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('ETAIN_EQU_THERMO_nx_1e4_t_112.mat');