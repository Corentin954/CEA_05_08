clear all

a=0; b=1;

gamma=1.4;


%///////////////             RUNGE-KUTTA               ///////////////////
%///////////////    Plot  schema + sol analytique        ///////////////////


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choix de la loi d'etat et du cas test
%  tst : 0 (gaz parfait)   2
tst = 2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
elseif tst==4  % Sod symétrisé
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
end

if tst==0
  % test de Sod
  rhoG=1.0; rhoD=0.125;
  uG=0; uD=0;
  pG=1.0; pD=0.1;
  
  wL=[rhoG; uG; pG];
  wR=[rhoD; uD; pD];
  xdis=(b+a)/2;
  
  % Calcul de P et U
  system("sh compil_PU.sh");
  system("a.exe");
  resu=load("PU.txt");
  p12=resu(1)
  u12=resu(2)
endif




% =======================================================================
% Creation de l'executable
system("sh compil_main.sh");



%%%%%%%%%%%%%%%%%%%%%%%%  RKint
##Lname=['RK1 , q o1 vNR, so 2, kin fix CRAS 2016'; ];
##
##Lsch=[100];
##Lso=[2]; % spatial order
##Lz=[1];
##Lq=[0];
##Ldpi=[0];

%%%%%%%%%%%%%%%%%%%%%%%%  RKint
##Lname=['RK1 , q o1 vNR, so 3, kin fix CRAS 2016'; ];
##
##Lsch=[100];
##Lso=[3]; % spatial order
##Lz=[1];
##Lq=[0];
##Ldpi=[0];

%%%%%%%%%%%%%%%%%%%%%%%%  Gtot
##Lname=['Godunov Despres ';
##       'Godunov Jaouen ';
##       'Godunov Sol acous o1 '; ];
##
##Lsch=[000,001,002];
##Lso=[-1,-1,-1]; % spatial order
##Lz=[-1,-1,-1];
##Lq=[-1,-1,-1];
##Ldpi=[-1,-1,-1];

%%%%%%%%%%%%%%%%%%%%%%%%  Solveur Acoustique ordre 2
##Lname=['Godunov Sol acous o2 sans limiteur';
##       'Godunov Sol acous o2 lim MinMod';
##       'Godunov Sol acous o2 lim Superbee'; ];
##
##Lsch=[010,011,012];
##Lso=[-1,-1,-1]; % spatial order
##Lz=[-1,-1,-1];
##Lq=[-1,-1,-1];
##Ldpi=[-1,-1,-1];

####%%%%%%%%%%%%%%%%%%%%%%%%  BBC et vNR
##Lname=['BBC JCP 2009, q o1 vNR';
##       'BBC Pred Corr, q o1 vNR';
##       'BBC RK2 average, q o1 vNR';];
##
##Lsch=[200,201,202];
##Lso=[-1,-1,-1]; % spatial order
##Lz=[-1,-1,-1];
##Lq=[00,00,00];
##Ldpi=[-1,-1,-1];

%%%%%%%%%%%%%%%%%%%%%%%%  Kovaleskaya
##Lname=['Cauchy-Kovaleskaya, kin fix'];
##
##Lsch=[400];
##Lso=[-1]; % spatial order
##Lz=[1,];
##Lq=[-1];
##Ldpi=[-1];

%%%%%%%%%%%%%%%%%%%%%%%%  BBC et vNR
##Lname=['BBC JCP 2009, q o1 vNR';
##       'Godunov Sol acous o1';];
##
##Lsch=[200,002];
##Lso=[-1,-1]; % spatial order
##Lz=[-1,-1];
##Lq=[00,00];
##Ldpi=[-1,-1];


%%%%%%%%%%%%%%%%%%%%%%%%  
##Lname=['Godunov Sol acous o1';];
##
##Lsch=[002];
##Lso=[-1]; % spatial order
##Lz=[-1];
##Lq=[-1];
##Ldpi=[-1];

####%%%%%%%%%%%%%%%%%%%%%%%%  BBC
##Lname=['BBC JCP 2009, q o1 vNR';];
##
##Lsch=[200];
##Lso=[-1]; % spatial order
##Lz=[-1];
##Lq=[00];
##Ldpi=[-1];

%%%%%%%%%%%%%%%%%%%%%%%%  Solveur Acoustique ordre 2
##Lname=['Godunov Sol acous o1';
##       'RK1, q o1 vNR, so 2, kin fix CRAS 2016'; 
##       'BBC JCP 2009, q o1 vNR';
##       'vNR, q o1 vNR';
##       'RK5 Cash-Karp, q o1 vNR, so 5, kin fix CRAS 2016';];
##
##Lsch=[002,100,200,300,104]; % scheme
##Lso=[-1,2,-1,-1,3];          % spatial order
##Lz=[-1,1,-1,-1,1];            % kinetic energy fix
##Lq=[-1,2,2,2,2];             % artificial viscosity
##Ldpi=[-1,0,-1,-1,0];         % dPI vs delta PI bar

%%%%%%%%%%%%%%%%%%%%%%%%  Presentation mi-stage
Lname=['vNR, q o1 vNR';
       'Godunov Sol acous o1';
       'BBC JCP 2009, q o1 vNR';
       'RK1, q o1 vNR, so 2, kin fix CRAS 2016'; 
       'RK5 Cash-Karp, q o1 vNR, so 5, kin fix CRAS 2016';];

Lsch=[300,002,200,100,104]; % scheme
Lso=[-1,-1,-1,2,3];         % spatial order
Lz=[-1,-1,-1,1,1];          % kinetic energy fix
Lq=[2,-1,2,2,2];            % artificial viscosity
Ldpi=[-1,-1,-1,0,0];        % dPI vs delta PI bar


ls=length(Lsch);


% Coeff de pseudo
Cq=1.5;
Cl=0.15;

% affichage
aff=1;

% Max iter
     Nmax=3E8;
     
% symoble pour la legende
Symleg=['.-'; '-.'; ':';'--', '*-'; 'x-'; '+-'; 'o-'; 'd-'; 's-'];


Etot=zeros(ls, 1e5);
     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:ls
  ["a.exe  -a ",num2str(aff)," -q ",num2str(Lq(j))," -d ",num2str(Cq)," -l ",num2str(Cl)," -p ",num2str(Ldpi(j))," -o ",num2str(Lso(j))," -z ",num2str(Lz(j))," -c ",num2str(Nmax)," -n ",num2str(nx)," -t ",num2str(tst)," -s ",num2str(Lsch(j))]
  system(["a.exe  -a ",num2str(aff)," -q ",num2str(Lq(j))," -d ",num2str(Cq)," -l ",num2str(Cl)," -p ",num2str(Ldpi(j))," -o ",num2str(Lso(j))," -z ",num2str(Lz(j))," -c ",num2str(Nmax)," -n ",num2str(nx)," -t ",num2str(tst)," -s ",num2str(Lsch(j))]);
  
  resu=load("sol.txt");
  [~,n]=size(resu);
  size(resu)
  Xd(j,:) = resu(1,:); % X décentrée (frontières)
  Xc(j,:) = resu(2,1:n-1);
  tau(j,:) = resu(3,1:n-1);
  u(j,:) = resu(4,:);
  e(j,:) = resu(5,1:n-1);
  p(j,:) = resu(6,1:n-1);
  epsilon(j,:) = resu(7,1:n-1);
  
  % energie du sytème
  NRJ=load("NRJ.txt");
  Nbiter(j)=length(NRJ(1,:));
  Etot(j,1:Nbiter(j))=NRJ(1,:);
  IMPUL(j,1:Nbiter(j))=NRJ(2,:);
  Mtot(j,1:Nbiter(j))=NRJ(3,:);
  VOL(j,1:Nbiter(j))=NRJ(4,:);
  
endfor



%% Cas du test de Sod
if tst==0 
  Nana=1e4;
  Xdana=a:(b-a)/Nana:b;
  for i=1:Nana
    Xcana(i)=(Xdana(i)+Xdana(i+1))/2;
  endfor 
  % sol analytique
  [W]=sol_ana_PU(Xcana,xdis,T,wL,wR,p12,u12);
  RHOana=W(1,:);
  Pana=W(3,:);
  [W]=sol_ana_PU(Xdana,xdis,T,wL,wR,p12,u12);
  Uana=W(2,:);

  figure
  hold on
  plot(Xcana,Pana./((gamma-1)*RHOana),'-');
  for j=1:ls
    plot(Xc(j,:),epsilon(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne spécifique ,  t=0.2s , nx=',num2str(nx),'']);
  xlabel("m");
  ylabel(["E  (J/kg)"]);
  if ls==1
    legend('ana',Lname(1,:),'location','northwest');
  elseif ls==2
    legend('ana',Lname(1,:),Lname(2,:),'location','northwest');
  elseif ls==3
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),'location','northwest');
  elseif ls==4
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','northwest');
  elseif ls==5
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','northwest');
  elseif ls==7
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),Lname(7,:),'location','northwest');
  endif
  hold off

##  figure
##  hold on
##  plot(Xcana,RHOana,'-','linewidth',1.5);
##  for j=1:ls
##    plot(Xc(j,:),1./tau(j,:),Symleg(j,:));
##  endfor
##  grid();
##  set(gca, "fontsize", 40);
##  title(['\rho nx=',num2str(nx)]);
##  if ls==1
##    legend('ana',Lname(1,:),'location','southwest');
##  elseif ls==2
##    legend('ana',Lname(1,:),Lname(2,:),'location','southwest');
##  elseif ls==3
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
##  elseif ls==4
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
##  elseif ls==5
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
##  elseif ls==6
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
##  endif
##  hold off
##
##  figure
##  hold on
##  plot(Xdana,Uana,'-');
##  for j=1:ls
##    plot(Xd(j,:),u(j,:),Symleg(j,:));
##  endfor
##  grid();
##  set(gca, "fontsize", 40);
##  title(['u nx=',num2str(nx)]);
##  if ls==1
##    legend('ana',Lname(1,:),'location','northwest');
##  elseif ls==2
##    legend('ana',Lname(1,:),Lname(2,:),'location','northwest');
##  elseif ls==3
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),'location','northwest');
##  elseif ls==4
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','northwest');
##  elseif ls==5
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
##  elseif ls==6
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
##  endif
##  hold off
##
##  figure
##  hold on
##  plot(Xcana,Pana,'-');
##  for j=1:ls
##    plot(Xc(j,:),p(j,:),Symleg(j,:));
##  endfor
##  grid();
##  set(gca, "fontsize", 40);
##  title(['p nx=',num2str(nx)]);
##  legend('ana',Lname(1,:));
##  if ls==1
##    legend('ana',Lname(1,:),'location','southwest');
##  elseif ls==2
##    legend('ana',Lname(1,:),Lname(2,:),'location','southwest');
##  elseif ls==3
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
##  elseif ls==4
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
##  elseif ls==5
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
##  elseif ls==6
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
##  endif
##  hold off
  

##  % Plot de l'energie du système au cours du temps
##  figure
##  hold on
##  for j=1:ls
##    t=0:T/(Nbiter(j)-1):T;
##    plot(t,Etot(j,1:Nbiter(j)),'.-');
##  end
##  grid();
##  set(gca, "fontsize", 40);
##  title(['energie totale Sod']);
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
##  endif
##  hold off
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
##  title(['impulsion Sod']);
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
##  endif
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
##  title(['masse totale Sod']);
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
##  endif
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
##  title(['volume total Sod']);
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
##  endif
##  hold off
  
  
end



%%  Cas du Bizarrium
if tst==2 || tst==5 || tst==6

  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),epsilon(j,:),Symleg(j,:));
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['Energie interne spécfique ,  t=',num2str(T),' s , nx=',num2str(nx)]);
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
  endif  
  hold off

##  figure
##  hold on
##  for j=1:ls
##    plot(Xc(j,:),tau(j,:),Symleg(j,:));
##  endfor
##  grid();
##  set(gca, "fontsize", 40);
##  title(['\tau nx=',num2str(nx)]);
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
##  endif  
##  hold off
##  
##  figure
##  hold on
##  for j=1:ls
##    plot(Xc(j,:),1./tau(j,:),Symleg(j,:));
##  endfor
##  grid();
##  set(gca, "fontsize", 40);
##  title(['\rho nx=',num2str(nx)]);
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
##  endif  
##  hold off
##
##  figure
##  hold on
##  for j=1:ls
##    plot(Xd(j,:),u(j,:),Symleg(j,:));
##  endfor
##  grid();
##  set(gca, "fontsize", 40);
##  title(['u nx=',num2str(nx)]);
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
##  endif  
##  hold off
##
##  figure
##  hold on
##  for j=1:ls
##    plot(Xc(j,:),p(j,:),Symleg(j,:));
##  endfor
##  grid();
##  set(gca, "fontsize", 40);
##  title(['p nx=',num2str(nx)]);
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
##  endif  
##  hold off
  

##  % Plot de l'energie du système au cours du temps
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
##    legend(Lname(1,:),'location','southeast');
##  elseif ls==2
##    legend(Lname(1,:),Lname(2,:),'location','southeast');
##  elseif ls==3
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southeast');
##  elseif ls==4
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southeast');
##  elseif ls==5
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
##  elseif ls==6
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
##  endif
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
##  title(['impulsion Bizarrium']);
##  if ls==1
##    legend(Lname(1,:),'location','southwest');
##  elseif ls==2
##    legend(Lname(1,:),Lname(2,:),'location','southwest');
##  elseif ls==3
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
##  elseif ls==4
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
##  elseif ls==5
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
##  elseif ls==6
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
##  endif
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
##  title(['masse totale Bizarrium']);
##  if ls==1
##    legend(Lname(1,:),'location','southwest');
##  elseif ls==2
##    legend(Lname(1,:),Lname(2,:),'location','southwest');
##  elseif ls==3
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
##  elseif ls==4
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
##  elseif ls==5
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
##  elseif ls==6
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
##  endif
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
##  title(['volume total Bizarrium']);
##  if ls==1
##    legend(Lname(1,:),'location','southwest');
##  elseif ls==2
##    legend(Lname(1,:),Lname(2,:),'location','southwest');
##  elseif ls==3
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
##  elseif ls==4
##    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
##  elseif ls==5
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
##  elseif ls==6
##    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
##  endif
##  hold off
##  
  
end

%%  Acoustic wave
if tst==3

  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),epsilon(j,:),'.-');
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\epsilon nx=',num2str(nx)]);
  if ls==1
    legend(Lname(1,:));
  elseif ls==2
    legend(Lname(1,:),Lname(2,:));
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  endif  
  hold off

  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),tau(j,:),'.-');
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
  endif  
  hold off

  figure
  hold on
  for j=1:ls
    plot(Xd(j,:),u(j,:),'.-');
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
  endif  
  hold off

  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),p(j,:),'.-');
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
  endif  
  hold off
  

   % Plot de l'energie du système au cours du temps
  figure
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,Etot(j,1:Nbiter(j)),'.-');
  end
  grid();
  set(gca, "fontsize", 40);
  title(['energie totale Sod']);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  endif
  hold off
  
  % Plot de l'impulsion au cours du temps
  figure
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,IMPUL(j,1:Nbiter(j)),'.-');
  end
  grid();
  set(gca, "fontsize", 40);
  title(['impulsion Sod']);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  endif
  hold off
  
  % Plot de la masse totale au cours du temps
  figure
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,Mtot(j,1:Nbiter(j)),'.-');
  end
  grid();
  set(gca, "fontsize", 40);
  title(['masse totale Sod']);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  endif
  hold off
  
  % Plot du volume total au cours du temps
  figure
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,VOL(j,1:Nbiter(j)),'.-');
  end
  grid();
  set(gca, "fontsize", 40);
  title(['volume total Sod']);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  endif
  hold off
  
  
end

%% Cas 3 etats
if tst==4
    figure
  hold on
  for j=1:ls
    plot(Xc(j,:),epsilon(j,:),'.-');
  endfor
  grid();
  set(gca, "fontsize", 40);
  title(['\epsilon nx=',num2str(nx)]);
  if ls==1
    legend(Lname(1,:));
  elseif ls==2
    legend(Lname(1,:),Lname(2,:));
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  endif  
  hold off

  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),tau(j,:),'.-');
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
  endif  
  hold off

  figure
  hold on
  for j=1:ls
    plot(Xd(j,:),u(j,:),'.-');
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
  endif  
  hold off

  figure
  hold on
  for j=1:ls
    plot(Xc(j,:),p(j,:),'.-');
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
  endif  
  hold off
  
  

  % Plot de l'energie du système au cours du temps
  figure
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,Etot(j,1:Nbiter(j)),'.-');
  end
  grid();
  set(gca, "fontsize", 40);
  title(['energie totale cas 3 etats']);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  endif
  hold off
  
  % Plot de l'impulsion au cours du temps
  figure
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,IMPUL(j,1:Nbiter(j)),'.-');
  end
  grid();
  set(gca, "fontsize", 40);
  title(['impulsion 3 etats']);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  endif
  hold off
  
  % Plot de la masse totale au cours du temps
  figure
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,Mtot(j,1:Nbiter(j)),'.-');
  end
  grid();
  set(gca, "fontsize", 40);
  title(['masse totale 3 etats']);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  endif
  hold off
  
  % Plot du volume total au cours du temps
  figure
  hold on
  for j=1:ls
    t=0:T/(Nbiter(j)-1):T;
    plot(t,VOL(j,1:Nbiter(j)),'.-');
  end
  grid();
  set(gca, "fontsize", 40);
  title(['volume total 3 etats']);
  if ls==1
    legend(Lname(1,:),'location','southwest');
  elseif ls==2
    legend(Lname(1,:),Lname(2,:),'location','southwest');
  elseif ls==3
    legend(Lname(1,:),Lname(2,:),Lname(3,:),'location','southwest');
  elseif ls==4
    legend(Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),'location','southwest');
  elseif ls==5
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),'location','southwest');
  elseif ls==6
    legend('ana',Lname(1,:),Lname(2,:),Lname(3,:),Lname(4,:),Lname(5,:),Lname(6,:),'location','southwest');
  endif
  hold off
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save('RKint_Biz_RK3_workspace_nx=16000_1105.mat');