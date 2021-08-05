clear


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


% Diagramme P,T
SA_PT=load("fichiers/pSA_PT.txt");
SB_PT=load("fichiers/pSB_PT.txt");
SC_PT=load("fichiers/pSC_PT.txt");



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             POLYGONES  
figure
hold on
% Un point 
plot([0.000137230684781]*1e6 ,[6.125506226728779],'marker','+', "linewidth",5)
V= 0.000158049; E= 671163;
plot([V]*1e6 ,[E],'marker','+', "linewidth",5)

% Triangle
plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'r');
% Ploynomes
S=SA;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
S=SB;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
S=SC;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');

xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
grid()
set(gca, "fontsize", 40)
title("Diagramme (V,E)");
hold off

PT= 4510421781.609182357788085*1e-9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             POLYGONES  
figure
hold on
% Un point 
plot([10100]*1e-9 ,[300],'marker','+', "linewidth",5)
plot([4.2494e+10]*1e-9 ,[1258.17],'marker','+', "linewidth",5)
plot([3.94641e+10]*1e-9 ,[2435.59],'marker','+', "linewidth",5)
% Courbes mesures
% BETA/GAMMA
Pf=14;
Pplot=PT:(Pf - PT)/100:Pf;
Tplot=621.36 - 4.88*Pplot - 3.14*Pplot.*Pplot;
plot(Pplot, Tplot)

% BETA/liquide
Pf=0;
Pplot=PT:(Pf - PT)/100:Pf;
Tplot=505 + 32.1*Pplot - 1.77*Pplot.*Pplot;
plot(Pplot, Tplot)

% GAMMA/LIQUIDE
Pf=20;
Pplot=PT:(Pf - PT)/100:Pf;
Tplot=581.5 + 59.1*(Pplot - 2.87);
plot(Pplot, Tplot)
% ploynomes
S=SA_PT;
plot(S(:,1)*1e-9,S(:,2), "linewidth",0.5,"color",'black','marker','.');
S=SB_PT;
plot(S(:,1)*1e-9,S(:,2), "linewidth",0.5,"color",'black','marker','.');
S=SC_PT;
plot(S(:,1)*1e-9,S(:,2), "linewidth",0.5,"color",'black','marker','.');

xlabel("P (GPa)");
ylabel("T (K)");
grid()
set(gca, "fontsize", 40)
title("Diagramme (P,T)");
hold off



