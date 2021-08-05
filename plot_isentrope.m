clear

% isentrope
Visen=load("fichiers/Visen.txt");
Eisen=load("fichiers/Eisen.txt");
Pisen=load("fichiers/Pisen.txt");
Tisen=load("fichiers/Tisen.txt");
XAisen=load("fichiers/XAisen.txt");
XBisen=load("fichiers/XBisen.txt");
XCisen=load("fichiers/XCisen.txt");
Cisen=load("fichiers/Cisen.txt");

% hugoniot
Vhug=load("fichiers/Vhug.txt");
Ehug=load("fichiers/Ehug.txt");
Phug=load("fichiers/Phug.txt");
Thug=load("fichiers/Thug.txt");




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


[na,~]=size(SA); % de taille NA+1 car SA[NA]=S[0] (code C)
[nb,~]=size(SB);
[nc,~]=size(SC);
[nab,~]=size(SAB); % de taille NBA+1 car SA[NA]=S[0] (code C)
[nac,~]=size(SAC);
[nbc,~]=size(SBC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             POLYGONES  
figure
hold on
% Un point 
plot([0.000137230684781]*1e6 ,[6.125506226728779],'marker','+', "linewidth",5)
plot([0.000135077]*1e6 ,[-53826.5],'marker','+', "linewidth",5)
% isentrope
plot(Visen(:,1)*1e6,Eisen(:,1));
plot(Visen(:,2)*1e6,Eisen(:,2));
%hugoniot
plot(Vhug*1e6,Ehug);
% triangle
plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'r');
% ploynomes
n=na-1;
S=SA;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
n=nb-1;
S=SB;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
n=nc-1;
S=SC;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
##n=nab-1;
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
##n=nac-1;
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');
##n=nbc-1;
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black','marker','.');

xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
grid()
set(gca, "fontsize", 40)
title("Diagramme (V,E)");
hold off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             POLYGONES  
figure
hold on
% Un point 
plot([10100]*1e-9 ,[300],'marker','+', "linewidth",5)
plot([4.2494e+10]*1e-9 ,[1258.17],'marker','+', "linewidth",5)
plot([3.94641e+10]*1e-9 ,[2435.59],'marker','+', "linewidth",5)
%isentrope
plot(Pisen(:,1)*1e-9,Tisen(:,1));
plot(Pisen(:,2)*1e-9,Tisen(:,2));
%hugoniot
plot(Phug*1e-9,Thug);
% ploynomes
n=na-1;
S=SA_PT;
plot(S(:,1)*1e-9,S(:,2), "linewidth",0.5,"color",'black','marker','.');
n=nb-1;
S=SB_PT;
plot(S(:,1)*1e-9,S(:,2), "linewidth",0.5,"color",'black','marker','.');
n=nc-1;
S=SC_PT;
plot(S(:,1)*1e-9,S(:,2), "linewidth",0.5,"color",'black','marker','.');

xlabel("P (GPa)");
ylabel("T (K°)");
grid()
set(gca, "fontsize", 40)
title("Diagramme (P,T)");
hold off





figure
hold on
%isentrope
plot(Visen(:,1)*1e6,Pisen(:,1)*1e-9);
plot(Visen(:,2)*1e6,Pisen(:,2)*1e-9);
%hugoniot
plot(Vhug*1e6,Phug*1e-9);
ylabel("P (GPa)");
xlabel("V (cm^3/kg)");
grid()
set(gca, "fontsize", 40)
title("Diagramme (P,T)");
hold off


