clear

% hugoniot
Vhug=load("fichiers/Vhug.txt");
Ehug=load("fichiers/Ehug.txt");
Phug=load("fichiers/Phug.txt");
Thug=load("fichiers/Thug.txt");
Dhug=load("fichiers/Dhug.txt");
Uhug=load("fichiers/Uhug.txt");

XAhug=load("fichiers/XAhug.txt");
XBhug=load("fichiers/XBhug.txt");
XChug=load("fichiers/XChug.txt");

% Première hugoniotVhug=load("fichiers/Vhug.txt");
Vhug1=load("fichiers/Vhug1.txt");
Ehug1=load("fichiers/Ehug1.txt");
Phug1=load("fichiers/Phug1.txt");
Thug1=load("fichiers/Thug1.txt");
Dhug1=load("fichiers/Dhug1.txt");
Uhug1=load("fichiers/Uhug1.txt");

XAhug1=load("fichiers/XAhug1.txt");
XBhug1=load("fichiers/XBhug1.txt");
XChug1=load("fichiers/XChug1.txt");


V0= 137.230419*1e-6; E0= -0.286078; P0= 0.000101*1e9; T0= 300.000000;
V1= 123.823506*1e-6; E1= 51442.887611; P1= 7.674026*1e9; T1= 404.145405; D1= 3283.180446; u1= 320.754782;
V2= 120.650923*1e-6; E2= 77064.866965; P2= 8.478099*1e9; T2= 354.878407;
V3= 117.394432*1e-6; E3= 112608.972466; P3= 11.353935*1e9; T3= 388.433657; D3= 3283.180446; u3= 474.567695;
V4= 97.581086*1e-6; E4= 889478.598190; P4= 44.867180*1e9; T4= 1965.687527; D4= 4616.324687; u4= 1333.772744;
V5= 94.259292*1e-6; E5= 1222618.419852; P5= 56.904095*1e9; T5= 2231.323266; D5= 4993.825841; u5= 1563.722707;

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
%             POLYGONES  V,E
figure
hold on
% Un point 
plot(V0*1e6,E0,'marker','+', "linewidth",5);
plot(V1*1e6,E1,'marker','+', "linewidth",5);
plot(V2*1e6,E2,'marker','+', "linewidth",5);
plot(V3*1e6,E3,'marker','+', "linewidth",5);
plot(V4*1e6,E4,'marker','+', "linewidth",5);
plot(V5*1e6,E5,'marker','+', "linewidth",5);
% isentrope
plot(Vhug*1e6,Ehug,"linewidth",1.5);
plot(Vhug1*1e6,Ehug1,"linewidth",1.5,':');
% triangle
plot([Va,Vb,Vc,Va]*1e6,[Ea,Eb,Ec,Ea],'r');
% ploynomes
S=SA;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
S=SB;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
S=SC;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
xlim([80 160]);
ylim([-5e3 2e6]);
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
plot(P0*1e-9,T0,'marker','+', "linewidth",5);
plot(P1*1e-9,T1,'marker','+', "linewidth",5);
plot(P2*1e-9,T2,'marker','+', "linewidth",5);
plot(P3*1e-9,T3,'marker','+', "linewidth",5);
plot(P4*1e-9,T4,'marker','+', "linewidth",5);
plot(P5*1e-9,T5,'marker','+', "linewidth",5);
%isentrope
plot(Phug*1e-9,Thug,"linewidth",1.5);
plot(Phug1*1e-9,Thug1,"linewidth",1.5);
% ploynomes
S=SA_PT;
plot(S(:,1)*1e-9,S(:,2), "linewidth",0.5,"color",'black');
S=SB_PT;
plot(S(:,1)*1e-9,S(:,2), "linewidth",0.5,"color",'black');
S=SC_PT;
plot(S(:,1)*1e-9,S(:,2), "linewidth",0.5,"color",'black');
xlim([-2 100]);
ylim([0 3000]);
xlabel("P (GPa)");
ylabel("T (K)");
grid()
set(gca, "fontsize", 40)
title("Diagramme (P,T)");
hold off





figure
hold on
% points
plot(V0*1e6,P0*1e-9,'marker','+', "linewidth",5);
plot(V1*1e6,P1*1e-9,'marker','+', "linewidth",5);
plot(V2*1e6,P2*1e-9,'marker','+', "linewidth",5);
plot(V3*1e6,P3*1e-9,'marker','+', "linewidth",5);
plot(V4*1e6,P4*1e-9,'marker','+', "linewidth",5);
plot(V5*1e6,P5*1e-9,'marker','+', "linewidth",5);

%VM=1./8971; PM=21.56*1e9;
%VH=1./9415; PH=26.12*1e9;
VM=1./10331.2; PM=57.21*1e9;
VH=1./10848; PH=64.15*1e9;
plot((VM)*1e6,(PM)*1e-9,'marker','+', "linewidth",5);
plot((VH)*1e6,(PH)*1e-9,'marker','+', "linewidth",5);
%isentrope
plot(Vhug*1e6,Phug*1e-9,"linewidth",1.5);
plot(Vhug1*1e6,Phug1*1e-9,"-.","linewidth",1);
vplot=V0:(V3-V0)/100:V3-1e-5;
plot(vplot*1e6,(P0+(P1-P0)*(vplot-V0)/(V1-V0))*1e-9,"--","linewidth",0.8);
vplot=V0:(V4-V0)/100:V4-1e-5;
plot(vplot*1e6,(P0+(PM-P0)*(vplot-V0)/(VM-V0))*1e-9,"--","linewidth",0.8);
% Polygones
##P=SA_PT(:,2);
##V=SA(:,1);
##plot(V*1e6,P*1e-9, "linewidth",0.5,"color",'black');
##P=SB_PT(:,2);
##V=SB(:,1);
##plot(V*1e6,P*1e-9, "linewidth",0.5,"color",'black');
##P=SC_PT(:,2);
##V=SC(:,1);
##plot(V*1e6,P*1e-9, "linewidth",0.5,"color",'black');
ylabel("P (GPa)");
xlabel("V (cm^3/kg)");
grid()
set(gca, "fontsize", 40)
title("Diagramme (P,V)");
hold off



figure
hold on
%point
plot(u1,D1,'marker','+',"linewidth",5);
plot(u3,D3,'marker','+',"linewidth",5);
plot(u4,D4,'marker','+',"linewidth",5);
plot(u5,D5,'marker','+',"linewidth",5);
%hugoniot
plot(Uhug,Dhug,"linewidth",1.5);
xlabel("u (m/s)");
ylabel("D (m/s)");
grid()
set(gca,"fontsize",40)
title("Diagramme (u,D)");
hold off



figure
hold on
% points
plot([V0 V0]*1e6,[0 1],"--","linewidth",0.5);
plot([V1 V1]*1e6,[0 1],"--","linewidth",0.5);
plot([V2 V2]*1e6,[0 1],"--","linewidth",0.5);
plot([V3 V3]*1e6,[0 1],"--","linewidth",0.5);
plot([V4 V4]*1e6,[0 1],"--","linewidth",0.5);
plot([V5 V5]*1e6,[0 1],"--","linewidth",0.5);
%isentrope
plot(Vhug*1e6,XAhug,"linewidth",1.5);
plot(Vhug*1e6,XBhug,"linewidth",1.5);
plot(Vhug*1e6,XChug,"linewidth",1.5);

xlabel("V (cm^3/kg)");
ylabel("X fraction massique ");
legend('BETA','GAMMA','LIQ');
grid()
set(gca, "fontsize", 40)
title("Fraction massique en focntion de V");
hold off
##
##% PHASE 1/A BETA
##figure
##hold on
##% points
##plot([V0 V0]*1e6,[0 1],"--","linewidth",0.5);
##plot([V1 V1]*1e6,[0 1],"--","linewidth",0.5);
##plot([V2 V2]*1e6,[0 1],"--","linewidth",0.5);
##plot([V3 V3]*1e6,[0 1],"--","linewidth",0.5);
##plot([V4 V4]*1e6,[0 1],"--","linewidth",0.5);
##plot([V5 V5]*1e6,[0 1],"--","linewidth",0.5);
##%isentrope
##plot(Vhug*1e6,XAhug,"linewidth",1.5);
##
##xlabel("V (cm^3/kg)");
##ylabel("X fraction massique ");
##grid()
##set(gca, "fontsize", 40)
##title("Fraction massique BETA (1/A)");
##hold off
##
##
##% PHASE 1/A BETA
##figure
##hold on
##% points
##plot([V0 V0]*1e6,[0 1],"--","linewidth",0.5);
##plot([V1 V1]*1e6,[0 1],"--","linewidth",0.5);
##plot([V2 V2]*1e6,[0 1],"--","linewidth",0.5);
##plot([V3 V3]*1e6,[0 1],"--","linewidth",0.5);
##plot([V4 V4]*1e6,[0 1],"--","linewidth",0.5);
##plot([V5 V5]*1e6,[0 1],"--","linewidth",0.5);
##%isentrope
##plot(Vhug*1e6,XBhug,"linewidth",1.5);
##
##xlabel("V (cm^3/kg)");
##ylabel("X fraction massique ");
##grid()
##set(gca, "fontsize", 40)
##title("Fraction massique GAMMA (2/B)");
##hold off
##
##
##% PHASE 3/C LIQUIDE
##figure
##hold on
##% points
##plot([V0 V0]*1e6,[0 1],"--","linewidth",0.5);
##plot([V1 V1]*1e6,[0 1],"--","linewidth",0.5);
##plot([V2 V2]*1e6,[0 1],"--","linewidth",0.5);
##plot([V3 V3]*1e6,[0 1],"--","linewidth",0.5);
##plot([V4 V4]*1e6,[0 1],"--","linewidth",0.5);
##plot([V5 V5]*1e6,[0 1],"--","linewidth",0.5);
##%isentrope
##plot(Vhug*1e6,XChug,"linewidth",1.5);
##
##xlabel("V (cm^3/kg)");
##ylabel("X fraction massique ");
##grid()
##set(gca, "fontsize", 40)
##title("Fraction massique LIQUIDE(3/C)");
##hold off
##
##
##
