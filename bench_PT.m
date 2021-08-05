clear



P=load("fichiers/Ptest.txt");
T=load("fichiers/Ttest.txt");
X=load("fichiers/Xtest.txt");
V=load("fichiers/Vtest.txt");
E=load("fichiers/Etest.txt");
C=load("fichiers/Ctest.txt");
%G=load("fichiers/Gtest.txt");

ZONE_test=load("fichiers/ZONEtest.txt");


%%%%%%%%%%%%%%%%%%%%%%%%

[~,ne]=size(E);
[~,nv]=size(V);

P=reshape(P,ne,nv);
T=reshape(T,ne,nv);
X=reshape(X,ne,nv);
ZONE_test=reshape(ZONE_test,nv,ne);
C=reshape(C,nv,ne);
%G=reshape(G,nv,ne);

##Nb_ZONE=load("fichiers/Nb_ZONEbench.txt");
##V_dis=load("fichiers/Vbench_dis.txt");
##E_dis=load("fichiers/Ebench_dis.txt");
##[nv_dis,~]=size(V_dis);
##[ne_dis,~]=size(E_dis);
##Nb_ZONE=reshape(Nb_ZONE,nv_dis,ne_dis);


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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%          P(V,E)  
figure
hold on
%surf(V*1e6,E,P*1e-9);
%%contour(V*1e6,E,P,50);
%%pcolor(V*1e6,E,P);
%%contourf(V*1e6,E,P*1e-9,[0:2:70]);
contourf(V*1e6,E,P*1e-9,100);
% Un point 
plot([0.000137230684781]*1e6 ,[ 6.125506226728779],'marker','+', "linewidth",5)
%% ploynomes
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
grid()
set(gca, "fontsize", 40)
title("P(V,E)");
hold off


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##               T(V,E)
figure
hold on
%%surf(V*1e6,E,T);
%%contour(V*1e6,E,T,50);
%%pcolor(V*1e6,E,T);
%contourf(V*1e6,E,T,[300:50:2500]);
contourf(V*1e6,E,T,100);
%% ploynomes
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
grid()
set(gca, "fontsize", 40)
title("T(V,E)");
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%m�               C(V,E)
figure
hold on
surf(V*1e6,E,C);
%%contour(V*1e6,E,C,50);
%%pcolor(V*1e6,E,C);
%contourf(V*1e6,E,C,100);
% ploynomes
##n=nab-1;
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nac-1;
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nbc-1;
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
zlabel("vitesse du son (m/s)");
grid()
set(gca, "fontsize", 40)
title("C(V,E)");
hold off


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##%%%%m�               G(V,E)
##figure
##hold on
##surf(V*1e6,E,G);
##%%contour(V*1e6,E,G,50);
##%%pcolor(V*1e6,E,G);
##%contourf(V*1e6,E,G,100);
##% ploynomes
##n=nab-1;
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nac-1;
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nbc-1;
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##zlabel(" Derivee fondamentale");
##grid()
##set(gca, "fontsize", 40)
##title(" Derivee fondamentale G(V,E)");
##hold off


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####%               X fraction massique
figure
hold on
%surf(V*1e6,E,X);
%contour(V*1e6,E,X,50);
%pcolor(V*1e6,E,X);
contourf(V*1e6,E,X,[0:0.1:1]);
% ploynomes
S=SAB;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
S=SAC;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
S=SBC;
plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
grid()
set(gca, "fontsize", 40)
title("fraction massique dans les zones de m�lange");
hold off





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             Nombre de zones par mailles   
##figure
##hold on
##
##%%%% ploynomes
##n=na-1;
##S=SA;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nb-1;
##S=SB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nc-1;
##S=SC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nab-1;
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nac-1;
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nbc-1;
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');

%%% grille
##he=E_dis(2)-E_dis(1);
##Edeb=E_dis(1)-he/2;
##Efin=E_dis(ne_dis)+he/2;
##
##hv=V_dis(2)-V_dis(1);
##Vdeb=V_dis(1)-hv/2;
##Vfin=V_dis(nv_dis)+hv/2;
##
##for i=1:nv_dis
##  Vplot=Vdeb+i*hv;
##  plot([Vplot Vplot]*1e6,[Edeb Efin]);
##end
##for i=1:ne_dis+1
##  Eplot=Edeb+i*he;
##  plot([V_dis(1) V_dis(nv_dis)]*1e6,[Eplot Eplot]);
##end
##
##% plot
##surf(V_dis*1e6,E_dis,Nb_ZONE) %,'marker','o');
##%contour(V*1e6,E,P,50);
##%pcolor(V*1e6,E,P);
##%contourf(V*1e6,E,P*1e-9,[0:2:70]);
##%contourf(V*1e6,E,Nb_ZONE,100);
##
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##grid()
##set(gca, "fontsize", 40)
##title("Nb zone par case");
##hold off




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%             Zones par mailles   
##figure
##hold on
##surf(V*1e6,E,ZONE_test);
##contourf(V*1e6,E,Nb_ZONE,100);
##%%%%% ploynomes
##n=nab-1;
##S=SAB;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nac-1;
##S=SAC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##n=nbc-1;
##S=SBC;
##plot(S(:,1)*1e6,S(:,2), "linewidth",0.5,"color",'black');
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##grid()
##set(gca, "fontsize", 40)
##title("Nb zone par case");
##hold off


