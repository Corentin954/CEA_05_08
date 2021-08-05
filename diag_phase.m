clear;


% rectangle
 Vmin= 127.4949899800;
 Vmax= 127.7334669339;
 Emin= 97835.6713426854;
 Emax= 100951.9038076152;
 
 REC(1,1)=Vmin; REC(1,2)=Emin;
 REC(2,1)=Vmin; REC(2,2)=Emax;
 REC(3,1)=Vmax; REC(3,2)=Emax;
 REC(4,1)=Vmax; REC(4,2)=Emin;


Ea12=load("fichiers/Ea12.txt"); Eb12=load("fichiers/Eb12.txt");
Va12=load("fichiers/Va12.txt"); Vb12=load("fichiers/Vb12.txt");
P12=load("fichiers/P12.txt");   T12=load("fichiers/T12.txt");
Ea12=load("fichiers/Ea12.txt"); Va12=load("fichiers/Va12.txt"); 
Eb12=load("fichiers/Eb12.txt"); Vb12=load("fichiers/Vb12.txt"); 

Ea13=load("fichiers/Ea13.txt"); Eb13=load("fichiers/Eb13.txt");
Va13=load("fichiers/Va13.txt"); Eb13=load("fichiers/Vb13.txt");
P13=load("fichiers/P13.txt"); T13=load("fichiers/T13.txt");
Ea13=load("fichiers/Ea13.txt"); Va13=load("fichiers/Va13.txt"); 
Eb13=load("fichiers/Eb13.txt"); Vb13=load("fichiers/Vb13.txt"); 

Ea23=load("fichiers/Ea23.txt"); Eb23=load("fichiers/Eb23.txt");
Va23=load("fichiers/Va23.txt"); Vb23=load("fichiers/Vb23.txt");
P23=load("fichiers/P23.txt");   T23=load("fichiers/T23.txt");
Ea23=load("fichiers/Ea23.txt"); Va23=load("fichiers/Va23.txt"); 
Eb23=load("fichiers/Eb23.txt"); Vb23=load("fichiers/Vb23.txt"); 



% Frontières (P,T)
frontBAS=load("fichiers/fbas.txt");
frontHAUT=load("fichiers/fhaut.txt");
frontGAUCHE=load("fichiers/fgauche.txt");
frontDROITE=load("fichiers/fdroite.txt");

% Frontières (V,E)
frontVE1=load("fichiers/VEfront1.txt");
frontVE2=load("fichiers/VEfront2.txt");
frontVE3=load("fichiers/VEfront3.txt");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DIAGRAMME PHASE (P,T)
figure
hold on
plot(P12*1e-9,T12, "linewidth",2);
plot(P13*1e-9,T13, "linewidth",2);
plot(P23*1e-9,T23, "linewidth",2);
% frontières
plot(frontBAS(:,1)*1e-9, frontBAS(:,2), "linewidth",2);
plot(frontHAUT(:,1)*1e-9, frontHAUT(:,2), "linewidth",2);
plot(frontGAUCHE(:,1)*1e-9, frontGAUCHE(:,2), "linewidth",2);
plot(frontDROITE(:,1)*1e-9, frontDROITE(:,2), "linewidth",2);
title("DIAGRAMME PHASE (P,T)");
xlabel("P (GPa)");
ylabel("T (K)");
grid()
set(gca, "fontsize", 40);
hold off



%   DIAGRAMME PHASE (V,E)
figure
hold on
% 12
plot(Va12*1e6,Ea12, "linewidth",2);
plot(Vb12*1e6,Eb12, "linewidth",2);
% 13
plot(Va13*1e6,Ea13, "linewidth",2);
plot(Vb13*1e6,Eb13, "linewidth",2);
% 23
plot(Va23*1e6,Ea23, "linewidth",2);
plot(Vb23*1e6,Eb23, "linewidth",2);
% TRIANGLE
plot([Va12(1) Vb12(1)]*1e6,[Ea12(1) Eb12(1)],"linewidth",2,"color","b"); 
plot([Va13(1) Vb13(1)]*1e6,[Ea13(1) Eb13(1)],"linewidth",2,"color","b"); 
plot([Va23(1) Vb23(1)]*1e6,[Ea23(1) Eb23(1)],"linewidth",2,"color","b"); 
% Frontieres
plot(frontVE1(:,1)*1e6, frontVE1(:,2), "linewidth",2);
plot(frontVE2(:,1)*1e6, frontVE2(:,2), "linewidth",2);
plot(frontVE3(:,1)*1e6, frontVE3(:,2), "linewidth",2);
% Plot rectagnle maillage
plot([REC(1,1) REC(2,1)],[REC(1,2) REC(2,2)],'r')
plot([REC(2,1) REC(3,1)],[REC(2,2) REC(3,2)],'r')
plot([REC(3,1) REC(4,1)],[REC(3,2) REC(4,2)],'r')
plot([REC(4,1) REC(1,1)],[REC(4,2) REC(1,2)],'r')
title("DIAGRAMME PHASE (V,E)");
xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
%legend("phase A","phase B");
grid()
set(gca, "fontsize", 40);
hold off



