clear, close all;

Ea12=load("Ea12.txt"); Eb12=load("Eb12.txt");
Va12=load("Va12.txt"); Vb12=load("Vb12.txt");
P12=load("P12.txt");   T12=load("T12.txt");
Ea12=load("Ea12.txt"); Va12=load("Va12.txt"); 
Eb12=load("Eb12.txt"); Vb12=load("Vb12.txt"); 

Ea13=load("Ea13.txt"); Eb13=load("Eb13.txt");
Va13=load("Va13.txt"); Eb13=load("Vb13.txt");
P13=load("P13.txt"); T13=load("T13.txt");
Ea13=load("Ea13.txt"); Va13=load("Va13.txt"); 
Eb13=load("Eb13.txt"); Vb13=load("Vb13.txt"); 

Ea23=load("Ea23.txt"); Eb23=load("Eb23.txt");
Va23=load("Va23.txt"); Vb23=load("Vb23.txt");
P23=load("P23.txt");   T23=load("T23.txt");
Ea23=load("Ea23.txt"); Va23=load("Va23.txt"); 
Eb23=load("Eb23.txt"); Vb23=load("Vb23.txt"); 


% Maillage
ptsAPT=load("ptsAPT.txt");
ptsAVE=load("ptsAVE.txt");
[na,~]=size(ptsAPT);

ptsBPT=load("ptsBPT.txt");
ptsBVE=load("ptsBVE.txt");
[nb,~]=size(ptsBPT);

ptsCPT=load("ptsCPT.txt");
ptsCVE=load("ptsCVE.txt");
[nc,~]=size(ptsCPT);

maillA=load("maillA.txt");
[ma,~]=size(maillA);
maillB=load("maillB.txt");
[mb,~]=size(maillB);
maillC=load("maillC.txt");
[mc,~]=size(maillC);


% Frontières (P,T)
frontBAS=load("fbas.txt");
frontHAUT=load("fhaut.txt");
frontGAUCHE=load("fgauche.txt");
frontDROITE=load("fdroite.txt");

% Frontières (V,E)
frontVE1=load("VEfront1.txt");
frontVE2=load("VEfront2.txt");
frontVE3=load("VEfront3.txt");


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DIAGRAMME PHASE (P,T)
figure
hold on
plot(P12*1e-9,T12);
plot(P13*1e-9,T13);
plot(P23*1e-9,T23);
% frontières
plot(frontBAS(:,1)*1e-9, frontBAS(:,2), "linewidth",2);
plot(frontHAUT(:,1)*1e-9, frontHAUT(:,2), "linewidth",2);
plot(frontGAUCHE(:,1)*1e-9, frontGAUCHE(:,2), "linewidth",2);
plot(frontDROITE(:,1)*1e-9, frontDROITE(:,2), "linewidth",2);
% points
##n=na;
##for i=1:n
##  plot(ptsAPT(i,1)*1e-9,ptsAPT(i,2),'marker','+','color','black');
##end
##n=nb;
##for i=1:n
##  plot(ptsBPT(i,1)*1e-9,ptsBPT(i,2),'marker','+','color','black');
##end 
n=nc;
for i=1:n
  plot(ptsCPT(i,1)*1e-9,ptsCPT(i,2),'marker','+','color','black');
end 
% mailles
##m=ma;
##for i=1:m
##  a=maillA(i,1);
##  b=maillA(i,2);
##  c=maillA(i,3);
##  plot([ptsAPT(a,1), ptsAPT(b,1)]*1e-9,[ptsAPT(a,2), ptsAPT(b,2)],'--');
##  plot([ptsAPT(b,1), ptsAPT(c,1)]*1e-9,[ptsAPT(b,2), ptsAPT(c,2)],'--');
##  plot([ptsAPT(c,1), ptsAPT(a,1)]*1e-9,[ptsAPT(c,2), ptsAPT(a,2)],'--');
##end 
##m=mb;
##for i=1:m
##  a=maillB(i,1);
##  b=maillB(i,2);
##  c=maillB(i,3);
##  plot([ptsBPT(a,1), ptsBPT(b,1)]*1e-9,[ptsBPT(a,2), ptsBPT(b,2)],'--');
##  plot([ptsBPT(b,1), ptsBPT(c,1)]*1e-9,[ptsBPT(b,2), ptsBPT(c,2)],'--');
##  plot([ptsBPT(c,1), ptsBPT(a,1)]*1e-9,[ptsBPT(c,2), ptsBPT(a,2)],'--');
##end 
m=mc;
for i=1:m
  a=maillC(i,1);
  b=maillC(i,2);
  c=maillC(i,3);
  plot([ptsCPT(a,1), ptsCPT(b,1)]*1e-9,[ptsCPT(a,2), ptsCPT(b,2)],'--');
  plot([ptsCPT(b,1), ptsCPT(c,1)]*1e-9,[ptsCPT(b,2), ptsCPT(c,2)],'--');
  plot([ptsCPT(c,1), ptsCPT(a,1)]*1e-9,[ptsCPT(c,2), ptsCPT(a,2)],'--');
end 
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
plot(Va12*1e6,Ea12);
plot(Vb12*1e6,Eb12);
% 13
plot(Va13*1e6,Ea13);
plot(Vb13*1e6,Eb13);
% 23
plot(Va23*1e6,Ea23);
plot(Vb23*1e6,Eb23);
% TRIANGLE
plot([Va12(1) Vb12(1)]*1e6,[Ea12(1) Eb12(1)],"linewidth",2,"color","b"); 
plot([Va13(1) Vb13(1)]*1e6,[Ea13(1) Eb13(1)],"linewidth",2,"color","b"); 
plot([Va23(1) Vb23(1)]*1e6,[Ea23(1) Eb23(1)],"linewidth",2,"color","b"); 
% Frontieres
plot(frontVE1(:,1)*1e6, frontVE1(:,2));
plot(frontVE2(:,1)*1e6, frontVE2(:,2));
plot(frontVE3(:,1)*1e6, frontVE3(:,2));
% POINTS
##n=na;
##for i=1:n
##  plot(ptsAVE(i,1)*1e6,ptsAVE(i,2),'marker','+','color','black');
##end 
##n=nb;
##for i=1:n
##  plot(ptsBVE(i,1)*1e6,ptsBVE(i,2),'marker','+','color','black');
##end 
n=nc;
for i=1:n
  plot(ptsCVE(i,1)*1e6,ptsCVE(i,2),'marker','+','color','black');
end 
% mailles
##m=ma;
##for i=1:m
##  a=maillA(i,1);
##  b=maillA(i,2);
##  c=maillA(i,3);
##  plot([ptsAVE(a,1), ptsAVE(b,1)]*1e6,[ptsAVE(a,2), ptsAVE(b,2)],'--');
##  plot([ptsAVE(b,1), ptsAPVE(c,1)]*1e6,[ptsAVE(b,2), ptsAVE(c,2)],'--');
##  plot([ptsAVE(c,1), ptsAVE(a,1)]*1e6,[ptsAVE(c,2), ptsAVE(a,2)],'--');
##end
##m=mb;
##for i=1:m
##  a=maillB(i,1);
##  b=maillB(i,2);
##  c=maillB(i,3);
##  plot([ptsBVE(a,1), ptsBVE(b,1)]*1e6,[ptsBVE(a,2), ptsBVE(b,2)],'--');
##  plot([ptsBVE(b,1), ptsBVE(c,1)]*1e6,[ptsBVE(b,2), ptsBVE(c,2)],'--');
##  plot([ptsBVE(c,1), ptsBVE(a,1)]*1e6,[ptsBVE(c,2), ptsBVE(a,2)],'--');
##end
m=mc;
for i=1:m
  figure
  hold on
  for i=1:10
  a=maillC(i,1);
  b=maillC(i,2);
  c=maillC(i,3);
  plot([ptsCVE(a,1), ptsCVE(b,1)]*1e6,[ptsCVE(a,2), ptsCVE(b,2)],'--');
  plot([ptsCVE(b,1), ptsCVE(c,1)]*1e6,[ptsCVE(b,2), ptsCVE(c,2)],'--');
  plot([ptsCVE(c,1), ptsCVE(a,1)]*1e6,[ptsCVE(c,2), ptsCVE(a,2)],'--');
  end
  hold off
  
  
end 
title("DIAGRAMME PHASE (V,E)");
xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
%legend("phase A","phase B");
grid()
set(gca, "fontsize", 40);
hold off



