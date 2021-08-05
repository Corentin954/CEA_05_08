clear, close al;

E=load("E.txt");
V=load("V.txt");

P1=load("P1.txt"); P2=load("P2.txt"); P3=load("P3.txt");
T1=load("T1.txt"); T2=load("T2.txt"); T3=load("T3.txt");
S1=load("S1.txt"); S2=load("S2.txt"); S3=load("S3.txt");
G1=load("G1.txt"); G2=load("G2.txt"); G3=load("G3.txt");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   CRITERE
CRITERE=abs(G1-G3) + abs(G1-G2);

figure
hold on
surf(E*1e-3,V*1e6,CRITERE);
title("CRITERE=abs(G1-G3) + abs(G1-G2)");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off

figure
hold on
surf(E,V*1e6,G1);
surf(E,V*1e6,G2);
surf(E,V*1e6,G3);
title("G1, G2, G3");
ylabel("V (cm^3/kg)");
xlabel("E (J/kg)");
grid()
set(gca, "fontsize", 40);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   P PHASE 1
figure
hold on
surf(E*1e-3,V*1e6,P1*1e-9);
title("P phase 1 (GPa)");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off

%   P PHASE 2
figure
hold on
surf(E*1e-3,V*1e6,P2*1e-9);
title("P phase 2 (GPa)");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off

%   P PHASE 3
figure
hold on
surf(E*1e-3,V*1e6,P3*1e-9);
title("P phase 3(GPa)");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T PHASE 1
figure
hold on
surf(E*1e-3,V*1e6,T1);
title("T phase 1(K)");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off

% T PHASE 2
figure
hold on
surf(E*1e-3,V*1e6,T2);
title("T phase 2(K)");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off

% T PHASE 3
figure
hold on
surf(E*1e-3,V*1e6,T3);
title("T phase 3(K)");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S PHASE 1
figure
hold on
surf(E*1e-3,V*1e6,S1);
title("S phase 1");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off

% S PHASE 2
figure
hold on
surf(E*1e-3,V*1e6,S2);
title("S phase 2");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off

% S PHASE 3
figure
hold on
surf(E*1e-3,V*1e6,S3);
title("S phase 3");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G PHASE 1
figure
hold on
surf(E*1e-3,V*1e6,G1);
title("G phase 1");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off

% G PHASE 2
figure
hold on
surf(E*1e-3,V*1e6,G2);
title("G phase 2");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off

% G PHASE 3
figure
hold on
surf(E*1e-3,V*1e6,G3);
title("G phase 3");
ylabel("V (cm^3/kg)");
xlabel("E (kJ/kg)");
grid()
set(gca, "fontsize", 40);
hold off

