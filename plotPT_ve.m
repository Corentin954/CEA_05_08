clear;

ZONE=load("MAPS.txt");
Vzone=load("Vmaps.txt");
Ezone=load("Emaps.txt");
[~,ne]=size(Ezone);
[~,nv]=size(Vzone);

P=ZONE(:,1);
T=ZONE(:,2);

P=reshape(P,ne,nv);
T=reshape(T,ne,nv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   P(V,E)
figure
hold on
%surf(Vzone*1e6,Ezone,P);
contour(Vzone*1e6,Ezone,P,50);
%pcolor(Vzone*1e6,Ezone,ZONE');
title("P(V,E)");
xlabel("V (cm^3/kg)");
ylabel("E (J/kg)");
grid()
set(gca, "fontsize", 40);
hold off


##%   T(V,E)
##figure
##hold on
##surf(Vzone*1e6,Ezone,T);
##%pcolor(Vzone*1e6,Ezone,ZONE');
##title("T(V,E)");
##xlabel("V (cm^3/kg)");
##ylabel("E (J/kg)");
##grid()
##set(gca, "fontsize", 40);
##hold off




