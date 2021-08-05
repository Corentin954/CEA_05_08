
% Solution analytique pour une 1-onde de détente puis une 3-onde de choc 
% x : points de calculs
% xdis : position de la discontinuité
% T : temps de la sol finale
% wL, wR : donnée initiale
function [W]=sol_ana_PU(x,xdis,T,wL,wR,p12,u12)

format long % permet d'avoir une dizaine de chiffre significatif


gamma=1.4;
mu2 = (gamma-1)/(gamma+1);


rhoL=wL(1);
uL=wL(2); 
pL=wL(3);

rhoR=wR(1);
uR=wR(2); 
pR=wR(3);

cL=sqrt(gamma*pL/rhoL); cR=sqrt(gamma*pR/rhoR);


%    PRE STAGE
% Etats 1 et 2
##rho1=rhoL*(p12/pL)^(1/gamma);
##rho2=rhoR*(p12+mu2*pR)/(pR+mu2*p12);
##c1=sqrt(gamma*p12/rho1);

% vitesse de l'onde de choc
##sigma = uR + sqrt( (p12+mu2*pR)/((1-mu2)*rhoR) ); % PRE-STAGE


%   LIVRE
% vitesse de l'onde de choc   
aR=sqrt( rhoR*((gamma+1)*p12/2 + (gamma-1)*pR/2) );% LIVRE
sigma = uR + aR/rhoR;

% Etat 1 et 2 provenant du livre p119
c1=cL+(gamma-1)*(uL-u12)/2;
rho1=gamma*(p12/c1^2);
rho2=rhoR*aR/(aR+rhoR*(uR-u12));

w1=[rho1; u12; p12];
w2=[rho2; u12; p12];


% Calcul des frontieres pour epsilon
eps1=uL-cL;
eps2=u12-c1;
eps3=u12;
eps4=sigma;

% Fonction pour la 1-onde de detente
udet=@(eps) (gamma-1)/(gamma+1)*uL + 2/(gamma+1)*(cL+eps) ;
cdet=@(eps) (gamma-1)/(gamma+1)*(uL-eps) + 2/(gamma+1)*cL ;
pdet=@(eps) pL*(cdet(eps)/cL).^(2*gamma/(gamma-1));
rhodet=@(eps) gamma*pdet(eps)./cdet(eps).^2;

wdet=@(eps) [rhodet(eps); udet(eps); pdet(eps) ];

% Solution finale
w=@(eps) wL*(eps<=eps1) + wdet(eps)*(eps1<eps).*(eps<=eps2) + w1*(eps2<eps).*(eps<=eps3) + w2*(eps3<eps).*(eps<=eps4) + wR*(eps4<=eps);


%  T le temps auquel on veut la sol
eps=(x-xdis)/T;

for i=1:length(eps)
  W(:,i)=w(eps(i));
end


% Plot solution finale 
##
##figure
##subplot(2,2,1)
##hold on
##plot(x,W(1,:));
##plot([eps1*T eps1*T]+xdis,[0 max(W(1,:))]);
##plot([eps2*T eps2*T]+xdis,[0 max(W(1,:))]);
##plot([eps3*T eps3*T]+xdis,[0 max(W(1,:))]);
##plot([eps4*T eps4*T]+xdis,[0 max(W(1,:))]);
##title('densite \rho')
##grid()
##hold off
##
##subplot(2,2,2)
##hold on
##plot(x,W(3,:));
##title('pression p')
##grid()
##hold off
##
##subplot(2,2,3)
##hold on
##plot(x,W(2,:));
##title('vitesse u')
##grid()
##hold off

endfunction