% Solution exacte AVERAGE
% n : nx+1

function [ RHObar,Ubar,Pbar]=solbar(tst,m,Xd,Xu,n,xdis,T,wL,wR,p12,u12, fac_eps)
     
    format long % permet d'avoir une dizaine de chiffre significatif 
     
    ip=3; % simpson 3/8  NE PAS MODIF ip=3
    
    if tst==0  % Sod 
      gamma=1.4;
      mu2 = (gamma-1)/(gamma+1);

      rhoL=wL(1);  rhoR=wR(1);
      uL=wL(2);    uR=wR(2); 
      pL=wL(3);    pR=wR(3);

      cL=sqrt(gamma*pL/rhoL); cR=sqrt(gamma*pR/rhoR);

      %   LIVRE
      % vitesse de l'onde de choc   
      aR=sqrt( rhoR*((gamma+1)*p12/2 + (gamma-1)*pR/2) );% LIVRE
      sigma = uR + aR/rhoR;

      % Etat 1 et 2 provenant du livre p119
      c1=cL+(gamma-1)*(uL-u12)/2;
      rho1=gamma*(p12/c1^2);
      rho2=rhoR*aR/(aR+rhoR*(uR-u12));

      % fonctions  sol
      w1=[rho1; u12; p12];   w2=[rho2; u12; p12];
     
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

      % Solution finale
      rhosol=@(eps) wL(1).*(eps<=eps1) + rhodet(eps).*(eps1<eps).*(eps<=eps2) + w1(1).*(eps2<eps).*(eps<=eps3) + w2(1).*(eps3<eps).*(eps<=eps4) + wR(1)*(eps4<=eps);
      usol=@(eps) wL(2).*(eps<=eps1) + udet(eps).*(eps1<eps).*(eps<=eps2) + w1(2).*(eps2<eps).*(eps<=eps3) + w2(2).*(eps3<eps).*(eps<=eps4) + wR(2).*(eps4<=eps);
      psol=@(eps) wL(3).*(eps<=eps1) + pdet(eps).*(eps1<eps).*(eps<=eps2) + w1(3).*(eps2<eps).*(eps<=eps3) + w2(3).*(eps3<eps).*(eps<=eps4) + wR(3).*(eps4<=eps);
   
   elseif tst==2  % Planar Kidder 
     kappa=2*pi;
     rho0=1.0; p0=5./7;
     
     rhosol=@(eps) rho0 + fac_eps*sin(kappa*T*(eps-1));
     usol=@(eps) fac_eps*sin(kappa*T*(eps-1));
     psol=@(eps) p0 + fac_eps*sin(kappa*T*(eps-1));
   end
    
    %%%%%%%%%%%%%%%%%%%%   rho et p : centrée aux mailles      %%%%%%%%%%%%%%%%%%%%%%
    % Points de calcul centrée aux mailles
    for i=1:n-1
      ind=(i-1)*m*ip+1;
      Xsol(ind)=Xd(i);
      Xtemp=Xd(i):(Xd(i+1)-Xd(i))/(ip*m):Xd(i+1);
      Xsol(ind+1:ind+ip*m-1)=Xtemp(2:ip*m);
    endfor
    ip
    m
    i
    Xsol((n-1)*ip*m+1)=Xd(n);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eps=(Xsol-xdis)/T;
    RHOana=rhosol(eps);
    Pana=psol(eps);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Methode d'intégration (simpsons 3/8) sur chaque maille divisé en m sous-intervalles
    
    %%  RHO
    RHObar=zeros(1,n-1); 
    for i=1:n-1 
      h=(Xd(i+1)-Xd(i))/(m);
      indm=(i-1)*m*ip+1;
      RHObar(i)=RHOana(indm); % f(a) ref maille i
      %% 1er maille
      for lk=1:ip-1
          RHObar(i)=RHObar(i)+3*RHOana(indm+lk);
      endfor
      %%
      for lj=indm+ip:ip:indm+ip*(m-1)
        RHObar(i)=RHObar(i)+2*RHOana(lj);
        for lk=1:ip-1
          RHObar(i)=RHObar(i)+3*RHOana(lj+lk);
        endfor
      endfor
      RHObar(i)=RHObar(i)+RHOana(i*m*ip+1); % f(b) ref maille
      RHObar(i)=(1/(8*m))*RHObar(i); % On ne veut pas l'intégrale mais la valeur moyennes
    endfor

    %%  P
    Pbar=zeros(1,n-1); 
    for i=1:n-1 
      h=(Xd(i+1)-Xd(i))/(m);
      indm=(i-1)*m*ip+1;
      Pbar(i)=Pana(indm); % f(a) ref maille i
      %% 1er maille
      for lk=1:ip-1
          Pbar(i)=Pbar(i)+3*Pana(indm+lk);
      endfor
      %%
      for lj=indm+ip:ip:indm+ip*(m-1)
        Pbar(i)=Pbar(i)+2*Pana(lj);
        for lk=1:ip-1
          Pbar(i)=Pbar(i)+3*Pana(lj+lk);
        endfor
      endfor
      Pbar(i)=Pbar(i)+Pana(i*m*ip+1); % f(b) ref maille
      Pbar(i)=(1/(8*m))*Pbar(i); % On ne veut pas l'intégrale mais la valeur moyennes
    endfor
    

    
%%%%%%%%%%%%%%%%%%%%   u : frontieres aux mailles      %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Points de calcul centrée aux mailles
    for i=1:n
      ind=(i-1)*m*ip+1;
      XsolU(ind)=Xu(i);
      Xtemp=Xu(i):(Xu(i+1)-Xu(i))/(ip*m):Xu(i+1);
      XsolU(ind+1:ind+ip*m-1)=Xtemp(2:ip*m);
    endfor
    XsolU(n*ip*m+1)=Xu(n+1);
    size(Xsol)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eps=(XsolU-xdis)/T;
    Uana=usol(eps);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Methode d'intégration (simpsons 3/8) sur chaque maille divisé en m sous-intervalles
    %%  U
    Ubar=zeros(1,n); 
    for i=1:n
      h=(Xu(i+1)-Xu(i))/(m);
      indm=(i-1)*m*ip+1;
      Ubar(i)=Uana(indm); % f(a) ref maille i
      %% 1er maille
      for lk=1:ip-1
          Ubar(i)=Ubar(i)+3*Uana(indm+lk);
      endfor
      %%
      for lj=indm+ip:ip:indm+ip*(m-1)
        Ubar(i)=Ubar(i)+2*Uana(lj);
        for lk=1:ip-1
          Ubar(i)=Ubar(i)+3*Uana(lj+lk);
        endfor
      endfor
      Ubar(i)=Ubar(i)+Uana(i*m*ip+1); % f(b) ref maille
      Ubar(i)=(1/(8*m))*Ubar(i); % On ne veut pas l'intégrale mais la valeur moyennes
    endfor
    
    
endfunction