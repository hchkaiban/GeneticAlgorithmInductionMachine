clear all

Tsim=0.0001; %sample time
%Inv_nb=4;   %1:SPWM 2:DPWM 3:VPWM 4:SuperCos

%Moment d’inertie 
J = 0.0153;% kg.m2
%Nombre de paires de pôles 
p = 4;
%Inductance mutuelle cyclique 
Msr = 0.172; %H
%Inductance rotorique cyclique 
Lr = 0.188; %H
%Inductance statorique cyclique 
Ls = 0.188; %H
%Résistance rotorique 
Rr = 1.26; %Ohm
%Résistance statorique 
Rs = 1.34; %Ohm
%Coefficient de Blondel
sigma=1-(Msr^2/(Lr*Ls)); 
%Constante de temps rotorique
Tr=Lr/Rr;
%friction
f=0.0038; %Nm/rad/s
%Resistance equivalente des pertes fer (Cf 8)
Rfe=0.3;    %Ohm
%Inductance de fuite statorique (Cf8)
Lfs=0.01;
%Inductance de fuite rotorique (Cf8)
Lfr=0.008;

%Vitesse nominale
Vit_nom=150.79; %rad/s
%gamma Cf 11 p 65
gamma = (Rs+(Rr*Msr^2/Lr^2))/(sigma*Ls);
% Equation caracteristique vitesse boucle fermee: s^2+gamma*s+KF (Cf 11 p
% 66)
KF=Msr/(sigma*Ls*Tr);

%Observateur de Luenberger d'ordre 2 Cf 11 p76
%Phird_obs=0.05;
%A=[-f/J -p/J; 0 0];
%B=[(p/J)*(p*Msr*Phird_obs/Lr); 0];
%C=[1 0];
%D=0;
%sys_ss = ss(A, B, C, D);
%Plant = tf(sys_ss)

%Observabilite
%O=[C; C*A];
%rank(O);
%obs_pole = [-58; -50];
        
%[k_obs,prec,message] = place(A', C', obs_pole);
%l1 = k_obs * [1 0]';
%l2 = k_obs * [0 1]';

%Tmec = J/f;
% for AG initial estimation:
%K1 = -(Rs/(sigma*Ls*Tmec))
%K2 = (1/(sigma*Ls*Tmec))
%K3 = -Rs /(sigma*Ls)
%K4 = 1/(sigma*Ls)
%K5 = (-Rr*Ls - Rs*Lr)/(sigma*Ls*Lr)


