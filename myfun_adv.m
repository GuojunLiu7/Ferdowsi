function F = myfun_adv(x)
global param_myfun_adv;
C_a    = param_myfun_adv(1);
R0_a   = param_myfun_adv(2);
DeltaV = param_myfun_adv(3);
IHR    = param_myfun_adv(4);
tau    = param_myfun_adv(5);
V_H    = param_myfun_adv(6);
Betta_H= param_myfun_adv(7);
P_init = param_myfun_adv(8);
HR_init= param_myfun_adv(9);
Alpha  = param_myfun_adv(10);
gamma  = param_myfun_adv(11);
Delta_h= param_myfun_adv(12);
Pfd    = param_myfun_adv(13);
Hfd    = param_myfun_adv(14);

fcnsig=@(x,y,z)1/(1+exp(-1*x*(z-y)));
% F=[fcnsig(x(1),x(2),Pfd)+(1/Alpha)*(Pfd/(Hfd*DeltaV*R0_a) - 1 - Alpha);
%     fcnsig(x(1),x(2),Pfd) - (Delta_h*IHR + Betta_H - Delta_h*Hfd)/(V_H+Betta_H)];
F=(fcnsig(x(1),x(2),Pfd)+(1/Alpha)*(Pfd/(Hfd*DeltaV*R0_a) - 1 - Alpha))^2+...
    (fcnsig(x(1),x(2),Pfd) - (Delta_h*IHR + Betta_H - Delta_h*Hfd)/(V_H+Betta_H))^2;
%**************************************************************************
%% 1-d equation; Iterative k & c
% F = abs(x(1) - (-1/(Alpha*Pfd*R0_a ))*( Pfd/(DeltaV*Hfd)+ x(2)* Alpha*R0_a - R0_a*(1+Alpha)))+ ...
%     abs(x(2)- (1/(V_H+Betta_H))*(Betta_H-Delta_h*Hfd+Delta_h*IHR - x(1)*(V_H*Pfd+Betta_H*Pfd)));

%% 2-d equation; Iterative k & c
% F = [x(1) - (-1/(Alpha*Pfd*R0_a ))*( Pfd/(DeltaV*Hfd)+ x(2)* Alpha*R0_a - R0_a*(1+Alpha));
%      x(2)- (1/(V_H+Betta_H))*(Betta_H-Delta_h*Hfd+Delta_h*IHR - x(1)*(V_H*Pfd+Betta_H*Pfd))];

%% 1-d equation; (Hf - Hfd)^2+(Pf - Pfd)^2
% input3=[C_a,R0_a,DeltaV,IHR,tau,V_H,Betta_H,P_init,HR_init,Alpha,gamma,Delta_h,x(2),x(1)];
% output = Fcn_Fsolve_Sol(input3);
% Pf  = cell2mat(output(1));
% Hf  = cell2mat(output(2));
% % outlyap=Fcn_Lyap_EigVal([input3,Pf]);
% % rleigvec=cell2mat(outlyap(1));
% % imeigvec=cell2mat(outlyap(2));
% % F = ((Hf - Hfd)/Hfd)^2+ ( (Pf - Pfd)/Pfd)^2;
% 
% F = (((-1/Delta_h)*( V_H * fcnsig(x(1),x(2),Pf) - Betta_H * (1-fcnsig(x(1),x(2),Pf)) - Delta_h * IHR ) - Hfd)/Hfd)^2+...
%     (( (1 + Alpha*(1-fcnsig(x(1),x(2),Pf))) * R0_a * DeltaV * Hf - Pfd)/Pfd)^2;

end
