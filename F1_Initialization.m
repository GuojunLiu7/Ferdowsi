function F1_Initialization(vargin)
global param;
bound.ul=cell2mat(vargin(1,1));
bound.ll=cell2mat(vargin(1,2));
param.title= [{'C_{a}'},{'R^{0}_{a}'},{'\DeltaV'},{'\tau'},{'V_H'},...
              {'\beta_{H}'},{'\alpha'},{'\gamma'},{'\delta_{H}'},{'H_0'},...
              {'P_{sp}'},{'\alpha_{sp}'}];
param.nom  = [1.55,    0.6, 50,  3, 1.17,  0.84,  1.3,  0.2,   1.7,  1.66,  110, 0.05];
param.lb=zeros(1,12);
param.ub=zeros(1,12);

param.ub(1)  = param.nom(1);            param.lb(1) = param.nom(1);
param.ub(2)  = param.nom(2);            param.lb(2) = param.nom(2);
param.ub(3)  = param.nom(3)*bound.ul;   param.lb(3) = param.nom(3)*bound.ll;
param.ub(4)  = param.nom(4);            param.lb(4) = param.nom(4);
param.ub(5)  = param.nom(5)*bound.ul;   param.lb(5) = param.nom(5)*bound.ll;
param.ub(6)  = param.nom(6)*bound.ul;   param.lb(6) = param.nom(6)*bound.ll;
param.ub(7)  = param.nom(7)*bound.ul;   param.lb(7) = param.nom(7)*bound.ll;
param.ub(8)  = param.nom(8);            param.lb(8) = param.nom(8);
param.ub(9)  = param.nom(9);            param.lb(9) = param.nom(9);
param.ub(10) = param.nom(10);           param.lb(10)= param.nom(10);
param.ub(11) = param.nom(11)*bound.ul;  param.lb(11)= param.nom(11)*bound.ll;
param.ub(12) = param.nom(12);           param.lb(12)= param.nom(12);

end