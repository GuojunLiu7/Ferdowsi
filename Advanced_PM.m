clear all; 
clc; 
close all;
global param_myfun_adv;
%**************************************************************************
% GUI to get input
%**************************************************************************
profile on
options.Resize      ='on';
options.WindowStyle ='modal';
options.Interpreter ='tex';
prompt = {'C_{a}','R^0_a','\DeltaV','IHR','\tau','V_H','\beta_H','P_{init}', ...
    'HR_{init}','\alpha','\gamma','\delta_{H}','Pf','Hf'};
dlg_title = 'Cardiovascular System Parameters';
num_lines= 1;
def  = {'1.55','0.6', '50','1.66','3','1.17','0.84','160','2','1.3','0.2','1.7','91','1.7'};
input1  = inputdlg(prompt,dlg_title,num_lines,def,options);

%**************************************************************************
input2=zeros(length(def),1);
for i=1:length(def)
    input2(i)=eval(cell2mat(input1(i)));
end
C_a       = input2(1);
R0_a      = input2(2);
DeltaV    = input2(3);
IHR       = input2(4);
tau       = input2(5);
V_H       = input2(6);
Beta_H    = input2(7);
P_init    = input2(8);
HR_init   = input2(9);
Alpha     = input2(10);
gamma     = input2(11);
Delta_h   = input2(12);
Pfd       = input2(13);
Hfd       = input2(14);

%**************************************************************************
% Builit-in MATLAB iterative function
%**************************************************************************
param_myfun_adv=input2;
x0 = [0.05;100];              % Make a starting guess at the solution
% options=optimset('Algorithm','active-set');
% options.NonlEqnAlgorithm='dogleg';
options=optimset('Display','iter');
options.MaxFunEvals=10^5;
options.MaxIter=10^5;
[x,fval] = fmincon(@myfun_adv,x0,[-1,0;0,-1],[0;0],[],[],[0.005;70],[0.2;180],[],options); % Call optimizer
% [x,fval,exitflag,output] = fsolve(@myfun_adv,x0,options); % Call optimizer
% global outfsolve;
Result=[input2;x];

%**************************************************************************
% a set of equation related to the advanced pacemaker
%**************************************************************************
% cte(1)=-1/(Alpha*DeltaV*R0_a);
% cte(2)=(1+Alpha)/Alpha;
% cte(3)=(-Delta_h)/(V_H+Beta_H);
% cte(4)=(Delta_h*IHR+Beta_H)/(V_H+Beta_H);
% 
% cte(5)=cte(3)*Hfd+cte(4);
% cte(6)=cte(1)*Pfd/Hfd + cte(2);
% 
% cte(7)=log(1/cte(6) - 1);
% cte(8)=log(1/cte(5) - 1);
% kdl=(cte(5)-cte(6))/Pfd
% cdl=cte(5)-kdl*Pfd
 

