%% Name: PEDRAM ATAEE             -            UBC Student Number: 32120073
%**************************************************************************
%                         Fitness Function 
%**************************************************************************
function out=F1_Fitness(Vec)
global param;
param.t=param.t+1;
LoW=param.LoW;

y1 = F1_Sim([Vec,param.BP_Msrd_init,param.HR_Msrd_init]);
BP_Sim = cell2mat(y1(1));
HR_Sim = cell2mat(y1(2));

BP_Msrd = param.BP_Msrd';
HR_Msrd = param.HR_Msrd';

%**************************************************************************
% Fitness Function Type
%**************************************************************************
out= mean( 2*abs(BP_Msrd(1:end)-BP_Sim(30-LoW+1:30))./BP_Msrd(1:end)  ).^2  + mean(abs(HR_Msrd(1:end)-HR_Sim(30-LoW+1:30))./ HR_Msrd(1:end).^2 );

% s1=F2_DisInd(BP_Msrd(1:end),BP_Sim(30-LoW+1:30),'L2','LC2');
% s2=F2_DisInd(HR_Msrd(1:end),HR_Sim(30-LoW+1:30),'L2','LC2');
% out= s1 + s2;
end