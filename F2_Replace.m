%% Name: PEDRAM ATAEE             -            UBC Student Number: 32120073
%**************************************************************************
% Replace
%**************************************************************************
function output=F2_Replace(vargin)
y = cell2mat(vargin(1));
k = cell2mat(vargin(2));

yy=y;
yy(k)=[];
yymean = mean(yy);

y(k)=yymean;
output=y;
end