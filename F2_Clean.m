%% Name: PEDRAM ATAEE             -            UBC Student Number: 32120073
%**************************************************************************
% Cleaning HR and BP simultaneously
%**************************************************************************
function out=F2_Clean(vargin)
HR=cell2mat(vargin(1));
BP=cell2mat(vargin(2));
if length(vargin)==3
    CO=cell2mat(vargin(3));
end
%**************************************************************************
% Not a Number (NaN) Removal + Physiology Constraints Removal
%**************************************************************************
k1=find(isnan(HR));
k2=find(isnan(BP));
k3=find(HR>200);
k4=find(HR<30);
k5=find(BP>220);
k6=find(BP<30);
k=[k1; k2; k3; k4; k5; k6];

if length(vargin)==2
    HR=F2_Replace([{HR},{k}]);
    BP=F2_Replace([{BP},{k}]);
end

if length(vargin)==3
    k7=find(isnan(CO));
k=[k1; k2; k3; k4; k5; k6; k7];
    HR=F2_Replace([{HR},{k}]);
    BP=F2_Replace([{BP},{k}]);
    CO=F2_Replace([{CO},{k}]);
end

if length(vargin)==2
    out={HR BP k};
end

if length(vargin)==3
    out={HR BP CO k};
end

end