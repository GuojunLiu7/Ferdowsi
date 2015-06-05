%% Name: PEDRAM ATAEE             -            UBC Student Number: 32120073
%**************************************************************************
% Windowing - Input=(signal[n*1],window width,overlap) 
% Output=(windowed data, number of windows)
%**************************************************************************
function [win n]=F2_Window(sig,ww,ov)
x=size(sig,1);
win=[];

for i=1:floor((x-ww)/(ww-ov))+1
    q=sig((i-1)*(ww-ov)+1:(i-1)*(ww-ov)+ww,1);
    win=[win; q'];
end
win=win';
[m n]=size(win);
end