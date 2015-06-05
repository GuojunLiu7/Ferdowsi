function out=F1_Sim(nargin)

global param;

Vec=nargin;
%**************************************************************************
% Parameters
%**************************************************************************
p.C_a       = Vec(1);
p.R0_c      = Vec(2);
p.DeltaV    = Vec(3);
p.tau       = Vec(4);
p.V_H       = Vec(5);
p.Betta_H   = Vec(6);
p.Alpha     = Vec(7);
p.Gamma     = Vec(8);
p.Delta_h   = Vec(9);
p.IHR       = Vec(10);
p.sig_sp    = Vec(11);
p.sig_Alpha = Vec(12);

%**************************************************************************
% Initial condition
%**************************************************************************
p.P_init    = Vec(13);
p.HR_init   = Vec(14);

%**************************************************************************
% differential equation
%**************************************************************************
lags=[p.tau];
tspan=0:30;

sol = dde23(@dde_pedde,lags,[p.P_init;p.HR_init],tspan);
tint = linspace(0,30,30);
yint = deval(sol,tint);

% yint=sol.y;
% tint=15*sol.x/max(sol.x);

param.BP_Sim=yint(1,:);
param.HR_Sim=yint(2,:);

out=[{param.BP_Sim}, {param.HR_Sim}];


    function dydt = dde_pedde(t,y,Z)
        ylag1 = Z(:,1);
        p.Rc  = p.R0_c*(1+p.Alpha*(1-sig(ylag1(1))));
        p.Ts  = 1 - sig(ylag1(1));
        p.Tp  = sig(y(1));
        
        %% First diff eqn
        dpadt = -y(1) / (p.Rc * p.C_a) + y(2) * p.DeltaV / p.C_a;
        %dpadt = -y(1) / (1+p.Alpha * sig(ylag1)) / p.eps_a  +  p.mu * y(2) / p.eps_a ;

        %% Second diff eqn
        dhdt = p.Betta_H * p.Ts / (1+p.Gamma*p.Tp) - p.V_H * p.Tp  + p.Delta_h*( p.IHR - y(2) );
        % dhdt =  p.betta * sig(ylag1) / ( 1 + p.Gamma * (1-sig(y(1))) ) / p.eps_h - ...
        %     p.nu * (1-sig(y(1))) / p.eps_h + p.Delta * (1-y(2)) / p.eps_h;

        %% state vector
        dydt = [dpadt;dhdt];
    end

    function out=sig(x)
        out=1./(1+exp(-p.sig_Alpha.*(x-p.sig_sp)));
    end
end