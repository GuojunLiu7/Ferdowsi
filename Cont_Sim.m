%clear all;
global param

%**************************************************************************
% Data Loading
%**************************************************************************
file_name=[{'486'},{'484'},{'477'},{'476'},{'474'},{'289'}]; %{'485'}{'456'}{'458'}
casenum=1;
plotATM1=['../../../[Education] PhD/Dataset/MIMIC/Data-MAT/',cell2mat(file_name(casenum)),'nm.mat'];
plotATM2=['../../../[Education] PhD/Dataset/MIMIC/Data-RAW/',cell2mat(file_name(casenum)),'nm.info'];
x=plotATM(plotATM1,plotATM2);

dl = min(3600,length(x));

BP = x(1:dl,1);
CO = x(1:dl,5);
HR = x(1:dl,8);


clear casenum; clear x; clear plotATM1; clear plotATM2; clear file_name;
%**************************************************************************
% Cleaning
%**************************************************************************
out=F2_Clean([{HR},{BP},{CO}]);
HR=cell2mat(out(1));
BP=cell2mat(out(2));
CO=cell2mat(out(3));

HR=HR./60;
SV =CO*1000./(HR.*60);
clear out;

%**************************************************************************
% Filtering
%**************************************************************************
c=5;
a = 1;
b = 1/c.*ones(1,c);
HR = [mean(HR)*ones(c-1,1);HR];
BP = [mean(BP)*ones(c-1,1);BP];
SV = [mean(SV)*ones(c-1,1);SV];

HR = filter(b,a,HR);
BP = filter(b,a,BP);
SV = filter(b,a,SV);

HR = HR(c:end);
BP = BP(c:end);
SV = SV(c:end);

dl_F = length(HR);
clear a; clear b; clear c;
%clear CO;

% % % %**************************************************************************
% % % % Windowing
% % % %**************************************************************************
L=5;
[HR_W N_W]=F2_Window(HR,L,0);
BP_W=F2_Window(BP,L,0);
SV_W=F2_Window(SV,L,0);

%% Initializing
param.LoW=L;
param.Fs=1;
param.BP_Msrd_init=BP(1);
param.HR_Msrd_init=HR(1);
bound.ul=1.5; bound.ll=0.5;
F1_Initialization({bound.ul,bound.ll});
rpt= 0.9 + (1.1-0.9).*rand;
param.nom  = [1.55,    0.6, 50,  3, 1.17*rpt,  0.84*rpt,  1.3*rpt,  0.2,   1.7,  1.66,  110*rpt, 0.05];
clear rpt;
%% Estimation
p(1,:)=param.nom;
Q=F1_Sim([p(1,:),param.BP_Msrd_init,param.HR_Msrd_init]);
x=cell2mat(Q(1));
y=cell2mat(Q(2));
x_sim(1:param.LoW,1) = x(30-param.LoW+1:30);
y_sim(1:param.LoW,1) = y(30-param.LoW+1:30);
clear Q;

for i=2:N_W
    p(i,3)=mean(SV_W(:,i));
    param.ub(3)=p(i,3);    param.lb(3)=p(i,3);
    out=F1_Runopt(BP_W(:,i),HR_W(:,i),p(i-1,:));
    p(i,:)=cell2mat(out(1));
    
    param.BP_Msrd_init=x(end);
    param.HR_Msrd_init=y(end);
    
    Q=F1_Sim([p(i,:),param.BP_Msrd_init,param.HR_Msrd_init]);
    
    x=cell2mat(Q(1));
    y=cell2mat(Q(2));
    clear Q;
    x_sim(param.LoW*(i-1)+1:param.LoW*i,1) = x(30-param.LoW+1:30);
    y_sim(param.LoW*(i-1)+1:param.LoW*i,1) = y(30-param.LoW+1:30);
end
%%
param.ub(3)=50*bound.ul;   param.lb(3)=50*bound.ll;
clear t; clear pp;
t=ones(1,N_W*param.LoW);
pp=zeros(N_W*param.LoW,12);
for i=1:12
    for j=1:N_W
        for k=1:param.LoW
            pp(param.LoW*(j-1)+k,i)=p(j,i);
        end
    end
end


cond=2;
time=1:dl_F;
time=time';




HeartRate.signals.values=[HR,HR];
HeartRate.signals.dimensions=2;
HeartRate.time=time;

BloodPressure.signals.values=[BP,BP];
BloodPressure.signals.dimensions=2;
BloodPressure.time=time;

LungVolume.signals.values=1*ones(dl,1);
LungVolume.time=time;

Ca.signals.values=1.55*ones(dl,1);
Ca.time=time;

Tau.signals.values=3*ones(dl,1);
Tau.time=time;
    
Ra.signals.values=0.6*ones(dl,1);
Ra.time=time;

AlphaSig.signals.values=0.05*ones(dl,1);
AlphaSig.time=time;

time_resp(1)=0;
y_resp(1)=1;
for i=2:length(AlphaSig.signals.values)
time_resp(i)=time_resp(i-1)+(5*rand);
y_resp(i)=1+4*rand;
end

if cond==2
    StrokeVolume.signals.values=[SV,SV];
    StrokeVolume.signals.dimensions=2;
    StrokeVolume.time=time;
    
    CardiacVagal.signals.values=pp(:,5);
    CardiacVagal.time=time;
    
    CardiacSymp.signals.values=pp(:,6);
    CardiacSymp.time=time;
    
    Alpha.signals.values=pp(:,7);
    Alpha.time=time;
    
    P0.signals.values=pp(:,11);
    P0.time=time;
end

if cond==1
    dl=length(AlphaSig.signals.values);
    dl_F=length(AlphaSig.signals.values);
    
    time=1:dl_F;
    time=time';
    
    StrokeVolume.signals.values=[50*ones(dl,1),50*ones(dl,1)];
    StrokeVolume.signals.dimensions=2;
    StrokeVolume.time=time;
    
    Alpha.signals.values=1.3*ones(dl,1);
    Alpha.time=time;
    
    P0.signals.values=100*ones(dl,1);
    P0.time=time;
    
    CardiacSymp.signals.values=0.84*ones(dl,1);
    CardiacSymp.time=time;
    
    CardiacVagal.signals.values=1.17*ones(dl,1);
    CardiacVagal.time=time;
end




simplot(ScopeData)



tarseem=0;
if tarseem==1
    
simplot(ScopeData)

%**************************************************************************
%HR and BP - Sim vs. Msrd
figure
subplot(2,1,2)
plot(ScopeData.time,ScopeData.signals.values(:,1),'b',ScopeData.time,ScopeData.signals.values(:,2),'r');
grid on
xlabel('Time [s]');
ylabel('HR [bps]');
legend('Measurement','APGC with a revised form')
axis tight
subplot(2,1,1)
plot(ScopeData1.time,ScopeData1.signals.values(:,1),'b',ScopeData1.time,ScopeData1.signals.values(:,2),'r');
grid on
xlabel('Time [s]');
ylabel('BP [mmHg]');
legend('Measurement','APGC with a revised form')
axis tight

%**************************************************************************


figure
subplot(2,1,2)
plot(ScopeData.time,ScopeData.signals.values(:,1),'b',atMBC,aMBCSigHR,'r',atAPGA,aAPGCSigHR,'g',atAPGArev,aAPGCrevSigHR,'k');
grid on
xlabel('Time [s]');
ylabel('HR [bps]');
legend('Measurement','MBC','APGC','APGC with a revised form')
axis tight
subplot(2,1,1)
plot(ScopeData1.time,ScopeData1.signals.values(:,1),'b',atMBC,aMBCSig,'r',atAPGA,aAPGCSig,'g',atAPGArev,aAPGCrevSig,'k');
grid on
xlabel('Time [s]');
ylabel('BP [mmHg]');
legend('Measurement','MBC','APGC','APGC with a revised form')
axis tight




%Model-Based Controller
%Adaptive Proportional Gain Controller
%**************************************************************************
%HR and BP - Sim
figure
subplot(2,1,1)
plot(ScopeData.time,ScopeData.signals.values(:,1),'b');
grid on
xlabel('Time [s]');
ylabel('Heart Rate [bps]');
axis tight
subplot(2,1,2)
plot(ScopeData1.time,ScopeData1.signals.values(:,1),'b');
grid on
xlabel('Time [s]');
ylabel('Blood Pressure [mmHg]');
axis tight
%**************************************************************************
%P0
figure
plot(ScopeData2.time,ScopeData2.signals.values(:,2),'b',ScopeData2.time,ScopeData2.signals.values(:,1),'r');
grid on
xlabel('Time [s]');
ylabel('P_0 [mmHg]');
legend('Identification','APGC with a revised form')
axis tight


%P0
figure
plot(ScopeData2.time,ScopeData2.signals.values(:,2),'b',atMBC,aMBCP0,'r',atAPGA,aAPGCP0,'g',atAPGArev,aAPGCrevP0,'k');
grid on
xlabel('Time [s]');
ylabel('P_0 [mmHg]');
legend('Identification','MBC','APGC','APGC with a revised form')
axis tight


% figure
% plot(ScopeData8.time,ScopeData8.signals.values(:,1),'r',ScopeData8.time,ScopeData8.signals.values(:,2),'b',...
%     ScopeData8.time,ScopeData8.signals.values(:,3),'k');
% grid on
% xlabel('Time [s]');
% ylabel('P_0 [mmHg]');
% legend('Proportional Controller','Advanced Controller (HR:0.9, BP:0.1)','Brain')
% title('P_0: Control Variable');
% axis tight
%**************************************************************************

figure
subplot(2,2,1)
%y=CardiacVagal.signals.values;
y=0.91*ones(1,length(CardiacVagal.time));
x=CardiacVagal.time;
plot(x,y)
grid on
xlabel('Time [s]');
ylabel('V_H');
axis([min(x),max(x),0.95*min(y),1.05*max(y)])

subplot(2,2,2)
%y=CardiacSymp.signals.values;
y=0.72*ones(1,length(CardiacSymp.time));
x=CardiacSymp.time;
plot(x,y)
grid on
xlabel('Time [s]');
ylabel('\beta_H');
axis([min(x),max(x),0.95*min(y),1.05*max(y)])

subplot(2,2,3)
%y=Alpha.signals.values;
y=1.37*ones(1,length(Alpha.time));
x=Alpha.time;
plot(x,y)
grid on
xlabel('Time [s]');
ylabel('\alpha');
axis([min(x),max(x),0.95*min(y),1.05*max(y)])

subplot(2,2,4)
y=P0.signals.values;
x=P0.time;
plot(x,y)
grid on
xlabel('Time [s]');
ylabel('P_0');
axis([min(x),max(x),0.95*min(y),1.05*max(y)])



annotation('ellipse',[0.68  0.25  0.15  0.12],...
    'EdgeColor','red',...
    'LineWidth',3);
annotation('textarrow',[.5 .68],[.2 .3],...
    'TextEdgeColor','none',...
    'TextLineWidth',3,...
    'TextBackgroundColor',[1 .9 .8],...
    'FontSize',12,...
    'Color','red',...
    'LineWidth',3,...
    'String','I want to show this');

end