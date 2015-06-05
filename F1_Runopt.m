function out_main=F1_Runopt(BP,HR,x0)
global param
param.t=1;

param.BP_Msrd=BP;
param.HR_Msrd=HR;

options=optimset('MaxIter',10^2,'Algorithm','active-set','MaxFunEvals',10^5,'Display','off'); %'TolFun','TolCon',1e-4
x_est=fmincon(@F1_Fitness,x0,[],[],[],[],param.lb,param.ub,[],options);
y1=F1_Sim([x_est,param.BP_Msrd_init,param.HR_Msrd_init]);
BP_Sim= cell2mat(y1(1));
HR_Sim= cell2mat(y1(2));

out_main(1)={x_est};
out_main(2)={BP_Sim};
out_main(3)={HR_Sim};

%% Plot
% figure(2)
% hold on
% plot(BP_Sim,'r','LineWidth',2)
% hold on
% plot(param.BP_Msrd,'b','LineWidth',2)
% xlabel('time (s)','fontsize',11,'fontweight','b');
% ylabel('Blood Pressure (mmHg)','fontsize',11,'fontweight','b');
% legend('Estimated','Measured');
% grid on

% figure(3)
% hold on
% plot(HR_Sim,'r','LineWidth',2)
% hold on
% plot(param.HR_Msrd,'b','LineWidth',2)
% xlabel('time (s)','fontsize',11,'fontweight','b');
% ylabel('Hear Rate (beat/sec)','fontsize',11,'fontweight','b');
% legend('Estimated','Measured');
% grid on

% figure(4)
% bar(param.y)
% legend('Initial values','Estimated values','Nominal values');
% title('1: C_{a}, 2: R^{0}_{a}, 3: \DeltaV, 4:\tau, 5: V_H, 6: \beta_{H}, 7: \alpha, 8: \gamma, 9: \delta_{h}, 10: IHR, 11: sig_{sp}, 12: sig_{\alpha}');
%
% figure(5)
% for i=1:length(param.title)
% subplot(3,4,i)
% plot(param.y(i,:),'*-b','LineWidth',2)
% grid
% title(param.title(i));
% set(gca,'XTick',[1 2 3]);
% set(gca,'XTickLabel',{'Init','Est','Nom'})
% end

clear x0; clear x_ss; clear x;

    function out=gensig_est(X0,Y0)
        x1=F1_Sim([param.nom,X0,Y0]);

        BP_Gen= cell2mat(x1(1));
        HR_Gen= cell2mat(x1(2));

        %         figure(1)
        %         plot(BP_Gen,HR_Gen,'k','LineWidth',2);
        %         xlabel('Blood Pressure (mmHg)','fontsize',11,'fontweight','b');
        %         ylabel('Heart Rate (beat/sec)','fontsize',11,'fontweight','b');
        %         grid on

        out=[{BP_Gen},{HR_Gen}];
    end

end