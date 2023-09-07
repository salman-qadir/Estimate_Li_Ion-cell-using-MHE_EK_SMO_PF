a1=load('C:\Users\salma\MATLAB Drive\T11 MHE\1csmos\dsfhm11t1zsscec.mat');
a2=load('C:\Users\salma\MATLAB Drive\T11 MHE\1csmos\dsfhm11t1zd2c.mat');
dsfhm11t1zssce=a1.dsfhm11t1zsscec;
dsfhm11t1zd2=a2.dsfhm11t1zd2c;
% initial11t1c;
% load('C:\Users\salma\MATLAB Drive\T11 MHE\onecdata')
% g1=lpf(g1,7e-1);g2=lpf(g2,7e-1);%g1=[1,g1];g2=[g2(1),g2];
% g3=1:g1(end); g4=interp1(g1,g2,g3); %%
% p.gt=g3;p.gv=g4+ra(1:length(g4));
% tr=min([length(dsfhm11t1cssce.y),length(g3)]);
% figure
% yyaxis left
% plot(dsfhm11t1cssce.y,'--')
% %t=1:length(dsfhm11t1cd2.v);
% hold on
% plot(p.gv(1:tr),'-.')
% ylabel('Voltage[V]')
% yyaxis right
% error=(p.gv(1:tr)-dsfhm11t1cssce.y(1:tr));
% rmsv=(sqrt(sum(error.^2)/length(error)))*(100/mean(dsfhm11t1cssce.y(1:tr)))
% 
% plot(error(1,4:end),':')
% legend('Estimate','Experiment','Error')
% xlabel('Time[s]')
% ylabel('Error[V]')
%%
figure
yyaxis left
plot(dsfhm11t1zssce.y,'-.')
%t=1:length(dsfhm11t1cd2.v);
hold on
plot(dsfhm11t1zd2.v(1:length(dsfhm11t1zd2.v)-1))
ylabel('Voltage[V]')
yyaxis right
error=(dsfhm11t1zd2.v(1:length(dsfhm11t1zd2.v)-1)-...
    dsfhm11t1zssce.y(1:length(dsfhm11t1zd2.v)-1));
rmsv=(sqrt(sum(error.^2)/length(error)))*(100/mean(dsfhm11t1zssce.y))

plot(error(1,4:end),'--')
 legend('Estimate','Plant','Error')
 xlabel('Time[s]')
 ylabel('Error[V]')
%  ylim([-6e-3,6e-3])
% create a new pair of axes inside current figure
% axes('position',[.65 .195 .25 .25])
% box on % put box around new pair of axes
% it = (t < 430)&(t >400); % range of t near perturbation
% plot(t(it),dsfhm11t1d2.v(it)) % plot on new axes
%  hold on
%  plot(t(it),dsfhm11t1ssce.y(it)) % plot on new axes
% %axis tight
% xlim([400,430])
% ylim([3.5,3.85])
%%
figure
yyaxis left
plot(dsfhm11t1zssce.socp,'--')
t=1:length(dsfhm11t1zd2.v);
hold on
plot(dsfhm11t1zd2.socp(:,p.N:end),'--')
hold on
plot(p.socr,'-')
%legend('Estimate','Plant','Reference')
xlabel('Time[s]')
ylabel('SoC[%]')
yyaxis right%length(dsfhm11t1bd2.socp)-1
error1=(dsfhm11t1zd2.socp-...
    dsfhm11t1zssce.socp);
rmsc=(sqrt(sum(error1.^2)/length(error1)))*(100/mean(dsfhm11t1zd2.socp(:,3:end)))

% plot(error1(1,p.N:end),'-.')
% legend('Estimate','Plant','Reference','Error')
% xlabel('Time[s]')
% ylabel('Error[%]')
%%
 figure
 plot(abs(dsfhm11t1zssce.u(1,3:end)*p.a),'-.')
% %legend('Estimate','Plant')
 xlabel('Time[s]')
 ylabel('Input current [A]')
% a1=load('C:\Users\salma\MA
% %%TLAB Drive\T10_SEI\sfhmmpc\dsfhm10t1ss.mat');
% a2=load('C:\Users\salma\MATLAB Drive\T10_SEI\sfhmmpc\dsfhm10t1sscc.mat');
% sqq=load('C:\Users\salma\MATLAB Drive\T9_ss_mpc\New Folder\sqq.mat');
% dsfhm10t1ss=a1.dsfhm10t1ss;
% dsfhm10t1sscc=a2.dsfhm10t1sscc;
% sq=sqq.sqq;
% initial10t1
% t1=linspace(1,dsfhm10t1ss.ij*p.t1,length(dsfhm10t1ss.v));
% t2=linspace(1,dsfhm10t1sscc.ij*p.t1,length(dsfhm10t1sscc.v));
% t3=linspace( 1,length(p.socr)*p.t1,length(p.socr));
% %subplot(5,1,1)
% %plot(sq.v(:,3),'-.','LineWidth',2);hold on;
% % figure
% % plot(t1,dsfhm10t1ss.v,'-.','LineWidth',2);hold on;
% % plot(t2,dsfhm10t1sscc.v,'--','LineWidth',2);hold on;
% % legend('proposed controller','CCCV')
% % xlabel('Time[s]')
% % ylabel('Voltage[v]')
% figure
% plot(t1,dsfhm10t1ss.x(:,6),'-.','LineWidth',2);hold on;
% plot(t2,dsfhm10t1sscc.x(:,6),'--','LineWidth',2);hold on;
% xlabel('Time[s]')
% ylabel('R film[\ohm m^{2}]')
% legend('proposed controller','CCCV')
% % figure
% % plot(t1,abs(p.a*dsfhm10t1ss.u),'-.','LineWidth',2);hold on;
% % plot(t2,abs(p.a*dsfhm10t1sscc.u),'--','LineWidth',2);hold on;
% % xlabel('Time[s]')
% % ylabel('Input current[A]')
% % legend('proposed controller','CCCV')
% %  figure
% %  plot(t1,dsfhm10t1ss.opns(1:end),'-.','LineWidth',2);hold on;
% %  plot(t2,dsfhm10t1sscc.opns (1:end),'--','LineWidth',2);hold on;
% % xlabel('Time[s]')
% % ylabel('Over potential[V]')
% % legend('proposed controller','CCCV')
% figure
% plot(t1,abs(dsfhm10t1ss.socn),'-.','LineWidth',2);hold on;
% plot(t2,abs(dsfhm10t1sscc.socn),'--','LineWidth',2);hold on;
% plot(t2,99.5*ones(1,length(t2)),'-','LineWidth',2);hold on;
% xlabel('Time')
% ylabel('SoC[%]')
% legend('proposed controller','CCCV','Reference')
% figure
% plot(t1,100*(dsfhm10t1ss.soh),'-.','LineWidth',2);hold on;
% %plot(abs(p.soh),'-','LineWidth',2);hold on;
% plot(t2,100*(dsfhm10t1sscc.soh),'--','LineWidth',2);hold on;
% xlabel('Time[s]')
% ylabel('SoH[%]')
% % legend('proposed controller','CCCV')
% % figure
% % plot(dsfhm10t1ss.x(:,6),'-.','LineWidth',2);hold on;
% % plot(dsfhm10t1ss1.x(:,6),'-.','LineWidth',2);hold on;
% % legend('rfilm')
% % figure
% % plot(abs(dsfhm10t1ss.u),'-.','LineWidth',2);hold on;
% % plot(abs(dsfhm10t1ss1.u),'-.','LineWidth',2);hold on;
% % legend('Input 1', 'Input 2')
% % xlabel('time')
% % ylabel('Input Current[A]')
% % figure
% % plot(dsfhm10t1ss.opns(3:end),'-.','LineWidth',2);hold on;
% % plot(dsfhm10t1ss1.opns(3:end),'-.','LineWidth',2);hold on;
% % % plot(dsfhm10t1ss.opn,'-.','LineWidth',2);hold on;
% % legend('\eta_{nsr,1}','\eta_{nsr,2}')
% % xlabel('Time')
% % ylabel('\eta[Volt]')
% % figure
% % plot(abs(dsfhm10t1ss.socn),'-.','LineWidth',2);hold on;
% % plot(abs(dsfhm10t1ss1.socn),'-.','LineWidth',2);hold on;
% % legend('SoC 1', 'SoC 2')
% % xlabel('Time')
% % ylabel('SoC')
