%% For  sir khurram 

%% system matrices

A=[0,1;5,-4];

B=[0;1];

C=[1,4];

ALF=[];
AVG=[];

%%   simulation variables

ts=.001; % sampling time

fs=1/ts;   % sampling frequency

tt=2;   %  total time
var=tt/ts;

tv=0:ts:tt; % time variable

%% pre allocation of variables

x= zeros(length(B),length(tv));
xe=x;

y= zeros(1,length(tv));
ye=y;
y1=y;
er=y;
u=y;
lf1=y;lf=lf1;

%% initializations

x1(:,1)=[0;.1]; % initial condition of system
x2(:,1)=[0;.1]; % initial condition of system
x3(:,1)=[0;.1]; % initial condition of system
x4(:,1)=[0;.1]; % initial condition of system




xe1(:,1)=[0;0]; % initial condition of observer
xe2(:,1)=[0;0]; % initial condition of observer
xe3(:,1)=[0;0]; % initial condition of observer
xe4(:,1)=[0;0]; % initial condition of observer


pc=[-5,-7]; % desired eigen values of system

ph=[-10,-13];% desired eigen values of observer

k= place(A,B,pc);  % controller design

h= place(A',C',ph); % observer design

q=5*eye(2); % symmetric matrix

p = lyap(A',q);
% filename='RESULT.xlsx';
% sheet
% RATE1=xlsread(filename,sheet,xlRange);
% RATE2=xlsread(filename,sheet,xlRange);
% RATE3=xlsread(filename,sheet,xlRange);
% RATE4=xlsread(filename,sheet,xlRange);


for i=1:length(tv)-1 % time loop

    %% system

    u1(i)=-k*xe1(:,i);  % control law
%     u2(i)=-k*xe2(:,i);  % control law
%     u3(i)=-k*xe3(:,i);  % control law
%     u4(i)=-k*xe4(:,i);  % control law
%    
    x1(:,i+1)=x1(:,i)+(ts*A*x1(:,i))+(ts*B*u1(i));  % state equation
%     x2(:,i+1)=x2(:,i)+(ts*A*x2(:,i))+(ts*B*u2(i));  % state equation
%     x3(:,i+1)=x3(:,i)+(ts*A*x3(:,i))+(ts*B*u3(i));  % state equation
%     x4(:,i+1)=x4(:,i)+(ts*A*x4(:,i))+(ts*B*u4(i));  % state equation
%     

    sn=.001*randn(1); % sensor noise
% 
%     if(SAMPLE==0)
%         y(i)=C*OLD+sn;
%     else
%         y(i)=C*x(:,i)+sn; %output equation
%     end


    y1(i)=C*x1(:,i)+sn; % output equation
%     y2(i)=C*x2(:,i)+sn; % output equation
%     y3(i)=C*x3(:,i)+sn; % output equation
%     y4(i)=C*x4(:,i)+sn; % output equation
   



    %% observer

    ye1(i)=C*xe1(:,i); % estimate of  output
%     ye2(i)=C*xe2(:,i); % estimate of  output
%     ye3(i)=C*xe3(:,i); % estimate of  output
%     ye4(i)=C*xe4(:,i); % estimate of  output
    


    er1(i)=(y1(i)-ye1(i)); % estimation error
%     er2(i)=(y2(i)-ye2(i)); % estimation error
%     er3(i)=(y3(i)-ye3(i)); % estimation error
%     er4(i)=(y4(i)-ye4(i)); % estimation error
    

%     AVE=(er1(i)+er2(i)+er3(i)+er4(i))/4;
    AVE=er1(i)/1;
    AVG=[AVG AVE];

    ct1 =h'*er1(i); % correction term of observer
%     ct2 =h'*er2(i); % correction term of observer
%     ct3 =h'*er3(i); % correction term of observer
%     ct4 =h'*er4(i); % correction term of observer
    

    RATE1=1000;
%     ;    RATE2=1;    RATE3=1;    RATE4=1;
   
    SAMPLE1=RATE1*ts;
%     SAMPLE2=RATE2*ts;
%     SAMPLE3=RATE3*ts;
%     SAMPLE4=RATE4*ts;
   
    
    if(SAMPLE1<1)
        xe1(:,i+1)=xe1(:,i);    
    else
        xe1(:,i+1)= xe1(:,i)+ts*A*xe1(:,i)+ts*B*u1(i)+ts*ct1; % observer state 
    end
%     if(SAMPLE2<1)
%         xe2(:,1)=xe2(:,1);
%     else
%         xe2(:,i+1)= xe2(:,i)+ts*A*xe2(:,i)+ts*B*u2(i)+ts*ct2;  % observer state 
%     end
%     if(SAMPLE3<1)
%         xe3(:,1)=xe3(:,1);
%     else
%         xe3(:,i+1)= xe3(:,i)+ts*A*xe3(:,i)+ts*B*u3(i)+ts*ct3;  % observer state 
%     end
%     if(SAMPLE4<1)
%         xe4(:,1)=xe4(:,1);
%     else
%         xe4(:,i+1)= xe4(:,i)+ts*A*xe4(:,i)+ts*B*u4(i)+ts*ct4;  % observer state 
%     end
    

%     xe(:,i+1)= xe(:,i)+ts*A*xe(:,i)+ts*B*u(i)+ts*ct;  % observer state 

    lf11(i)= xe1(:,i)'*p*xe1(:,i); %lyapunov function
    lf1(i)=-xe1(:,i)'*q*xe1(:,i);

%     lf22(i)= xe2(:,i)'*p*xe2(:,i); %lyapunov function
%     lf2(i)=-xe2(:,i)'*q*xe2(:,i);
% 
%     lf33(i)= xe3(:,i)'*p*xe3(:,i); %lyapunov function
%     lf3(i)=-xe3(:,i)'*q*xe3(:,i);
% 
%     lf44(i)= xe4(:,i)'*p*xe4(:,i); %lyapunov function
%     lf4(i)=-xe4(:,i)'*q*xe4(:,i);

    

%     ALF=[ALF ((lf4(i)+lf3(i)+lf2(i)+lf1(i))/4)];
    ALF=[ALF (lf1(i)/1)];

end

%% output plot

figure;

plot(y1); hold on;

plot(ye1);

xlabel('Samples')

ylabel('Amplitude')

legend('Output','Estimate of Output')
xlim([0,var]);

%% states plot

figure

subplot(2,1,1)

plot(x1(1,:)); hold on

plot(xe1(1,:));

xlabel('Samples')

ylabel('Amplitude')

legend('State 1','Estimate of  State 1')
xlim([0,var]);

subplot(2,1,2)

plot(x1(2,:)); hold on

plot(xe1(2,:));

xlabel('Samples')

ylabel('Amplitude')

legend('State 2','Estimate of  State 2')
xlim([0,var]);

%%  lyapunov function plot

figure;

plot(lf11);

xlabel('Samples')

ylabel('Amplitude')

legend(' Lyapunov Function')
xlim([0,var]);

%% derivative lyapunov function plot

figure;

plot(lf1);

xlabel('Samples')

ylabel('Amplitude')

legend('Derivative of Lyapunov Function')
xlim([0,var]);

%%  average-derivative lyapunov function plot

figure;

plot(ALF);

xlabel('Samples')

ylabel('amplitude')

legend('Average of Lyapunov Function (All Plants)')
xlim([0,var]);
%%  average-ERROR plot

figure;

plot(AVG);

xlabel('Samples')

ylabel('Amplitude')

legend('Average Output Error (All Plants)')
xlim([0,var]);
