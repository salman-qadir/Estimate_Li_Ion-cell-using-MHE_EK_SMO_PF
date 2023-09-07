p=struct;
p=initial(p);
%p=initialk(p);
p=initialek(p);
p=simt(p);
plots(p);
function p =simt(p)
    for i=1:(p.n-1)
       p.i=i;
       p=nsis(p);
       p=ekalm(p);
    end
end
%%
function [p] = initial(p)
p.ts=.1;
p.n=1000;p.t=0:p.ts:(p.n-1)*p.ts;
p.q1=1e-1;p.r1=1e-1;
end
function [p] = initialek(p)
p.cb=0.5*0;p.cc=0.4;p.ca=1;
p.f=@(x,p)[x(2);-p.ca*sin(x(1))-p.cb*x(2)+p.cc*p.ua];
p.fx=@(x,p) [1,0;-cos(x(1)),-p.cb];
p.fu=@(p) [0;p.cc];p.smd.C=[1,0];p.smd.D=0;
p.u=0*ones(1,length(p.t));
p.x0=[pi/2;1];
p.xe0=[pi/2-.1;1-.1];
p.ta=0;p.x=p.x0;
p.smd.A=p.fx(p.xe0,p); p.smd.B=p.fu(p);
p.bn=length(p.smd.B);p.smd.D=0;
p.dn=length(p.smd.D);
p.pq=zeros(length(p.smd.B));
p.xp=zeros(length(p.smd.B),p.n);
p.xp(:,1)=p.x0;
% plot(y)1
% p.sm=ss(p.A,p.B,p.C,p.D);
% p.smd=c2d(p.sm,ts);
% p.rs=lsim(p.smd,p.u,p.t,p.x(:,1));
end
function [p] = initialk(p)
p.A=[-.1,0;.7,-.8];p.B=[.5;.7];
p.C=[.1,.1];p.D=0;
p.u=ones(1,p.n);
p.sm=ss(p.A,p.B,p.C,p.D);
p.smd=c2d(p.sm,p.ts);
p.x(:,1)=[10;10];
p.rs=lsim(p.smd,p.u,p.t,p.x(:,1));
p.pq=zeros(length(p.B));
p.xp=zeros(length(p.B),p.n);
p.bn=length(p.smd.B);
p.dn=length(p.smd.D);
end
%%
function [p] =sis(p)
x=p.smd.A*p.x(:,p.i)+p.smd.B*p.u(:,p.i)+...
    p.q1*randn(p.bn,1);
p.x(:,p.i+1)=x;
y=p.smd.C*p.x(:,p.i)+p.smd.D*p.u(:,p.i)+...
    p.r1*randn(p.dn,1);
p.y(:,p.i)=y;
end
function [p] =nsis(p)
p.ua=p.u(p.i);
[ta,ya]=ode15s(@(ta,ya,u)p.f(ya,p),...
    [p.t(p.i) p.t(p.i+1)], p.x0 );
p.ta=[p.ta,ta(end)]; 
p.x=[p.x,fm(ya(end,:))'+p.q1*randn(p.bn,1)];
p.y(:,p.i)=p.smd.C*fm(ya(end,:))'+p.r1*randn(p.dn,1);
p.x0=ya(end,:);
end
function [p] = kalm(p)
xp=p.smd.A*p.x(:,p.i)+p.smd.B*p.u(:,p.i);
% predicted state KAlMAN
pp=p.smd.A*p.pq*p.smd.A'+p.q1*eye(p.bn);
% predicted error covariance
yp=p.smd.C*p.x(:,p.i)+p.smd.D*p.u(:,p.i);
% predicted output
ye=p.y(:,p.i)-yp;% predicted output error
s=p.smd.C*pp*p.smd.C'+p.r1*eye(p.dn);
%correction term
kg=(pp*p.smd.C')/s; %kalman gain
p.xp(:,p.i+1)=fm(xp)'+kg*ye; %updated state
p.pq=pp-kg*p.smd.C*pp'; %updated covariance
p.yp(:,p.i)=p.smd.C*p.xp(:,p.i);  %updated output
p.ye(:,p.i)=p.y(:,p.i)-p.yp(:,p.i); %updated output error
end
function [p] = ekalm(p)
% predicted state extended Kalman
p.ua=p.u(p.i);
[~,ya]=ode15s(@(ta,ya) p.f(ya,p), [p.t(p.i) p.t(p.i+1)], p.xe0 );
xp=ya(end,:);
p.smd.A=p.fx(p.xe0,p); p.smd.B=p.fu(p);
pp=p.smd.A*p.pq*p.smd.A'+p.q1*eye(p.bn);
% predicted error covariance
yp=p.smd.C*fm(xp)';
% predicted output
ye=p.y(:,p.i)-yp;% predicted output error
s=p.smd.C*pp*p.smd.C'+p.r1*eye(p.dn);
%correction term
kg=(pp*p.smd.C')/s; %kalman gain
p.xp(:,p.i)=fm(xp)'+kg*ye; %updated state
p.pq=pp-kg*p.smd.C*pp'; %updated covariance
p.yp(:,p.i)=p.smd.C*p.xp(:,p.i);  %updated output
p.ye(:,p.i)=p.y(:,p.i)-p.yp(:,p.i); %updated output error
p.xe0=p.xp(:,p.i);
end
function plots(p)
plot(p.y','-.')
hold on
% plot(p.ya,'--')
% hold on
plot(p.yp,'-')
legend('pos','speed','estimate')
end
function y= fm(x)
y=reshape(x,1,length(x));
end