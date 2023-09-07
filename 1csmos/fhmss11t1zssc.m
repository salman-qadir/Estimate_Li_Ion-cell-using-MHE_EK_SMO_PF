function[y] = fhmss11t1zssc(u,d,p)
 ce=d.cel(:,p.ij-p.N+1:p.ij);
 p.cen=ce(1:p.n,:);
 p.cep=ce(p.n+p.p+1:p.x,:);
 y.cel=ce;y.qf=d.qf;
%%
t=p.tc;
u=reshape(u,length(u),1);
[a,b]=size(d.x1);
xa=reshape(d.x1,a,b);
cpp=xa(1);
r0=d.rfilm(end);q0=d.soh(end)*p.c;
jsn=d.jsn1(end);
jn11=u/(p.ln);y.jp1=-u/(p.lp);
% xpp=cpp/p.csp;
% yp=mean(xpp);
% y.socp=(100*(yp-p.xp0)/(p.xp1-p.xp0));
% y.socn=y.socp;
% xnn=((y.socn/100)*(p.xn1-p.xn0)+p.xn0);
% cnn=xnn'*p.csn;
% [uns,~,~,~] = ocp11t1(cnn/p.csn,cpp/p.csn);
% ecdn=real(mean(p.kn*sqrt((mean(p.cen).*cnn).*(1-cnn/p.csn))));
% ajs=-(p.an*p.ios)*exp(-p.kb*(uns-p.uref-.05));
% cjs=(1)./(2*ecdn);
% bjs=(-(u)/(p.f*p.ln*2*ecdn));
% y.jsn=((bjs+sqrt( bjs.^2+1-2.*cjs.*ajs )  )./(1./ajs-2.*cjs));
y.jn1=jn11-fm(jsn)';
%% anode  
% c1nn1=cnn+(p.xm*t)*(p.f\(-y.jn1/(3*p.nsn) - (u.*p.mk)/p.ln));
% c1nn=c1nn1(end);
% b=-(jn1(:,end)*p.ln^2)/(p.f*p.dsn*6*p.nsn);
% a= c1nn-b/3;
% csn=(a+b*(p.zn./p.ln).^2);
% yn=mean(c1nn/p.csn);
% socn=(100*(yn-p.xn0)/(p.xn1-p.xn0));
% xpp=((socn/100)*(p.xp1-p.xp0)+p.xp0);
% c1pp=xpp'*p.csp;
% b1=-(jp1*p.lp^2)/(p.f*p.dsp*6*p.nsp);
% a1= c1pp-b1/3;
% csp=fliplr(a1+b1*((p.zp/p.lp).^2));
% yp=mean(fliplr(csp')/p.csp);
% socp=(100*(yp-p.xp0)/(p.xp1-p.xp0));
%% cathode
c1pp =cpp+(p.xxm.*t).*(-(3*y.jp1)/(p.rs*p.ap*p.f));
%c1pp=fm(cpp)'+(p.xm*t)*(p.f\(-y.jp1/(3*p.nsp) + (u.*p.mk)/p.lp));
%xpp=((socn/100)*(p.xp1-p.xp0)+p.xp0);
%c1pp=xpp'*p.csp;
% b1=-(jp1*p.lp^2)/(p.f*p.dsp*6*p.nsp);
% a1= c1pp-b1/3;
% csp=fliplr(a1+b1*((p.zp/p.lp).^2));
xpp=c1pp/p.csp;yp=xpp;
y.socp=(100*(yp-p.xp0)/(p.xp1-p.xp0));
y.socn=y.socp;
xnn=((y.socn/100)*(p.xn1-p.xn0)+p.xn0);
cnn=xnn'*p.csn;
% b=-(jn1*p.ln^2)/(p.f*p.dsn*6*p.nsn);
% a= c1nn-b/3;
% csn=(a+b*(p.zn./p.ln).^2);
% yn=mean(xnn);
% socn=(100*(yn-p.xn0)/(p.xn1-p.xn0));

%% SEI
[y.un,y.up,~,~] = ocp11t1z(xnn,xpp);%qb(1)/(p.nen*p.ln),qb(3)/(p.nep*p.lp)
ecdn=real((p.an*p.knp.*sqrt((mean(p.cen)'.*cnn').*(1-cnn'/p.csn))));
ecdp=real((p.ap*p.kpp.*sqrt((mean(p.cep)'.*cpp').*(1-cpp'/p.csp))));
%ecdn=real((p.kn.*sqrt((fm(mean(p.ce))'.*cnn').*(1-xnn))));
%ecdp=real((p.kp*sqrt((fm(mean(p.ce))'.*c1pp).*(1-xpp))) );
% ecdn1=real((p.kn.*sqrt(((mean(y.cen))'.*cnn').*(1-xnn))));
% ecdp1=real((p.kp*sqrt(((mean(y.cep))'.*cpp').*(1-xpp))) );
y.opn=p.kb\asinh(y.jn1./(2*ecdn)); 
y.opp=p.kb\asinh(y.jp1./(2*ecdp));
%%
opns1=y.opn+y.un-p.uref+0.36; y.opns1=opns1;
y.jsn1=-p.an*p.ios*exp(-p.kb*y.opns1); %%1e8
qtt=2.1e-1*(p.a*p.ln*y.jsn1);y.qtt=qtt;%3e-4
rfilmt=-2e-3*(p.mp*y.jsn1.*p.mk)/(p.pp*p.kps*p.f);
%%
%y.opns1=y.opn+y.un-p.uref+0; 
%y.jsn1=-p.an*p.ios*exp(-p.kb*y.opns1); 
%qtt=3e-4*(p.a*p.ln*y.jsn);

y.qt=fm(q0)'+t*p.xm*fm(qtt)';
%rfilmt=-4e-6*(p.mp*fm(y.jsn)'.*p.mk)/(p.pp*p.kps*p.f);%
y.rfilm=fm(r0)'+t*p.xm*fm(rfilmt)';
%y.x=[c1pp(end),y.rfilm(end),y.qt(end),y.jsn1(end)];
%y.x=[c1pp(end),y.rfilm(end),y.qt(end)];
%  phied=real(((p.ln+p.lp+2*p.ls)*u)/(2*p.ke) +...
%      p.kb\p.tp*p.ke*( log(y.cel(end,:))-log(y.cel(1,:))  )');
phied=real(((p.ln+p.lp+2*p.ls)*u)/(2*p.ke) +...
    p.kb\p.tp*p.ke*( log(p.ce(1)')-log(p.ce(end)')  )');
%y.opns=y.opn+y.un+1*u*p.rc*p.a;
y.soh=y.qt/p.c;
y.v=(y.opp-y.opn+phied+y.up-y.un-u*p.rc*p.a)';
y.y=y.v;
end
function o=fm(i)
o=reshape(i,1,length(i));
end
