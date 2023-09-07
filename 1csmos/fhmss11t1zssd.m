function[d, yt] = fhmss11t1zssd(d,qf,p)
x0=d.xq;
da = fhmss11t1zsslyte(qf,p);u=p.ua;
p.cen=da.cel(1:p.n);
p.cep=da.cel(p.n+p.p+1:end);
yt.cel=da.cel(:);
yt.qf=fm(da.qf)';
%qb=p.x2(1:p.qbl);
%%
u=reshape(u,length(u),1);
y=reshape(x0,length(x0),1);
cpp=y(3);
%yt.rfilm=y(1+1);yt.qt=y(1+2);
xpp=cpp/p.csp;
yp=mean(xpp);
yt.socp=(100*(yp-p.xp0)/(p.xp1-p.xp0));
yt.socn=yt.socp;
jsn=d.jsn1(end);
xpp=cpp/p.csp;
yp=mean(xpp);
yt.socp=(100*(yp-p.xp0)/(p.xp1-p.xp0));
yt.socn=yt.socp;
xnn=((yt.socn/100)*(p.xn1-p.xn0)+p.xn0);
cnn=xnn'*p.csn;
% [uns,~,~,~] = ocp11t1b(cnn/p.csn,cpp/p.csn);
% ecdn=real(mean(p.kn*sqrt((mean(p.cen).*cnn).*(1-cnn/p.csn))));
% ajs=-(p.an*p.ios)*exp(-p.kb*(uns-p.uref-.05));
% cjs=(1)./(2*ecdn);
% bjs=(-(u)/(p.f*p.ln*2*ecdn));
% jsn=((bjs+sqrt( bjs.^2+1-2.*cjs.*ajs )  )./(1./ajs-2.*cjs));
  jn11=u/(p.ln);yt.jp1=-u/(p.lp);
yt.jn1=jn11-fm(jsn)';
%% electrolyte simplified
%  g=p.ca*([0;0;0;0;fm(qb)']);
%  yt.cen=(g(2)'.*((p.zn).^2)+g(1)'.*ones(length(p.zn),1)');
%  yt.ces=(g(5)'.*((p.zs).^2)+g(4)'.*p.zs+g(3)'.*ones(length(p.zs),1)');
%  yt.cep=fliplr(g(7)'.*((p.zp).^2)+g(6)'.*ones(length(p.zp),1)');
% %cel=[yt.cen';yt.ces';(yt.cep)'];
% yt.cen=p.ce*ones(length(p.zp),1)';
% yt.cep=p.ce*ones(length(p.zp),1)';
% yt.ces=p.ce*ones(length(p.zp),1)';
%  cel=[yt.cen';yt.ces';(yt.cep)'];
% yt.cel=cel;
%% Volts
[yt.un,yt.up,~,~] = ocp11t1z(xnn,xpp);
ecdn=real((p.an*p.knp.*sqrt((mean(p.cen)'.*cnn').*(1-cnn'/p.csn))));
ecdp=real((p.ap*p.kpp.*sqrt((mean(p.cep)'.*cpp').*(1-cpp'/p.csp))));
%ecdn=real((p.kn.*sqrt((mean(p.cen)'.*cnn').*(1-xnn))));
%ecdp=real((p.kp*sqrt((mean(p.cep)'.*cpp').*(1-xpp))) );
% ecdn=real((p.kn.*sqrt(((mean(yt.cen))'.*cnn').*(1-xnn))));
% ecdp=real((p.kp*sqrt(((mean(yt.cep))'.*cpp').*(1-xpp))) );
yt.opn=p.kb\asinh(yt.jn1./(2*ecdn')); 
yt.opp=p.kb\asinh(yt.jp1./(2*ecdp'));
% phied=real(((p.ln+p.lp+2*p.ls)*u)/(2*p.ke) +...
%     p.kb\p.tp*p.ke*( log(cel(end,:))-log(cel(1,:))  )');
phied=real(((p.ln+p.lp+2*p.ls)*u)/(2*p.ke) +...
    p.kb\p.tp*p.ke*( log(p.cen(1))-log(p.cep(end))  )');
%%
opns1=yt.opn+yt.un-p.uref+0.36; yt.opns1=opns1;
yt.jsn1=-p.an*p.ios*exp(-p.kb*yt.opns1); %%1e8
% %%1e8
qtt=2.1e-1*(p.a*p.ln*yt.jsn1);yt.qtt=qtt(end);%3e-4
rfilmt=-2e-3*(p.mp*yt.jsn1.*p.mk)/(p.pp*p.kps*p.f);
qt=(d.soh(end)*p.c)+p.tc*p.xxm.*fm(qtt)';
yt.qt=qt(end);
rfilm=(d.rfilm(end))+p.tc*p.xxm.*fm(rfilmt)';
yt.rfilm=rfilm(end);
yt.x=cpp;
d.x0=repmat(cpp,1,p.N);d.x1=cpp;
%%
%yt.jsn1=jsn;%-p.an*p.ios*exp(-p.kb*yt.opns1);
%yt.opns=yt.opn+yt.un+1*u*p.rc*p.a;
yt.soh=yt.qt/p.c;
yt.v=(yt.opp-yt.opn+phied+yt.up-yt.un-u*p.rc*p.a)';
yt.y=yt.v;
 end
function o=fm(i)
o=reshape(i,1,length(i));
end