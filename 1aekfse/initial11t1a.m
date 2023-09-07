addpath(genpath('C:\Users\salma\MATLAB Drive\T11 MHE\1aekfse') );
% clearvars; close all; clc
p=param11t1a;
ab=p.c;bc=p.rsei;soh=1;rs=p.rsei;
% p.c=ab; p.rsei=bc;
% p.cr=1.525;
p.cr=-1;p.tme= 0.0;
p.u1=(((p.cr*1*p.c))/p.a);
p.uu=p.u1*ones(1,4e3);%p.u=p.uu;
p.soci=0;p.socr=103;
% p.y0.soci=df.soci;p.y0.socf=df.socf;
% cnn=p.csn*p.xn1;cpp=p.csp*p.xp1;
xnn=((p.soci/100)*(p.xn1-p.xn0)+p.xn0);
xpp=((p.soci/100)*(p.xp1-p.xp0)+p.xp0);
cnn=xnn*p.csn;cpp=xpp*p.csp;
cn=cnn*ones(p.n,1);cp=cpp*ones(p.p,1);
p.y0.cel=p.ce*ones(p.x,1); p.cen=p.ce*ones(p.n,1);
p.vl=[2.5,4.7];p.snl=[1e-4,p.csn];p.spl=[1e-4,p.csp];
p.ecl=[1e-6, 5*p.ce];
 [phisn1,phisp1]=ocp11t1a(xnn,xpp);
 p.y0.v=phisp1-phisn1;p.y0.jsn=0;
% p.jn=(p.u/p.ln)*ones(length(p.zn),1);
% p.jp=-(p.u/p.lp)*ones(length(p.zp),1);
 p.ts=0:p.t1:abs(p.cr)\(3.780e3); %[0 3786];
% p.ts=0:p.t:(p.cr)\(1e4*abs((p.socf-p.soci)*.01));
%%
p.qb=[p.nen*p.ce*p.ln;p.nes*p.ce*p.ls;p.nep*p.ce*p.lp];
p.g=[p.ce,0,p.ce,0,0,p.ce,0]';
p.qbl=length(p.qb);p.gl=length(p.g);
% yd=[cn;cp;p.ce1];
% y1=[phisn1;phisp1;p.jn;p.jp];
% y0 = consist(y1,yd,p);
% p.y0=[p.qb;cnn;cpp;y0];
% v=phisp1-phisn1;
% p.y0=[p.qb;cnn;cpp;v;p.soci;p.soci;cn;cp;p.ce1];
p.y0.socp=p.soci;p.y0.opns=phisn1;p.y0.un=phisn1;
p.y0.opn=0;p.y0.opns1=phisn1;p.y0.jsn=0;p.y0.jsn1=0;p.y0.up=phisp1;
%p.y0.x=[p.qb;cpp;bc;ab;p.y0.jsn1]';p.y0.y=phisp1-phisn1;
p.y0.x=[cpp;bc;ab;p.y0.jsn1]';p.y0.y=phisp1-phisn1;p.nx=length(p.y0.x);
p.y0.socn=p.soci;p.y0.soh=ab/p.c;p.y0.rfilm=p.rsei;
p.y0.u0= p.u1*ones(1,p.M);
p.socr = [linspace(p.soci,p.socr,round(length(p.ts))),p.socr*ones(1,2e3)];
p.soh=[linspace(1,0.99998,round(length(p.ts))),0.998*ones(1,3e1)];
p.iter=length(p.socr)-1;
%%
p.lb1=[p.xp1*p.csp,p.rsei,0.7*p.c,-1e6];
p.ub1=[p.xp0*p.csp,p.rsei*1e2,p.c,1e6];
p.lb=repmat(p.lb1,1,p.N);
p.ub=repmat(p.ub1,1,p.N);
p.no=length(p.y0.x);
%% ekf
p.smd.A=eye(4);
p.hx=@(xnn,xpp) ocp11t1aa(xnn,xpp);p.smd.D=0;
p.xe0=p.y0.x; 
p.bn=4;
p.dn=1;
p.xp(:,1)=p.xe0;
p.pq=0;
load('C:\Users\salma\MATLAB Drive\T11 MHE\1aekfse\noise.mat');
p.ra=ra;p.qa=qa;p.r1=r1;p.q1=q1;
function o=fm(i)
o=reshape(i,1,length(i));
end